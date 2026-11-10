"""Lossless composition of the audited S13 canvas; called by the module R script."""
from pathlib import Path
import csv
import hashlib
import io
import sys
import tempfile
import zipfile

import numpy as np
from PIL import Image, ImageDraw
from pypdf import PdfReader
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


def need(condition, message):
    if not condition:
        raise ValueError(message)


# The two audited exports differ in scale and palette. Coordinates are native
# pixels with half-open bounds. Each stripe is sampled at categorical cell centers.
PALETTE = [(119, 119, 119), (230, 159, 0), (204, 121, 167)]
GEOMETRY = {
    "embed": ((2284, 3599), (1692, 3134), (1597, 1632, 1667), (2190, 2225, 2260)),
    "fix": ((2919, 4230), (2074, 3394), (1940, 1990, 2040), (2785, 2835, 2885)),
}


def centers(bounds, n):
    left, right = bounds
    return (left + (np.arange(n) + 0.5) * (right - left) / n).astype(int)


def nearest_color(pixels, colors):
    delta = pixels.astype(float)[:, None, :] - np.asarray(colors)[None, :, :]
    return np.argmin(np.sum(delta * delta, axis=2), axis=1)


def check_tracks(pixels, tracks, kind):
    xb, yb, rows, cols = GEOMETRY[kind]
    ix, iy = centers(xb, len(tracks[0])), centers(yb, len(tracks[0]))
    for j, (name, expected) in enumerate(zip(("Detect", "PCA", "QC"), tracks)):
        colors = PALETTE if name == "QC" else PALETTE[:2]
        if kind == "fix" and name == "Detect":
            colors = [(255, 146, 137), (255, 132, 195)]
        for axis, sample in (("columns", pixels[rows[j], ix]), ("rows", pixels[iy, cols[j]])):
            observed = nearest_color(sample, colors)
            need(np.array_equal(observed, expected), f"{kind}: {name} differs on {axis}")


def detect_band(flags, width, height):
    supersampling = 12
    idx = ((np.arange(width * supersampling) + 0.5)
           * len(flags) / (width * supersampling)).astype(int)
    line = Image.fromarray(np.asarray(PALETTE[:2], dtype=np.uint8)[flags[idx]][None, :, :])
    line = line.resize((width, 1), Image.Resampling.LANCZOS)
    return line.resize((width, height), Image.Resampling.NEAREST)


def build(orig, panel_a, flags):
    need(panel_a.size == (1890, 1210), "Unexpected panel-A dimensions")
    candidate = orig.copy()
    panel_box = (8, 117, 1898, 1327)
    upper_band = (2284, 1583, 3599, 1613)
    left_band = (2175, 1692, 2205, 3134)
    legend_box = (1527, 2432, 1702, 2469)
    candidate.paste(panel_a, panel_box[:2])
    candidate.paste(detect_band(flags, 1315, 30), upper_band[:2])
    vertical = detect_band(flags, 1442, 30).transpose(Image.Transpose.TRANSPOSE)
    candidate.paste(vertical, left_band[:2])

    # Keep the source typeface: copy the 99.5% label, substitute its own '9'
    # for the '5', and draw the third (double-width, solid) key one row below.
    label = orig.crop((1601, 2373, 1699, 2410))
    label.paste(orig.crop((1603, 2373, 1622, 2410)), (46, 0))
    candidate.paste(label, (1601, 2432))
    ImageDraw.Draw(candidate).line((1532, 2449, 1577, 2449), fill="black", width=5)

    changed = np.any(np.asarray(orig) != np.asarray(candidate), axis=2)
    for left, top, right, bottom in (panel_box, upper_band, left_band, legend_box):
        changed[top:bottom, left:right] = False
    need(not np.any(changed), "Pixels changed outside the four declared edit areas")
    return candidate


def main():
    original, ref, panel, tracks, out = map(Path, sys.argv[1:6])
    available, cut = int(sys.argv[6]), float(sys.argv[7])
    if original.suffix.lower() == ".docx":
        with zipfile.ZipFile(original) as z:
            data = z.read("word/media/image14.png")
    else:
        data = original.read_bytes()
    need(hashlib.sha256(data).hexdigest() ==
         "8b6656daef49eeb13cfb2d71923f2ba209f8fab384eb745be2da41ab7ed5183b",
         "The original canvas differs from the reviewed S13 image.")
    need(hashlib.sha256(ref.read_bytes()).hexdigest() ==
         "9d4d25ab39100292275fc3def47a95d1b7258736ff3d89e5c8e977fdd459cb9a",
         "The reference differs from the audited corrected plotting export.")
    with Image.open(io.BytesIO(data)) as im:
        orig = im.convert("RGB")
    with tracks.open(encoding="ascii", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    old, new, pca, qc = [np.array([int(r[k]) for r in rows]) for k in
                         ("historical", "corrected", "pca", "qc")]
    n = len(new)
    need(0 < available <= n, "Invalid availability")
    check_tracks(np.asarray(orig), [old, pca, qc], "embed")
    with Image.open(ref) as im:
        check_tracks(np.asarray(im.convert("RGB")), [new, pca, qc], "fix")
    with Image.open(panel) as im:
        candidate = build(orig, im.convert("RGB"), new)
    check_tracks(np.asarray(candidate), [new, pca, qc], "embed")
    caption=f"""**Supplementary Figure 13. RNA-sequencing quality control.** Technical QC classes (PASS, grey; FLAG, orange; FAIL, pink) are based on uniquely mapped read percentage and sequencing depth. Final QC inclusion required technical PASS or FLAG and no more than one strong expression-level review flag. Panels B-D include all {n} samples before final QC selection, including samples subsequently excluded.

**A.** Preseq-estimated distinct reads (millions) at an evaluation depth of 442 million total reads, by technical QC class. Estimates were available for {available} of {n} samples ({n-available} unavailable). Only available preseq estimates are shown; no values were substituted or imputed.

**B.** Relative log expression (RLE) pooled across gene-sample observations within each technical QC class. Each RLE value is a retained gene's variance-stabilized expression in one sample minus that gene's median across samples. Boxes summarize the pooled values, rather than per-sample summaries. The dashed line marks zero.

**C.** PC1 and PC2 from the 5,000 most variable VST genes. Thin solid, dotted and thick solid lines mark the 97.5%, 99.5% and 99.9% two-dimensional reference contours, respectively. Larger points denote samples flagged by the separate ten-PC rule: squared robust Mahalanobis distance above the 99.5th percentile of the chi-squared distribution with ten degrees of freedom.

**D.** Hierarchically clustered sample-sample Spearman correlations on the same variable VST genes. Row and column annotations show technical QC (QC), the ten-PC review flag (PCA), and descriptive low detection (Detect). Detected genes are counted on the full, unrounded tximport lengthScaledTPM count matrix at >=10 before filtering or VST. Detect flags samples strictly below the cohort 10th percentile (<{cut:,.0f} genes; {new.sum()} of {n}). The distinct 0.5th-percentile strong flag used for main-analysis inclusion is not shown.

**Abbreviations:** QC, quality control; VST, variance-stabilizing transformation; HVGs, highly variable genes; PCA, principal component analysis; PC, principal component; MD, Mahalanobis distance; RLE, relative log expression.
"""
    # Publish only completed, verified outputs in the caller's empty directory.
    names = ("S13_corrected.png", "S13_corrected.pdf", "S13_caption.md")
    need(out.is_dir() and not any(out.iterdir()), "Output directory must be empty")
    with tempfile.TemporaryDirectory(dir=out) as stage:
        png, pdf, cap = [Path(stage) / name for name in names]
        candidate.save(png, "PNG", dpi=(300, 300), compress_level=8)
        dw, dh = [v * 72 / 300 for v in orig.size]
        c = canvas.Canvas(str(pdf), pagesize=(dw, dh), invariant=1, pageCompression=1)
        c.setTitle("Supplementary Figure 13. RNA-sequencing quality control")
        c.drawImage(ImageReader(candidate), 0, 0, width=dw, height=dh)
        c.showPage()
        c.save()
        pages = PdfReader(pdf).pages
        need(len(pages) == 1, "Expected one PDF page")
        images = [o.get_object() for o in pages[0]["/Resources"]["/XObject"].values()
                  if o.get_object().get("/Subtype") == "/Image"]
        need(len(images) == 1 and images[0].get_data() == candidate.tobytes(),
             "PDF pixels differ from the figure")
        cap.write_text(caption, encoding="ascii")
        for path in (png, pdf, cap):
            need(not (out / path.name).exists(), "Refusing to replace an output")
            path.rename(out / path.name)


if __name__ == "__main__":
    main()
