# Sample flow

The available expression and technical-QC inputs contain 506 unique RNA-seq samples and 506 unique nonmissing patient identifiers. This starting point is the sequenced and quantified population, not all trial participants.

STAR-based technical QC classifies 365 samples as PASS, 79 as FLAG, and 62 as FAIL. Technical QC produces:

- Strict technical set: 365 samples and 365 patients.
- Main-analysis technical set: 444 samples and 444 patients.

Expression-level QC then produces:

- Strict final set: 301 samples and 301 patients.
- Main-analysis final set: 444 samples and 444 patients.

The strict final set requires technical QC PASS and no expression-level review flags. The main-analysis final set requires technical QC PASS or FLAG and no more than one expression-level review flag.

The ordered main-analysis inclusion file is the sole sample-set input to the SSPBC intrinsic subtyping script. The script requires exactly 444 unique samples and patients and preserves inclusion-file order.
