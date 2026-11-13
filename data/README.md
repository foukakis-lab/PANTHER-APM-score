# Data sources and example inputs

This guide lists the input locations for the APM analysis. All paths below are relative to this `data/` directory, which belongs at the repository root.
Create the listed subdirectories locally and preserve the expected filenames, including capitalization and spaces.

Public datasets and annotation files are downloaded from their original sources. 
Real PANTHER data and the non-public TEX clinical and sample mapping files are not distributed with the code.

## Public downloads

### Annotation

| Save as | Source | Download instructions |
| --- | --- | --- |
| `annotation/gannot-cip-junecelltype.csv` | [JITC supplementary ZIP archive](https://jitc.bmj.com/content/jitc/5/1/18/media-5.default.zip) | Download the archive, extract the relevant CSV, and save it at the indicated path. |
| `annotation/LBL-10541-03-Breast-Cancer-360-Gene-List.xlsx` | [Bruker Breast Cancer 360 V2 gene list](https://brukerspatialbiology.com/resources/breast-cancer-360-v2-gene-list/) | Select **Download File** and save the Excel workbook at the indicated path. |

### MATADOR

| Save as | Source | Download instructions |
| --- | --- | --- |
| `MATADOR/20240507 CF32 MATADOR iscience mendeley data repository.csv` | [Mendeley Data, version 2](https://data.mendeley.com/datasets/n24jkp487r/2) | Download the named CSV from this version. |
| `MATADOR/20250401 CF32 MATADOR IHC manuscript Opdam 2025 clean survival data file.csv` | [Mendeley Data, version 3](https://data.mendeley.com/datasets/n24jkp487r/3) | Download the named CSV from this version. |
| `MATADOR/GSE167977_20210223_matadorCounts.csv.gz` | [GEO GSE167977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167977) | Download the named supplementary file and retain its `.gz` compression. |

Use the specified Mendeley versions for the two clinical files rather than substituting a later version.

### SCAN-B

Download the following files from [Mendeley Data, version 3](https://data.mendeley.com/datasets/yzxtxn4nmd/3).
Extract the files first if downloading them as an archive.

| Save as | Content |
| --- | --- |
| `SCAN-B/Gene.ID.ann.RData` | Gene annotation |
| `SCAN-B/SCANB.9206.genematrix_noNeg.RData` | Gene expression matrix |
| `SCAN-B/Supplementary_Data_Table_1_2023-01-13.xlsx` | Associated supplementary clinical and sample information (Original file name `Supplementary Data Table 1 - 2023-01-13.xlsx`) |

Preserve the R object names stored within the downloaded `.RData` files.

### TEX expression and probe annotation

Public microarray expression data are available from
[GEO GSE56493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56493).
The TEX analysis uses platform
[GPL10379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10379).

The main script downloads the GSE56493 series matrix and the GPL10379 annotation with `GEOquery` into TEX/GEO/ on the first run and reads them from there afterwards. 
The MCP-counter gene list is read from its GitHub repository on every run. 

The public GEO expression data do not replace the two non-public TEX inputs listed below.

## Synthetic PANTHER inputs

Synthetic versions of the two PANTHER inputs are included (made with `workflows/apm_analysis/scripts/PANTHER_make_synthetic_data.R`, seed 2027; see `README_SYNTHETIC_DATA.txt`):
| Path | Content |
| --- | --- |
| `PANTHER/PANTHER_clinical.rds` | 2,003 simulated randomized patients with the variables and coding used by the analysis; sample_id set for 444 |
| `PANTHER/PANTHER_expression_log2_TMM_TPM.rds` | Simulated log2(TPM + 1) matrix: 5,000 genes × the same 444 fictional samples |

## Non-public TEX inputs

The original versions of the following files cannot be distributed:

| Intended path | Role |
| --- | --- |
| `TEX/pd.unique.mets.RData` | Clinical input loaded as the R object `pd.unique.mets` |
| `TEX/TEX_clinical incl microarray SEB.xlsx` | Workbook used to match expression arrays to clinical sample records |

Synthetic versions of both files are included (made with workflows/apm_analysis/scripts/PANTHER_make_synthetic_data.R). 
They keep the object name, column names and workbook layout. 
Of the 140 synthetic records (SYN0001–SYN0140), 120 are linked at random to the public GSE56493 arrays; treatment labels, clinical characteristics and outcomes are simulated and say nothing about the GEO participants.

## Interpretation of synthetic examples

Synthetic examples demonstrate code execution and contain no real PANTHER or TEX clinical records. 
A treatment × (HLA6 − IFN5) interaction is deliberately introduced into the synthetic BCRFS data to exercise signature selection and downstream analyses. 
The resulting estimates, p-values and figures are illustrative and must not be interpreted as study results.

## Data handling and attribution

Keep downloaded datasets and any authorized real participant data outside version control. 
Only explicitly approved synthetic examples should be included with the code. 
Local data files may be present even when GitHub Desktop does not show them, because `.gitignore` excludes data inputs.

Each external dataset or annotation file remains subject to its source's terms of use and citation requirements. 
Refer to the linked source records for attribution. 
The repository's code license does not grant permission to redistribute external datasets or access controlled participant data.
