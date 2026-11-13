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

Preserve the R object names stored within the downloaded `.Rdata` files.

### TEX expression and probe annotation

Public microarray expression data are available from
[GEO GSE56493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56493).
The TEX analysis uses platform
[GPL10379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10379).

The expression-loading section uses `GEOquery` to retrieve the GEO data.
Use the local download/cache location configured in the script; the initial download requires an internet connection. 

The public GEO expression data do not replace the two non-public TEX inputs listed below.

## Synthetic PANTHER inputs

The following example files are being prepared and are not yet included:

| Intended path | Planned content |
| --- | --- |
| `PANTHER/PANTHER_clinical.rds` | Synthetic clinical records with the variables and coding required by the analysis |
| `PANTHER/PANTHER_expression_log2_TMM_TPM.rds` | Synthetic gene expression matrix in the input format expected by the analysis |

The synthetic clinical records and expression columns must use matching fictional identifiers. 
These files are intended for the downstream APM analysis and do not supply the raw sequencing, QC, or subtyping inputs.

## Non-public TEX inputs

The original versions of the following files cannot be distributed:

| Intended path | Role |
| --- | --- |
| `TEX/pd.unique.mets.RData` | Clinical input loaded as the R object `pd.unique.mets` |
| `TEX/TEX_clinical incl microarray SEB.xlsx` | Workbook used to match expression arrays to clinical sample records |

Synthetic replacements are under consideration and are not yet included.
They must preserve the required object name, column names, workbook layout, and relationships between sample identifiers. 
The exact schema and matching rules will be checked against the final TEX analysis before preparing them.

If synthetic clinical records are paired with public GEO expression samples, the sample links will be explicitly artificial. 
Assigned treatment labels, clinical characteristics, and outcomes must not be interpreted as information about the original GEO participants.

## Interpretation of synthetic examples

Synthetic examples are intended only to demonstrate code execution. 
They will contain no real PANTHER or TEX clinical records and will not reproduce the study's estimates, p-values, signature-selection results, or figures.
They must be clearly identified as synthetic in their accompanying documentation and outputs.

## Data handling and attribution

Keep downloaded datasets and any authorized real participant data outside version control. 
Only explicitly approved synthetic examples should be included with the code. 
Local data files may be present even when GitHub Desktop does not show them, because `.gitignore` excludes data inputs.

Each external dataset or annotation file remains subject to its source's terms of use and citation requirements. 
Refer to the linked source records for attribution. 
The repository's code license does not grant permission to redistribute external datasets or access controlled participant data.
