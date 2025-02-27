# AlleleDoser

**Germline variant analysis with copy number and phased genotype integration**

## Description
This project provides an R script for analyzing germline variants in tumor samples by integrating germline variant counts, copy number data, metadata (e.g., tumor purity), and phased genotypes. The script calculates expected allele frequencies, confidence intervals, and final allele dosages, incorporating both copy number variation (CNV) and phased genotype information.

## Installation
Clone the repository and ensure you have the necessary R packages:
```bash
# Clone the repository
git clone https://github.com/HautaniemiLab/AlleleDoser.git
cd AlleleDoser

# Install required R packages (if not already installed)
Rscript -e "install.packages(c('tidyverse', 'Rsamtools'))"
```

## Usage
Run the script as follows:
```bash
Rscript alleleDoser.R \
  path/to/germline_variants.tsv \
  path/to/segment_file.tsv \
  path/to/metadata.tsv \
  path/to/output_directory \
  patientID \
  path/to/normal_samples.txt \
  path/to/phased_genotypes.bcf
```

## Inputs
The script requires **seven command-line arguments**:
1. Germline in tumor genotype
2. Segment file
3. Metadata file
4. Output directory
5. Patient ID
6. Normal sample IDs
7. Phased genotypes

### Input File Descriptions

#### 1. Germline in Tumor Genotype File (TSV)

You can generate the input TSV file from your VCF file using the following **GATK command**:

```bash
gatk VariantsToTable \
  -V input.vcf \
  -O output.tsv \
  -F CHROM -F POS -F REF -F ALT -GF AD
```
**Note**: The VCF file should contain a normal sample along with the tumor samples.

##### **Example Output Table Format:**

| CHROM | POS      | REF | ALT | patient1_sample1 | patient1_sample2 | patient1_sample3 |
|-------|----------|-----|-----|------------------|------------------|------------------|
| chr1  | 12345678 | A   | G   | 10,5             | 8,7              | 12,3             |
| chr2  | 87654321 | C   | T   | 20,2             | 15,4             | 18,6             |

---

#### 2. Segment File (TSV)

The segment file contains CNV information for each sample.

##### **Required Columns:**
- **`sample`**: Sample ID.  
- **`chr`**: Chromosome identifier.  
- **`startpos`**: Start position of the CNV segment (genomic coordinate).  
- **`endpos`**: End position of the CNV segment (genomic coordinate).  
- **`nMinor`**: Minor allele copy number within the segment.  
- **`nMajor`**: Major allele copy number within the segment.  

##### **Example Segment File Format:**

| sample            | chr     | startpos  | endpos    | nMinor | nMajor |
|-------------------|---------|-----------|-----------|--------|--------|
| patient1_sample1  | chr1    | 1000000   | 2000000   | 1      | 2      |
| patient1_sample1  | chr1    | 2000001   | 3000000   | 0      | 3      |
| patient1_sample2  | chr2    | 500000    | 1500000   | 1      | 1      |
| patient1_sample3  | chr2    | 1500001   | 2500000   | 2      | 2      |

---

#### 3. Metadata File (TSV)

The metadata file provides sample-specific tumor purity values.

##### **Required Columns:**
- **`sample`**: Sample ID.  
- **`purity`**: Tumor purity value (a numeric value between 0 and 1, representing the proportion of tumor cells in the sample).  

##### **Example Metadata File Format:**

| sample            | purity |
|-------------------|--------|
| patient1_sample1  | 0.75   |
| patient1_sample2  | 0.60   |
| patient1_sample3  | 0.85   |

---

#### 4. Output Directory

The output directory specifies the location where the final results (e.g., CSV files containing allele dosages and frequency calculations) will be saved.

- The directory **must already exist** before running the script.  
- The path should be provided as an absolute or relative path.  
- Ensure the directory has write permissions for saving output files.  

---

#### 5. Patient ID  

The **Patient ID** is used to filter sample-specific data from the input files and name the output file (`<patientID>.csv`).  

- The provided ID **must be included** in the sample names associated with that patient.  
- Sample names are assumed to follow a structure like `patient1_sample1`, `patient1_sample2`, `patient1_sample3`, etc., where `patient1` is the patient ID.

--- 

#### 6. Normal Samples File (TXT)

The normal samples file provides a list of normal (non-tumor) sample IDs.

- Plain text file (`.txt`), with one normal sample ID per line.  
- Sample IDs should match those used in the germline variant and segment files.

---

#### 7. Phased Genotype File (BCF)

The phased genotype file contains phased genotype information for each variant.

- Binary Call Format (`.bcf`) file, which should be indexed. 

## Output
The script generates a CSV file named `<patientID>.csv` in the specified output directory.

### Output Column Descriptions

| **Column**               | **Description**                                                                                        |
|--------------------------|--------------------------------------------------------------------------------------------------------|
| `sample`                 | Sample ID.                                                                                             |
| `CHROM`                  | Chromosome identifier.                                                                                 |
| `POS`                    | Genomic position of the variant.                                                                       |
| `REF`                    | Reference allele.                                                                                      |
| `ALT`                    | Alternate allele.                                                                                      |
| `Mutation`               | Combined variant notation (`CHROM:POSREF>ALT`).                                                        |
| `AD.0`                   | Reference allele read count.                                                                           |
| `AD.1`                   | Alternate allele read count.                                                                           |
| `DP`                     | Total read depth (`AD.0 + AD.1`).                                                                      |
| `AF`                     | Allele frequency of the alternate allele (`AD.1 / DP`).                                                |
| `purity`                 | Tumor purity from the metadata file.                                                                   |
| `nMajor`                 | Major allele copy number.                                                                              |
| `nMinor`                 | Minor allele copy number.                                                                              |
| `totalCN`                | Total copy number (`nMajor + nMinor`).                                                                 |
| `expMajorAF`             | Expected alternate allele frequency if the major allele is mutated.                                    |
| `expMajorCI.lo`          | Lower bound of the 95% confidence interval for `expMajorAF`.                                           |
| `expMajorCI.hi`          | Upper bound of the 95% confidence interval for `expMajorAF`.                                           |
| `expMajor.pbinom.extreme`| Probability of observing the alternate allele count or more extreme under the major allele hypothesis. |
| `expMinorAF`             | Expected alternate allele frequency if the minor allele is mutated.                                    |
| `expMinorCI.lo`          | Lower bound of the 95% confidence interval for `expMinorAF`.                                           |
| `expMinorCI.hi`          | Upper bound of the 95% confidence interval for `expMinorAF`.                                           |
| `expMinor.pbinom.extreme`| Probability of observing the alternate allele count or more extreme under the minor allele hypothesis. |
| `ALTallele`              | Preliminary classification of which allele (major or minor) carries the mutation.                      |
| `normal.AD.0`            | Reference allele read count in normal samples.                                                         |
| `GT`                     | Phased genotype.                                                                                       |
| `germline`               | Germline status classification (`0homo`, `1homo`, `hetero`).                                           |
| `majorfirst`             | Boolean indicating if the major allele is likely inherited first based on phased data.                 |
| `excl_mean`              | Excluded rolling mean of `majorfirst_numeric` within a window.                                         |
| `phasesegmented`         | Final phased segment classification (`TRUE`/`FALSE`) indicating phase-consistent segments.             |
| `ALTallelefinal`         | Final classification of the mutated allele after incorporating phased genotype information.            |
| `ALTdosagefinal`         | Final dosage estimate of the alternate allele accounting for copy number and phasing.                  |


## Contributors
- Yilin Li
- [Samuel Leppiniemi](https://github.com/SamuelLepp)

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

