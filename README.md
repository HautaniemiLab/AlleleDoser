# AlleleDoser

**Germline variant analysis with copy number and phased genotype integration**

## Description
This project provides an R script for analyzing germline variants in tumor samples by integrating germline variant counts, copy number data, metadata (e.g., tumor purity), and phased genotypes. The script calculates expected allele frequencies, confidence intervals, and final allele dosages, incorporating both copy number variation (CNV) and phased genotype information.

## Inputs
The script requires **seven command-line arguments**:
1. **Germline in tumor genotype TSV file**: Variant counts with allelic depths.
2. **Segment file**: CNV data with major/minor allele copy numbers.
3. **Metadata file**: Contains tumor purity and other sample-specific information.
4. **Output directory**: Path where the results will be saved.
5. **Patient ID**: Identifier to filter data for the patient of interest.
6. **Normal sample IDs TXT file**: List of normal samples to reference allele frequencies.
7. **Phased genotypes BCF file**: Contains phased genotype information.

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
  path/to/segment_file.txt \
  path/to/metadata.txt \
  path/to/output_directory \
  patientID \
  path/to/normal_samples.txt \
  path/to/phased_genotypes.bcf
```

## Output
The script generates a CSV file named `<patientID>.csv` in the specified output directory.

### Output Column Descriptions

| **Column**               | **Description**                                                                 |
|--------------------------|---------------------------------------------------------------------------------|
| `sample`                 | Sample ID.                                                                     |
| `CHROM`                  | Chromosome identifier.                                                          |
| `POS`                    | Genomic position of the variant.                                                |
| `REF`                    | Reference allele.                                                               |
| `ALT`                    | Alternate allele.                                                               |
| `Mutation`               | Combined variant notation (`CHROM:POSREF>ALT`).                                 |
| `AD.0`                   | Reference allele read count.                                                     |
| `AD.1`                   | Alternate allele read count.                                                     |
| `DP`                     | Total read depth (`AD.0 + AD.1`).                                                |
| `AF`                     | Allele frequency of the alternate allele (`AD.1 / DP`).                         |
| `purity`                 | Tumor purity from the metadata file.                                             |
| `nMajor`                 | Major allele copy number.                                                        |
| `nMinor`                 | Minor allele copy number.                                                        |
| `totalCN`                | Total copy number (`nMajor + nMinor`).                                           |
| `expMajorAF`             | Expected alternate allele frequency if the major allele is mutated.              |
| `expMajorCI.lo`          | Lower bound of the 95% confidence interval for `expMajorAF`.                     |
| `expMajorCI.hi`          | Upper bound of the 95% confidence interval for `expMajorAF`.                     |
| `expMajor.pbinom.extreme`| Probability of observing the alternate allele count or more extreme under the major allele hypothesis. |
| `expMinorAF`             | Expected alternate allele frequency if the minor allele is mutated.              |
| `expMinorCI.lo`          | Lower bound of the 95% confidence interval for `expMinorAF`.                     |
| `expMinorCI.hi`          | Upper bound of the 95% confidence interval for `expMinorAF`.                     |
| `expMinor.pbinom.extreme`| Probability of observing the alternate allele count or more extreme under the minor allele hypothesis. |
| `ALTallele`              | Preliminary classification of which allele (major or minor) carries the mutation.|
| `normal.AD.0`            | Reference allele read count in normal samples.                                  |
| `GT`                     | Phased genotype.                                            |
| `germline`               | Germline status classification (`0homo`, `1homo`, `hetero`).                    |
| `majorfirst`             | Boolean indicating if the major allele is likely inherited first based on phased data. |
| `excl_mean`              | Excluded rolling mean of `majorfirst_numeric` within a window.                   |
| `phasesegmented`         | Final phased segment classification (`TRUE`/`FALSE`) indicating phase-consistent segments. |
| `ALTallelefinal`         | Final classification of the mutated allele after incorporating phased genotype information. |
| `ALTdosagefinal`         | Final dosage estimate of the alternate allele accounting for copy number and phasing. |


## Contributors
- [Yilin Li]
- [Samuel Leppiniemi](https://github.com/SamuelLepp)

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

