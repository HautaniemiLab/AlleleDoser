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
git clone https://github.com/yourusername/germline-variant-analysis.git
cd germline-variant-analysis

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
The script generates a CSV file named `<patientID>.csv` in the specified output directory. The output includes:
- Chromosome and position information
- Reference and alternate alleles
- Calculated allele frequencies and dosages
- Copy number states and purity adjustments
- Phased genotype segmentation

## Contributors
- [Yilin Li]
- [Samuel Leppiniemi](https://github.com/SamuelLepp)

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

