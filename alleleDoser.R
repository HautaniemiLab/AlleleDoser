#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(Rsamtools))

# Inputs as command line parameters
# 1. Germline in tumor genotype .csv file
# 2. Segment file
# 3. Meta-data file with e.g. tumor purities
# 4. Path to output directory
# 5. Patient ID
# 6. Normal sample IDs .txt file
# 7. Phased genotypes .bcf file

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=7) {
  stop("Provide seven arguments 
       \n germline variant count table 
       \n CN-data file 
       \n metadata file 
       \n path to output 
       \n patient_id 
       \n normal sample file
       \n phased genotype file
       \n", call.=FALSE)
}

germintumorfile <- args[1]
tumorcnfile     <- args[2]
metadatfile     <- args[3]
projectdir      <- args[4]
mypatient       <- args[5]
normalsfile     <- args[6]
phasefile       <- args[7]

# --------------------------------------------------
# Functions
# --------------------------------------------------

#Combine two lists into one by pasting element-wise
pasteByList = function(a_list, b_list, sep=":") {
  lapply(
    1:length(a_list),
    function(i) {
      paste(a_list[[i]], b_list[[i]], sep=sep)
    }
  )
}

#Function for computing the probability of drawing ALT read count or more extreme vs binomial distribution
pbinomExtreme = function(ad, dp, exp_af) {
  ifelse(
    ad / dp < exp_af,
    pbinom(ad, dp, exp_af),
    1 - pbinom(ad - 1, dp, exp_af)
  )
}

#Function for expected AF of germline heterozygous mutation given purity and tumour mutation and total copy number
expectedAF = function(mut_cn, total_cn, purity) {
  (mut_cn * purity + 1 * (1 - purity)) / (total_cn * purity + 2 * (1 - purity))   # NOTICE: One mutant copy present in normal!
}

rolling_sum <- function(x, window_size) {
  sapply(1:length(x), function(i) {
    start <- max(1, i - window_size)
    end <- min(length(x), i + window_size)
    sum(x[start:end], na.rm = TRUE)
  })
}

# --------------------------------------------------
# Inputs
# --------------------------------------------------

cnv_table = read.table(tumorcnfile, header=T, sep="\t", as.is=T) %>%
  filter(str_detect(sample, mypatient))

metadata_table = read.table(metadatfile, header=T, sep="\t", as.is=T) %>%
  filter(str_detect(sample, mypatient))

germline_in_tumour_table = read.table(germintumorfile, header = TRUE, sep = "\t", as.is = TRUE) 

germline_in_tumour_table = germline_in_tumour_table %>%
  mutate(samples = paste(gsub("\\.AD", "", colnames(.)[grepl("\\.AD", colnames(.))]),collapse = ";")) %>%
  unite("readCounts", ends_with("AD"), sep = ";") %>%
  select(-ends_with("AD"))

normal_samples <- readLines(normalsfile)

# --------------------------------------------------
# Processing with metadata and copy number data
# --------------------------------------------------

#Join the germline mutations
mutation_table = germline_in_tumour_table %>% 
  select(CHROM, POS, REF, ALT, readCounts, samples)

#Display mutations in desired format with AD mapped to AF and samples separated to rows
mutations_with_af_purity = mutation_table %>%
  mutate(Mutation = paste0(CHROM, ":", POS, REF, ">", ALT)) %>%
  mutate(AD = strsplit(readCounts, ";")) %>%
  mutate(sample = strsplit(samples, ";")) %>%
  mutate(sample.AD = pasteByList(sample, AD)) %>%
  unnest(sample.AD) %>%
  select(-sample, -AD, -readCounts, -samples) %>%
  separate(sample.AD, into=c("sample", "AD"), ":") %>%
  separate(AD, into=c("AD.0", "AD.1"), ",", convert=T) %>%
  mutate(DP = as.numeric(AD.0) + as.numeric(AD.1)) %>%
  mutate(AF = round(AD.1 / (AD.0 + AD.1), 2))

normals = mutations_with_af_purity[mutations_with_af_purity$sample %in% normal_samples, ]

#Merge purities and round numbers
mutations_with_af_purity = mutations_with_af_purity %>%
  merge(metadata_table %>% select(sample, purity) %>% add_row(sample="none", purity=0)) %>%
  mutate(purity = round(purity, 2),
         nMajor = NA, 
         nMinor = NA)

#Add copy number data
#Now in a for-loop, because looping over hundreds of segments is in any case faster than looping over millions of variants

for(i in 1:nrow(cnv_table)){
  variantindeces    <- with(mutations_with_af_purity, which(sample==cnv_table$sample[i] & 
                                                              CHROM==cnv_table$chr[i] & 
                                                              POS>=cnv_table$startpos[i] & 
                                                              POS<cnv_table$endpos[i]))
  if(length(variantindeces>0)){                                                          
    mutations_with_af_purity[variantindeces, c("nMajor", "nMinor")]  <-   cnv_table[i, c("nMajor", "nMinor")]
  }
  
  variantindeces    <- c()
}

mutations_with_cnv = mutations_with_af_purity %>%
  mutate(totalCN = nMajor + nMinor)

#Compute distributions for either cases of LOH
mutations_with_ci = mutations_with_cnv %>%
  mutate(expMajorAF = expectedAF(nMajor, totalCN, purity)) %>%
  mutate(expMajorCI.lo = qbinom(0.025, DP, expMajorAF)) %>%
  mutate(expMajorCI.hi = qbinom(0.975, DP, expMajorAF)) %>%
  mutate(expMajor.pbinom.extreme = pbinomExtreme(AD.1, DP, expMajorAF)) %>%
  mutate(expMinorAF = expectedAF(nMinor, totalCN, purity)) %>%
  mutate(expMinorCI.lo = qbinom(0.025, DP, expMinorAF)) %>%
  mutate(expMinorCI.hi = qbinom(0.975, DP, expMinorAF)) %>%
  mutate(expMinor.pbinom.extreme = pbinomExtreme(AD.1, DP, expMinorAF)) 

mutations_with_dosages <- mutations_with_ci %>%
  mutate(ALTallele = ifelse(nMajor == nMinor, "nMinor", 
                            ifelse(expMinor.pbinom.extreme > expMajor.pbinom.extreme, "nMinor", "nMajor")))

colnames(normals)[7] <- "normal.AD.0"

mutations_with_dosages <- merge(mutations_with_dosages, normals[, c("CHROM", "POS", "REF", "ALT", "normal.AD.0")], 
                                by = c("CHROM", "POS", "REF", "ALT"), 
                                all.x = TRUE)

mutations_with_dosages <- mutations_with_dosages %>%
  arrange(CHROM, POS, REF, ALT)

# --------------------------------------------------
# Processing with phased genotype data
# --------------------------------------------------

tempphase <- scanBcf(phasefile)

phase_table <- data.frame(
  CHROM = tempphase$CHROM,
  POS = tempphase$POS,
  REF = tempphase$REF,
  ALT = tempphase$ALT,
  GENO = tempphase$GENO
)

normals <- normals %>%
  mutate(flag=ifelse(AD.1/DP<0.2, "0homo", ifelse(AD.1/DP>0.8, "1homo", "hetero")))

mutations_with_dosages <- mutations_with_dosages %>%
  inner_join(phase_table, by = c("CHROM", "POS", "REF", "ALT")) %>%
  arrange(sample, CHROM, POS)

mutations_with_dosages <- mutations_with_dosages %>%
  merge(normals %>% select(CHROM, POS, REF, ALT, germline=flag), by=c("CHROM", "POS", "REF", "ALT"), all.x=T) %>%
  mutate(ALTallele = ifelse(germline=="0homo", "err", 
                            ifelse(germline=="1homo", "homoz", 
                                   ifelse(expMajor.pbinom.extreme>expMinor.pbinom.extreme, "nMajor",
                                          ifelse(expMajor.pbinom.extreme<expMinor.pbinom.extreme, "nMinor", "balanced"))))) %>%
  arrange(sample, CHROM, POS)

phasematrix <- mutations_with_dosages %>%
  filter(ALTallele %in% c("nMajor", "nMinor", "balanced")) %>% 
  select(sample, CHROM, POS, REF, ALT, nMajor, nMinor, GT, ALTallele) %>%
  arrange(sample, CHROM, POS) %>%
  mutate(majorfirst = ifelse(ALTallele=="balanced", NA, (GT=="0|1" & ALTallele=="nMinor") | (GT=="1|0" & ALTallele=="nMajor")))

phasematrix$majorfirst_numeric <- ifelse(is.na(phasematrix$majorfirst), 0, as.numeric(phasematrix$majorfirst))

phasematrix$phasesegmented <- NA
phasematrix$incl_mean <- NA
phasematrix$excl_mean <- NA

for(s in unique(phasematrix$sample)){
  
  tempmatrix    <- subset(phasematrix, sample==s)
  no_rows <- nrow(tempmatrix)
  
  tempmatrix$incl_sum <- rolling_sum(tempmatrix$majorfirst_numeric, 10)
  tempmatrix$incl_mean <- tempmatrix$incl_sum / 21
  
  tempmatrix$excl_sum <- tempmatrix$incl_sum - tempmatrix$majorfirst_numeric
  tempmatrix$excl_mean <- tempmatrix$excl_sum / 20
  
  tempmatrix$avgmajorfirst <- tempmatrix$incl_mean > 0.5
  tempmatrix$phasesegmented <- ifelse(is.na(tempmatrix$majorfirst) | 
                                        tempmatrix$excl_mean < 0.3 | 
                                        tempmatrix$excl_mean >= 0.7, 
                                      tempmatrix$avgmajorfirst, 
                                      tempmatrix$majorfirst)
  
  phasematrix$phasesegmented[which(phasematrix$sample==s)] <- tempmatrix$phasesegmented
  phasematrix$incl_mean[which(phasematrix$sample==s)] <- tempmatrix$incl_mean
  phasematrix$excl_mean[which(phasematrix$sample==s)] <- tempmatrix$excl_mean
  
}

# --------------------------------------------------
# Calculate final results and save them
# --------------------------------------------------

casefinal <- mutations_with_dosages %>%
  merge(phasematrix %>%
          select(sample, CHROM, POS, majorfirst, excl_mean, phasesegmented),
        by = c("sample", "CHROM", "POS"), all.x = TRUE) %>%
  mutate(ALTallelefinal = ifelse(is.na(ALTallele), NA, 
                                 ifelse(ALTallele %in% c("err", "homoz"), ALTallele,
                                        ifelse(phasesegmented & GT == "1|0", "nMajor",
                                               ifelse(phasesegmented & GT == "0|1", "nMinor",
                                                      ifelse(!phasesegmented & GT == "1|0", "nMinor", 
                                                             ifelse(!phasesegmented & GT == "0|1", "nMajor", NA))))))) %>%
  mutate(ALTdosagefinal = ifelse(is.na(ALTallelefinal), NA,
                                 ifelse(ALTallelefinal == "err", 0,
                                        ifelse(ALTallelefinal == "homoz", 1,
                                               ifelse(ALTallelefinal == "nMajor", nMajor/totalCN,
                                                      ifelse(ALTallelefinal == "nMinor", nMinor/totalCN, NA))))))

projectdir <- if (!grepl("/$", projectdir)) paste0(projectdir, "/") else projectdir
csv_name <- paste0(projectdir, mypatient, ".csv")
write_csv(casefinal, csv_name)

