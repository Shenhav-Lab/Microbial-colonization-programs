#!/bin/bash
################################################################################
# Summary of Nasal and Gut 16S rRNA gene data denoising and taxonomic assignment
################################################################################

# IMPORTANT NOTES:
## As-per recommendations for DADA2, this process was done on a per-run basis
## Nasal data was recieved intermittently as sequencing runs were completed over the course of ~1.5 years 
## DADA2 and initial QC checks were completed on a per-run basis intermittently as sequencing runs were being completed
## There were a total of 57 sequencing runs for nasal swab samples, and 5 runs for stool samples
### As an example, DADA2 and steps upstream of this are only shown here for one nasal sequencing run and all gut sequencing runs 
## activate QIIME2 (v2019.10, https://docs.qiime2.org/2019.1/) before running: source activate qiime2-2019.10
## cd 2-Datasets/Data/Preprocessing/raw_nasal_gut/

################################################################################
# Infant nasal swab preprocessing 
# Demultiplexed fastq files, V3 hypervariable region of the 16S rRNA gene
################################################################################

# Example per-run trimming, dada2 and QC check: Run LR389 (March 2021)
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path nasal_manifests/MANIFEST_LR389 \
  --output-path LR389-pa-demux.qza \
  --input-format PairedEndFastqManifestPhred33 && qiime cutadapt trim-paired \
  --i-demultiplexed-sequences LR389-pa-demux.qza \
  --p-front-f ATTACCGCGGCTGCTGG \
  --p-front-r CCTACGGGAGGCAGCAG \
  --p-discard-untrimmed \
  --o-trimmed-sequences 389-5trimmed.qza \
  --verbose \
&> 389_5trimming.log && qiime cutadapt trim-paired \
    --i-demultiplexed-sequences 389-5trimmed.qza \
    --p-adapter-f CTGCTGCCTCCCGTAGG \
    --p-adapter-r CCAGCAGCCGCGGTAAT \
    --p-discard-untrimmed \
    --o-trimmed-sequences 389-alltrimmed.qza \
    --verbose \
&> 389_alltrimming.log && qiime demux summarize --i-data 389-alltrimmed.qza --o-visualization 389-alltrimmed.qzv && qiime tools view 389-alltrimmed.qzv

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs 389-alltrimmed.qza \
    --p-trunc-len-f 130 --p-trunc-len-r 130 --p-trim-left-f 0 --p-trim-left-r 0 \
    --o-table table_389 --o-representative-sequences rep-seqs-389 --o-denoising-stats sample_data_389 --verbose \
&> 389_DADA2.log && qiime metadata tabulate --m-input-file sample_data_389.qza --o-visualization 389_denoising-stats.qzv && qiime tools view 389_denoising-stats.qzv


# Merge new runs with previously sequenced runs (already previously concatenated) - feature tables and representative sequences

qiime feature-table merge --i-tables new_nasal_table_NovSepJunAprFeb2020.qza \
    --i-tables table_375.qza --i-tables table_378.qza \
    --i-tables table_379.qza --i-tables table_380.qza \
    --i-tables table_382.qza --i-tables table_383.qza \
    --i-tables table_385.qza --i-tables table_386.qza \
    --i-tables table_387.qza --i-tables table_388.qza \
    --i-tables table_389.qza --o-merged-table new_nasal_table_March2021.qza  

qiime feature-table merge-seqs --i-data new_nasal_rep-seqs_NovSepJunAprFeb2020.qza \
    --i-data rep-seqs-375.qza --i-data rep-seqs-378.qza \
    --i-data rep-seqs-379.qza --i-data rep-seqs-380.qza \
    --i-data rep-seqs-382.qza --i-data rep-seqs-383.qza \
    --i-data rep-seqs-385.qza --i-data rep-seqs-386.qza \
    --i-data rep-seqs-387.qza --i-data rep-seqs-388.qza \
    --i-data rep-seqs-389.qza --o-merged-data new_nasal_rep-seqs_March2021.qza


# Clustering at 99% similarity using SILVA v138 reference database
## So nasal data can be merged with the gut data (uses different hypervariable region of 16S rRNA gene)
## also allows taxonomic assignment step to be done separately for each sample type and may improve consistency (more run to run variation with unique ASV level)
## uses SILVA v138 qza object from the QIIME2 website

qiime vsearch cluster-features-closed-reference \
    --i-sequences new_nasal_rep-seqs_March2021.qza \
    --i-table new_nasal_table_March2021.qza \
    --i-reference-sequences SILVA_138_Qiime2_FullLength_SSURef_NR99/silva-138-99-seqs.qza \
    --p-perc-identity 0.99 \
    --p-strand 'both' \
    --o-clustered-table table-cr-99_Mar2021.qza \
    --o-clustered-sequences rep-seqs-cr-99_Mar2021.qza \
    --o-unmatched-sequences unmatched-seqs_Mar2021 \
    --verbose

# Exporting table and sequences 
qiime tools export \
--input-path rep-seqs-cr-99_Mar2021.qza \
--output-path $PWD

qiime tools export \
--input-path table-cr-99_Mar2021.qza \
--output-path $PWD

# converting the qiime2 formatted SILVA taxonomy to txt format
qiime tools export \
    --input-path SILVA_138_Qiime2_FullLength_SSURef_NR99/silva-138-99-tax.qza \
    --output-path SILVA_138_Qiime2_FullLength_SSURef_NR99
## Get taxonomy.tsv output from this

# Add #OTUID and taxonomy as headers in SILVA taxonomy (taxonomy.tsv) before running this:
biom add-metadata -i feature-table.biom -o feature-table-tax.biom --observation-metadata-fp taxonomy.tsv --sc-separated taxonomy

biom convert --table-type="OTU table" -i feature-table-tax.biom -o feature-table-tax_nasal_March2021.txt --to-tsv --header-key taxonomy

## move main output data to be used as input for further preprocessing in R
mv dna-sequences.fasta ../nasal_input/
mv feature-table-tax_nasal_March2021.txt ../nasal_input/


################################################################################
# Infant stool preprocessing
# Demultiplexed fastq files, V4 hypervariable region of the 16S rRNA gene
################################################################################

#December run:
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path gut_manifests/MANIFEST_Dec_2 \
    --output-path Dec2-pa-demux.qza \
    --input-format PairedEndFastqManifestPhred33

#Quality check for each Run
qiime demux summarize --i-data Dec2-pa-demux.qza --o-visualization december-pa-demux.qzv 
qiime tools view december-pa-demux.qzv

#After viewing december-pa-demux.qzv and quality plots, decided to run dada2 using tuncation at 160 

#DADA2, the exact same parameters were used for each run 
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs Dec2-pa-demux.qza \
    --p-trunc-len-f 240 --p-trunc-len-r 160  --p-trim-left-f 10 --p-trim-left-r 10 \
    --o-table table_dec_160 --o-representative-sequences rep-seqs_dec_160 --o-denoising-stats sample_data_dec_160

#Checking the results of DADA2 for each run
qiime metadata tabulate --m-input-file sample_data_dec_160.qza --o-visualization denoising-stats_dec_160.qzv && qiime tools view denoising-stats_dec_160.qzv


#January run
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path gut_manifests/MANIFEST_Jan_2 \
    --output-path jan2-pa-demux.qza \
    --input-format PairedEndFastqManifestPhred33

qiime demux summarize --i-data jan2-pa-demux.qza --o-visualization jan-pa-demux.qzv

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs jan2-pa-demux.qza \
    --p-trunc-len-f 240 --p-trunc-len-r 160  --p-trim-left-f 10 --p-trim-left-r 10 \
    --o-table table_jan --o-representative-sequences rep-seqs_jan --o-denoising-stats sample_data_jan

qiime metadata tabulate --m-input-file sample_data_jan.qza --o-visualization denoising-stats_jan.qzv && qiime tools view denoising-stats_jan.qzv

#February run
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path gut_manifests/MANIFEST_Feb_2 \
  --output-path feb2-pa-demux.qza \
  --input-format PairedEndFastqManifestPhred33

qiime demux summarize --i-data feb2-pa-demux.qza --o-visualization feb-pa-demux.qzv

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs feb2-pa-demux.qza \
    --p-trunc-len-f 240 --p-trunc-len-r 160  --p-trim-left-f 10 --p-trim-left-r 10 \
    --o-table table_feb --o-representative-sequences rep-seqs_feb --o-denoising-stats sample_data_feb

qiime metadata tabulate --m-input-file sample_data_feb.qza --o-visualization denoising-stats_feb.qzv && qiime tools view denoising-stats_feb.qzv

#March 15 run
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path gut_manifests/MANIFEST_M15_2 \
  --output-path mar15_2-pa-demux.qza \
  --input-format PairedEndFastqManifestPhred33

qiime demux summarize --i-data mar15_2-pa-demux.qza --o-visualization mar15-pa-demux.qzv

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs mar15_2-pa-demux.qza \
    --p-trunc-len-f 240 --p-trunc-len-r 160  --p-trim-left-f 10 --p-trim-left-r 10 \
    --o-table table_mar15 --o-representative-sequences rep-seqs_mar15 --o-denoising-stats sample_data_mar15

qiime metadata tabulate --m-input-file sample_data_mar15.qza --o-visualization denoising-stats_mar15.qzv && qiime tools view denoising-stats_mar15.qzv

#March 20 run
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path gut_manifests/MANIFEST_M20_2 \
  --output-path mar20_2-pa-demux.qza \
  --input-format PairedEndFastqManifestPhred33


qiime demux summarize --i-data mar20_2-pa-demux.qza --o-visualization mar20-pa-demux.qzv

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs mar20_2-pa-demux.qza \
    --p-trunc-len-f 240 --p-trunc-len-r 160  --p-trim-left-f 10 --p-trim-left-r 10 \
    --o-table table_mar20 --o-representative-sequences rep-seqs_mar20 --o-denoising-stats sample_data_mar20

qiime metadata tabulate --m-input-file sample_data_mar20.qza --o-visualization denoising-stats_mar20.qzv && qiime tools view denoising-stats_mar20.qzv

#All Blanks (these were sent separately later on):
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path gut_manifests/MANIFEST_Blanks \
  --output-path blanks2-pa-demux.qza \
  --input-format PairedEndFastqManifestPhred33


qiime demux summarize --i-data blanks2-pa-demux.qza --o-visualization blanks-pa-demux.qzv

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs blanks2-pa-demux.qza \
    --p-trunc-len-f 240 --p-trunc-len-r 160  --p-trim-left-f 10 --p-trim-left-r 10 \
    --o-table table_blanks --o-representative-sequences rep-seqs_blanks --o-denoising-stats sample_data_blanks

qiime metadata tabulate --m-input-file sample_data_blanks.qza --o-visualization denoising-stats_blanks.qzv && qiime tools view denoising-stats_blanks.qzv


# Merge per-run tables are rep-seqs and visualize

qiime feature-table merge \
    --i-tables table_dec_160.qza --i-tables table_jan.qza \
    --i-tables table_feb.qza --i-tables table_mar15.qza \
    --i-tables table_mar20.qza --i-tables table_blanks.qza \
    --o-merged-table merged_table_gut.qza 

qiime feature-table merge-seqs \
    --i-data rep-seqs_dec_160.qza --i-data rep-seqs_jan.qza \
    --i-data rep-seqs_feb.qza --i-data rep-seqs_mar15.qza \
    --i-data rep-seqs_mar20.qza --i-data rep-seqs_blanks.qza \
    --o-merged-data rep-seqs_gut.qza

qiime feature-table summarize --i-table merged_table_gut.qza --o-visualization merged_table_gut.qzv


# Taxonomy assignment is closed reference so it can be merged with the nasal data (uses different hypervariable region of 16S rRNA gene)

qiime vsearch cluster-features-closed-reference \
    --i-sequences rep-seqs_gut.qza \
    --i-table merged_table_gut.qza \
    --i-reference-sequences 3_SILVA_138_Qiime2_FullLength_SSURef_NR99/silva-138-99-seqs.qza \
    --p-perc-identity 0.99 \
    --p-strand 'both' \
    --o-clustered-table gut-table-cr-99.qza \
    --o-clustered-sequences gut-rep-seqs-cr-99.qza \
    --o-unmatched-sequences unmatched-seqs \
    --verbose

# Exporting table and sequences 
qiime tools export \
    --input-path gut-rep-seqs-cr-99.qza \
    --output-path $PWD

qiime tools export \
    --input-path gut-table-cr-99.qza \
    --output-path $PWD

# converting the qiime2 formatted SILVA taxonomy to txt format
qiime tools export \
    --input-path SILVA_138_Qiime2_FullLength_SSURef_NR99/silva-138-99-tax.qza \
    --output-path SILVA_138_Qiime2_FullLength_SSURef_NR99
## Get taxonomy.tsv output from this

# Add #OTUID and taxonomy as headers in SILVA taxonomy:
biom add-metadata -i feature-table.biom -o gut-feature-table-tax.biom --observation-metadata-fp taxonomy.tsv --sc-separated taxonomy

biom convert --table-type="OTU table" -i gut-feature-table-tax.biom -o gut-feature-table-tax.txt --to-tsv --header-key taxonomy

## move main output data to be used as input for further preprocessing in R
mv dna-sequences.fasta ../gut_input/
mv gut-feature-table-tax.txt ../gut_input/