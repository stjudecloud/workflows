# Input Files used for testing

All the files in this directory were either randomly generated or sourced from publicly available genomic data. Many of the files have been downsampled to make them a manageable size for storing in a Git repository (using `git-lfs`) and to enable quick test execution in our CI. For the most part, these files are not representative of any meaningful scientific data and should be considered as mocked data only suitable for the purpose of generic testing.

The following list is sorted alphabetically:

## test.bwa_aln_pe.chrY_chrM.mutect2.vcf.gz

Unfiltered Mutect2 somatic VCF produced by running `gatk Mutect2` with `test/fixtures/bams/test.bwa_aln_pe.with_variants.chrY_chrM.bam` as the tumor sample and `test/fixtures/bams/test.bwa_aln_pe.chrY_chrM.bam` as the matched normal, against `test/fixtures/reference/GRCh38.chrY_chrM.fa`. Used as the `unfiltered_somatic_vcf` input for the `filter_mutect` task in `tools/mutect2.wdl`.

## test.bwa_aln_pe.chrY_chrM.mutect2.vcf.gz.tbi

Tabix index for `test.bwa_aln_pe.chrY_chrM.mutect2.vcf.gz`.

## test.bwa_aln_pe.chrY_chrM.mutect2.vcf.gz.stats

Mutect2 stats sidecar emitted alongside `test.bwa_aln_pe.chrY_chrM.mutect2.vcf.gz`. Used as the `unfiltered_somatic_vcf_stats` input for the `filter_mutect` task.

## test1.vcf.gz

Single position VCF file

## test1.vcf.gz.tbi

Index file for `test1.vcf.gz`.

## test2.vcf.gz

Single position VCF file

## test2.vcf.gz.tbi

Index file for `test2.vcf.gz`.

## testY.vcf.gz

Single position VCF file on `chrY`

## testY.vcf.gz.tbi

Index file for `testY.vcf.gz`.