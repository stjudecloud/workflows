# Mutect2 fixtures

Inputs and intermediate artifacts used by the per-task tests in `tools/test/mutect2.yaml` for `tools/mutect2.wdl`. Files are downsampled to `chrY` and `chrM` to match the reference at `test/fixtures/reference/GRCh38.chrY_chrM.fa` and the BAMs at `test/fixtures/bams/test.bwa_aln_pe*.chrY_chrM.bam`.

The following list is sorted alphabetically:

## af-only-gnomad.hg38.chrY_chrM.vcf.gz

GATK best-practices germline resource (`af-only-gnomad.hg38.vcf.gz` from the GATK resource bundle), subset to `chrY` and `chrM`. Used as the `germline_resource_vcf` input for the `mutect2` task and as both `intervals` and `variants` inputs for the `get_pileup_summaries` task.

## af-only-gnomad.hg38.chrY_chrM.vcf.gz.tbi

Tabix index for `af-only-gnomad.hg38.chrY_chrM.vcf.gz`.

## test.bwa_aln_pe.chrY_chrM_pileup_summaries.table

Output of `gatk GetPileupSummaries` run on `test/fixtures/bams/test.bwa_aln_pe.chrY_chrM.bam` (normal sample) over the sites in `af-only-gnomad.hg38.chrY_chrM.vcf.gz`. Used as the `normal_pileups` input for the `calculate_contamination` task.

## test.bwa_aln_pe.with_variants.chrY_chrM_pileup_summaries.contamination.table

Output of `gatk CalculateContamination` run on the tumor and normal pileup tables (matched-normal mode). Used as the `contamination_table` input for the `filter_mutect` task.

## test.bwa_aln_pe.with_variants.chrY_chrM_pileup_summaries.segments.table

Tumor segmentation output of the same `gatk CalculateContamination` invocation that produced the contamination table. Used as the `maf_segments` input for the `filter_mutect` task.

## test.bwa_aln_pe.with_variants.chrY_chrM_pileup_summaries.table

Output of `gatk GetPileupSummaries` run on `test/fixtures/bams/test.bwa_aln_pe.with_variants.chrY_chrM.bam` (tumor sample) over the sites in `af-only-gnomad.hg38.chrY_chrM.vcf.gz`. Used as the `tumor_pileups` input for the `calculate_contamination` task.
