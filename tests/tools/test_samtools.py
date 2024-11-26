import pytest
import pathlib
from collections import OrderedDict

import pysam
import fastq


@pytest.mark.workflow('samtools_split')
def test_samtools_split(workflow_dir):
    bam = pathlib.Path(workflow_dir, 'test-output/out/split_bams/0/test.1.bam')
    samfile = pysam.AlignmentFile(bam, "rb")
    bam_header = OrderedDict((k, v) for k, v in samfile.header.items())
    read_groups = [read_group['ID'] for read_group in bam_header.get('RG', []) if 'ID' in read_group]
    assert len(read_groups) == 1
    assert read_groups[0] == "1"

    second_bam = pathlib.Path(workflow_dir, 'test-output/out/split_bams/1/test.2.bam')
    second_samfile = pysam.AlignmentFile(second_bam, "rb")
    second_bam_header = OrderedDict((k, v) for k, v in second_samfile.header.items())
    second_read_groups = [read_group['ID'] for read_group in second_bam_header.get('RG', []) if 'ID' in read_group]
    assert len(second_read_groups) == 1
    assert second_read_groups[0] == "2"

@pytest.mark.workflow('samtools_merge')
def test_samtools_merge(workflow_dir):
    bam = pathlib.Path(workflow_dir, 'test-output/out/merged_bam/test.bam')
    samfile = pysam.AlignmentFile(bam, "rb")
    bam_header = OrderedDict((k, v) for k, v in samfile.header.items())
    read_groups = [read_group['ID'] for read_group in bam_header.get('RG', []) if 'ID' in read_group]
    assert len(read_groups) == 2
    assert read_groups[0] == "test2"
    assert read_groups[1] == "test.bwa_aln_pe"

@pytest.mark.workflow('samtools_collate', 'samtools_collate_to_fastq')
def test_samtools_collate(workflow_dir):
    bam = pathlib.Path(workflow_dir, 'test-output/out/collated_bam/test.bwa_aln_pe.chrY_chrM.collated.bam')
    samfile = pysam.AlignmentFile(bam, "rb")

    reads = list(samfile.fetch(until_eof=True))
    for c in range(0, 100, 2):
        assert reads[c].query_name == reads[c+1].query_name
        assert reads[c].is_read1 != reads[c+1].is_read1

@pytest.mark.workflow('samtools_bam_to_fastq', 'samtools_collate_to_fastq')
def test_samtools_bam_to_fastq(workflow_dir):
    fq1 = fastq.read(pathlib.Path(workflow_dir, 'test-output/out/read_one_fastq_gz/test.bwa_aln_pe.chrY_chrM.R1.fastq.gz'))
    fq2 = fastq.read(pathlib.Path(workflow_dir, 'test-output/out/read_two_fastq_gz/test.bwa_aln_pe.chrY_chrM.R2.fastq.gz'))

    for r1, r2 in zip(fq1, fq2):
        assert r1.head.removesuffix("/1") == r2.head.removesuffix("/2")
