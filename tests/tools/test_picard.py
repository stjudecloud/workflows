import pathlib
import pytest

@pytest.mark.workflow('picard_sort_queryname')
def test_picard_sort_queryname(workflow_dir):
    exists = pathlib.Path(workflow_dir, 'test-output/out/sorted_bam_index/test.bwa_aln_pe.chrY_chrM.sorted.bam.bai').exists()
    assert exists is False