import pytest
import pathlib
from collections import OrderedDict

@pytest.mark.workflow('hilow_filter')
def test_hilow_filter(workflow_dir):
    filtered_pairs = pathlib.Path(workflow_dir, 'test-output/out/filtered_pairs/00_JurkatCas9_ab3_050123.subset.allValidPairs.allValidPairs.filtered')
    removed_pairs = pathlib.Path(workflow_dir, 'test-output/out/removed_pairs/00_JurkatCas9_ab3_050123.subset.allValidPairs.allValidPairs.removed')
    stats = pathlib.Path(workflow_dir, 'test-output/out/filtered_stats/00_JurkatCas9_ab3_050123.subset.allValidPairs.allValidPairs.stats')
    # Read stats file
    with open(stats) as f:
        stats = f.read()
    assert stats == "978 read pairs are filtered out from a total of 350000, accounting for .27942857142857142800 %\n\n"
    
@pytest.mark.workflow('hilow_qcreport')
def test_hilow_qcreport(workflow_dir):
    qcreport = pathlib.Path(workflow_dir, 'test-output/out/qc_report/00_JurkatCas9_ab3_050123_QCreport.txt')
    # Read qcreport file
    with open(qcreport) as f:
        qcreport = f.read()
    assert qcreport == 'STAT\tCOUNTS\tPERCENTAGE\tVERDICT\nTotal_pairs_processed\t30000000\nR1_aligned\t20644728\t69\tBAD\nR2_aligned\t19932436\t66\tBAD\nvalid_interactionPairs\t7968247\t82\tGOOD\ncis_shortRange\t3244623\t41\tMARGINAL\ncis_longRange\t3367790\t42\tGOOD\n'