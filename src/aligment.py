import warnings
from typing import Any

from Bio import BiopythonDeprecationWarning

warnings.simplefilter('ignore', category=BiopythonDeprecationWarning)

import json
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.Seq import Seq


def find_best_alignment(seq: str) -> int | Any:
    with open('../data_absight/median_structure_cdr.json') as f:
        data_median = json.load(f)

    matrix = substitution_matrices.load("BLOSUM62")
    scores: dict = {}
    for key, value in data_median.items():
        seq1 = Seq(str(seq))
        seq2 = Seq(str(value))
        alignment = pairwise2.align.globalds(seq1, seq2, matrix, -26, -0.5)
        alignment_score = alignment[0].score
        if alignment_score > -24:
            scores[key] = alignment_score

    if not len(scores):
        warnings.warn("No alignment found")
        return

    best_key = max(scores, key=scores.get)
    return best_key
