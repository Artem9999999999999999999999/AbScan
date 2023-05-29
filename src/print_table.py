import warnings
from typing import Optional
from Bio import BiopythonDeprecationWarning

warnings.simplefilter('ignore', category=BiopythonDeprecationWarning)
import json
from prettytable import PrettyTable
from aligment import find_best_alignment


def get_freq_from_cdr(seq: str, family: str = "NA", number_alter: int = 3) -> Optional[None]:
    with open('../data_absight/freq_aa_in_canonical_form.json') as f:
        data = json.load(f)

    if family == "NA":
        family = find_best_alignment(seq)
        if family == "No alignment found":
            warnings.warn("Sorry, canonical forms corresponding to your CDR could not be found")
            return None

    print(f"\nCanonical form for your sequence: {family}")
    seq_len = len(seq)
    try:
        freq: dict = {}
        for i, aa in enumerate(seq):
            for d in data[family]:
                if str(i + 1) in d:
                    if aa in d[str(i + 1)]:
                        freq[i + 1] = d[str(i + 1)][aa]
                    else:
                        freq[i + 1] = None
                    break
            else:
                freq[i + 1] = None

        print()
        print("Amino acids frequencies:")
        table = PrettyTable()
        table.field_names = ["Position", "Amino Acid", "Frequency", "Alternatives"]
        for k, v in freq.items():
            alternatives: list = []
            if v is None:
                for d in data[family]:
                    if str(k) in d:
                        amino_acid = seq[k - 1]
                        row = ["{}*".format(k), amino_acid, "", ""]
                        table.add_row(row)
                        for aa_key, freq_value in d[str(k)].items():
                            if aa_key == amino_acid:
                                continue
                            alternatives.append((aa_key, freq_value))
                        break
            else:
                amino_acid = seq[k - 1]
                row = [k, amino_acid, "{}%".format(v), ""]
                table.add_row(row)
                for aa_key, freq_value in sorted(data[family][0][str(k)].items(), key=lambda x: x[1], reverse=True):
                    if aa_key == amino_acid:
                        continue
                    if len(alternatives) >= number_alter:
                        break
                    alternatives.append((aa_key, freq_value))
            alternatives_str = ", ".join(["{}({:.1f}%)".format(aa, freq) for aa, freq in alternatives])
            row = ["", "", "", alternatives_str]
            table.add_row(row)
        print(table)

        if seq_len > len(data[family][0]):
            print()
            if seq_len == len(data[family][0]) + 1:
                warnings.warn(
                    f"There is no information for the amino acid at the position {len(data[family][0]) + 1}, most "
                    f"likely something is wrong with your CDR =(")
            else:
                warnings.warn(
                    f"There is no information for amino acids in positions from {len(data[family][0]) + 1} to {seq_len}, most likely something is wrong with your CDR =(")

    except KeyError:
        warnings.warn("Error: key not found.")
