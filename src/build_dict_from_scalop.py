import warnings
from typing import Any, Dict, List

from Bio import BiopythonDeprecationWarning

warnings.simplefilter('ignore', category=BiopythonDeprecationWarning)

import subprocess


def get_dict_from_scalop(seq: str, flag: bool = False, chain: str = ""):
    args: List[str] = ["SCALOP", "-i", seq, "--numbering_scheme", "chothia", "--cdr_definition", "chothia"]

    try:
        output = subprocess.run(args, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        warnings.warn(f"An error occurred: {e}")

    my_dict: Dict[str, List[Dict[str, Any]]] = {}

    for line in output.stdout.split("\n")[3:]:
        fields = line.strip().split()
        if len(fields) == 5:
            input_id = fields[0]
            cdr_id = fields[1]
            cdr_sequence = fields[2]
            canonical_id = fields[3]
            median_id = fields[4]

            if flag:
                new_entry = {
                    "cdr_sequence": cdr_sequence,
                }

                if canonical_id not in my_dict:
                    my_dict[canonical_id] = []
                my_dict[canonical_id].append(new_entry)
            else:
                if canonical_id[:2] == chain:
                    new_entry = {
                        "cdr_sequence": cdr_sequence,
                    }

                    if canonical_id not in my_dict:
                        my_dict[canonical_id] = []
                    my_dict[canonical_id].append(new_entry)
    return my_dict
