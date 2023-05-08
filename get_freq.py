import warnings
from typing import Any, Optional, Dict, List, Tuple, Union

from Bio import BiopythonDeprecationWarning

warnings.simplefilter('ignore', category=BiopythonDeprecationWarning)

import json
import subprocess
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from prettytable import PrettyTable
import argparse


def find_best_alignment(seq: str) -> Union[int, Any]:
    with open('median_structure_cdr.json') as f:
        data_median = json.load(f)

    matrix = substitution_matrices.load("BLOSUM62")
    scores = {}
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


def get_freq_from_cdr(seq: str, family: str = "NA", number_alter: int = 3) -> Optional[None]:
    with open('freq_aa_in_canonical_form.json') as f:
        data = json.load(f)

    if family == "NA":
        family = find_best_alignment(seq)
        if family == "No alignment found":
            warnings.warn("Sorry, canonical forms corresponding to your CDR could not be found")
            return None

    print(f"\nCanonical form for your sequence: {family}")
    seq_len = len(seq)
    try:
        freq = {}
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
            alternatives = []
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
                    f"There is no information for the amino acid at the position {len(data[family][0]) + 1}, most likely something is wrong with your CDR =(")
            else:
                warnings.warn(
                    f"There is no information for amino acids in positions from {len(data[family][0]) + 1} to {seq_len}, most likely something is wrong with your CDR =(")

    except KeyError:
        warnings.warn("Error: key not found.")


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


def build_frame(seq: str, chain: str) -> str:
    seq_for_scalop = ""
    if chain == "L1":
        seq_for_scalop = "NFMLTQPHSVSESPGKTVTISC" + seq + "WYQRRPGSAPTTVIYEDNQRPSGVPDRFSASIDSSSNSASLTISGLKTEDEADYYCQSYDSSNWVFGGGTKLTVLGQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
    elif chain == "L2":
        seq_for_scalop = "SDVGTFDLVSWYQQYPGKAPKLIIY" + seq + "GVSDRFSGSKSGNTASLTISGLQAEDEADYYCSSYAGSVVFGGGTKLTVLGQPKGAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
    elif chain == "L3":
        seq_for_scalop = "SDVGTFDLVSWYQQYPGKAPKLIIYEGSRRPSGVSDRFSGSKSGNTASLTISGLQAEDEADYYC" + seq + "FGGGTKLTVLGQPKGAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
    elif chain == "H1":
        seq_for_scalop = "VQLQESGGCLVQAGGSLRLSCAAS" + seq + "TIGWFRQAPGKEREFVAAIHWDGGQTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAARGRRYFDFTYSDVYDYWGQGTQVTVS"
    elif chain == "H2":
        seq_for_scalop = "VQLQESGGCLVQAGGSLRLSCAASGSTFSTYTIGWFRQAPGKEREFVAAI" + seq + "TYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAARGRRYFDFTYSDVYDYWGQGTQVTVS"
    return seq_for_scalop


def get_all_mutant(seq: str, chain: str, output_file: str = "mutant_aa_result.txt") -> None:
    def find_true_family_and_cdr(seq: str, chain: str) -> Tuple[str, str]:
        seq_for_scalop = build_frame(seq, chain)
        my_dict = get_dict_from_scalop(seq_for_scalop, chain=chain)
        for key, value in my_dict.items():
            return key, value[0]['cdr_sequence']
        return '', ''

    true_family, true_cdr = find_true_family_and_cdr(seq, chain)
    result = set()

    def mutate_helper(seq: List[str], pos: int, new_seq: List[str]) -> None:
        if pos == len(seq):
            mutant_seq_for_scalop = ''.join(new_seq)
            mutant_dict = get_dict_from_scalop(mutant_seq_for_scalop, chain=chain)
            for mutant_key, mutant_value in mutant_dict.items():
                if mutant_key == true_family and mutant_value[0]['cdr_sequence'] != true_cdr and new_seq != seq:
                    result.add(''.join(new_seq))
            return

        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            if seq[pos] == aa:
                continue
            new_seq[pos] = aa
            mutate_helper(seq, pos + 1, new_seq)

    seq = list(seq)
    mutate_helper(seq, 0, seq)

    with open(output_file, 'w') as f:
        for item in result:
            f.write(item + '\n')



def get_freq_from_scalop(seq_for_scalop: str, number_alter: int = 3) -> None:
    my_dict = get_dict_from_scalop(seq_for_scalop, flag=True)

    for key, value in my_dict.items():
        print()
        print(f"Frequencies for CDR '{value[0]['cdr_sequence']}': ")
        get_freq_from_cdr(value[0]['cdr_sequence'], key, number_alter=number_alter)


def get_freq_from_scalop_with_chain(seq: str, chain: str, number_alter: int = 3) -> None:
    seq_for_scalop = build_frame(seq, chain)
    my_dict = get_dict_from_scalop(seq_for_scalop, chain=chain)

    for key, value in my_dict.items():
        get_freq_from_cdr(value[0]['cdr_sequence'], key, number_alter=number_alter)


def main():
    parser = argparse.ArgumentParser(
        description='Retrieve the canonical form and/or frequencies of amino acids in an antibody sequence. To work '
                    'with the amino acid sequence, you must have a local version of SCALOP installed on your device '
                    'and a PATH defined. if the chain parameter is not specified, then the canonical family will be '
                    'determined using alignment. Otherwise, the family will be defined with SCALOP')
    parser.add_argument('-s', '--sequence',
                        help='amino acid sequence for which the canonical form and frequencies will be obtained')
    parser.add_argument('-c', '--cdr',
                        help='CDR sequence for which frequency data and its corresponding canonical form will be '
                             'retrieved')
    parser.add_argument('-n', '--number',
                        help='the number of alternative amino acids to display (default=3)', type=int, default=3)
    parser.add_argument('-f', '--family',
                        help='enter the name of the circuit, the possible options are: L1, L2, L3, H1, H2')
    parser.add_argument('-m', '--mutant',
                        help='Search for mutations that do not affect the canonical form', action='store_true')
    parser.add_argument('-o', '--output',
                        help='File name for writing results', type=str, default='mutant_aa_result.txt')

    args = parser.parse_args()

    if not args.sequence and not args.cdr:
        parser.print_help()
        return

    if args.sequence:
        seq_for_scalop = args.sequence
        number = args.number
        get_freq_from_scalop(seq_for_scalop, number)

    if args.cdr and not args.family:
        seq = args.cdr
        number = args.number
        get_freq_from_cdr(seq, number_alter=number)

    if args.cdr and args.family and args.family:
        seq = args.cdr
        number = args.number
        chain = args.family
        get_freq_from_scalop_with_chain(seq, chain, number_alter=number)

    if args.mutant and args.cdr and args.family:
        seq = args.cdr
        output_filename = args.output
        chain = args.family
        get_all_mutant(seq, chain, output_filename)


if __name__ == "__main__":
    main()



