import os

from aligment import find_best_alignment
from build_dict_from_scalop import get_dict_from_scalop
from build_frame import build_frame
from csv_reader import read_csv_file
from json_writer import data_write_in_json
from typing import List, Dict
from write_fasta_file import write_fasta_file


def get_freq_in_file(file_input: str, file_output: str) -> None:
    sequences: List[Dict[str, str]] = read_csv_file(file_input)

    first_sequence = sequences[0]["sequence"]
    length_of_first_sequence = len(first_sequence)

    list_for_fasta = []
    for seq in sequences:
        sequence: str = seq["sequence"]
        if seq["chain"] is not None:

            chain: str = seq["chain"]

            print(f"Accepted sequences: {seq['sequence']}, chain: {seq['chain']}")

            seq_for_scalop: str = build_frame(sequence, chain)
            list_for_fasta.append(seq_for_scalop)

            write_fasta_file("sequences.fasta", list_for_fasta)

            my_dict: Dict[str, List[Dict[str, str]]] = get_dict_from_scalop("sequences.fasta", chain=chain)
            os.remove("sequences.fasta")
            if len(my_dict) == 0:
                print(
                    f"There is no frequency information for these sequences: {seq['sequence']} and chain: {seq['chain']}")
                print()

            for key, value in my_dict.items():
                true_family: str = key
                true_cdr: str = value[0]['cdr_sequence']

                print(f"Frequencies for CDR '{true_cdr} from {true_family} canonical form writing in output file'")

                data_write_in_json(true_cdr, true_family, file_output)
        elif length_of_first_sequence > 20:

            print(f"Accepted sequences: {seq['sequence']}")

            list_for_fasta.append(sequence)
            write_fasta_file("sequences.fasta", list_for_fasta)

            my_dict: Dict[str, List[Dict[str, str]]] = get_dict_from_scalop("sequences.fasta", flag=True)
            os.remove("sequences.fasta")
            for key, value in my_dict.items():
                true_family: str = key
                true_cdr: str = value[0]['cdr_sequence']
                print(f"Frequencies for CDR {true_cdr} from {true_family} canonical form writing in output file'")

                data_write_in_json(true_cdr, true_family, file_output)

        else:
            true_cdr: str = seq["sequence"]
            true_family: str = find_best_alignment(true_cdr)

            print(f"Accepted sequences: {seq['sequence']}")
            print(f"Frequencies for CDR {true_family} writing in output file")

            data_write_in_json(true_cdr, true_family, file_output)


