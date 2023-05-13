from aligment import find_best_alignment
from build_dict_from_scalop import get_dict_from_scalop
from build_frame import build_frame
from csv_reader import read_csv_file
from csv_writer import data_write_in_json
from typing import List, Dict


def get_freq_in_file_for_cdr(file_input: str, file_output: str) -> None:
    sequences: List[Dict[str, str]] = read_csv_file(file_input)

    for seq in sequences:
        if seq["chain"] is not None:
            sequence: str = seq["sequence"]
            chain: str = seq["chain"]

            print(f"Accepted sequences: {seq['sequence']}, chain: {seq['chain']}")

            seq_for_scalop: str = build_frame(sequence, chain)
            my_dict: Dict[str, List[Dict[str, str]]] = get_dict_from_scalop(seq_for_scalop, chain=chain)

            if len(my_dict) == 0:
                print(
                    f"There is no frequency information for these sequences: {seq['sequence']} and chain: {seq['chain']}")
                print()

            for key, value in my_dict.items():
                true_family: str = key
                true_cdr: str = value[0]['cdr_sequence']

                print(f"Frequencies for CDR '{true_cdr} from {true_family} canonical form writing in output file")

                data_write_in_json(true_cdr, true_family, file_output)
        else:
            true_cdr: str = seq["sequence"]
            true_family: str = find_best_alignment(true_cdr)

            print(f"Accepted sequences: {seq['sequence']}")
            print(f"Frequencies for CDR {true_family} writing in output file")

            data_write_in_json(true_cdr, true_family, file_output)


def get_freq_in_file_for_seq(file_input: str, file_output: str) -> None:
    sequences: List[Dict[str, str]] = read_csv_file(file_input)

    for seq in sequences:
        sequence: str = seq["sequence"]

        print(f"Accepted sequences: {seq['sequence']}")

        my_dict: Dict[str, List[Dict[str, str]]] = get_dict_from_scalop(sequence, flag=True)
        for key, value in my_dict.items():
            true_family: str = key
            true_cdr: str = value[0]['cdr_sequence']
            print(f"Frequencies for CDR {true_cdr} from {true_family} canonical form writing in output file'")

            data_write_in_json(true_cdr, true_family, file_output)
