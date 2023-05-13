from build_dict_from_scalop import get_dict_from_scalop
from print_table import get_freq_from_cdr
from build_frame import build_frame


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
