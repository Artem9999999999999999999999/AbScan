import warnings
from typing import Optional
from build_dict_from_scalop import get_dict_from_scalop
from build_frame import build_frame
from Bio import BiopythonDeprecationWarning

warnings.simplefilter('ignore', category=BiopythonDeprecationWarning)


def get_all_mutant(seq: str, chain: str, output_file: str = "mutant_aa_result.txt") -> Optional[None]:
    true_family: Optional[str] = None
    true_cdr: Optional[str] = None
    seq_for_scalop: str = build_frame(seq, chain)
    my_dict: dict = get_dict_from_scalop(seq_for_scalop, chain=chain)
    for key, value in my_dict.items():
        true_family = key
        true_cdr = value[0]['cdr_sequence']
        break

    result: set = set()

    def mutate(seq: str, pos: int, new_seq: str) -> None:
        # базовый случай: если обработали все позиции, добавляем мутантную последовательность в результат
        if pos == len(seq):
            mutant_seq_for_scalop = build_frame(new_seq, chain)
            mutant_dict = get_dict_from_scalop(mutant_seq_for_scalop, chain=chain)
            for mutant_key, mutant_value in mutant_dict.items():
                if mutant_key == true_family and mutant_value[0]['cdr_sequence'] != true_cdr and new_seq != seq:
                    result.add(new_seq)
            return

        # перебираем все возможные замены для текущей позиции
        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            if seq[pos] == aa:
                continue  # skip the same amino acid
            # рекурсивно вызываем функцию для следующей позиции
            mutate(seq, pos + 1, new_seq[:pos] + aa + new_seq[pos + 1:])

    mutate(seq, 0, seq)

    with open(output_file, 'w') as f:
        for item in result:
            f.write(item + '\n')

    return None
