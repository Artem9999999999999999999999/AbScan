from typing import List


def get_combinations_fast(seq: str) -> List[str]:
    alphabet = "ACDEFGHIKLMNPQRSTVWY" # 20 стандартных аминокислот
    valid_combos = [] # список для хранения допустимых комбинаций

    for i in range(len(seq)):
        for a in alphabet:
            if a != seq[i]:
                combo = seq[:i] + a + seq[i+1:] # создаем новую комбинацию, заменяя текущую аминокислоту на a
                valid_combos.append(combo) # добавляем допустимую комбинацию в список

    return valid_combos
