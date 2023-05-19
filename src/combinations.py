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


def get_combinations_with_exclude(seq: str, exclude_positions: List[int]) -> List[str]:
    alphabet = "ACDEFGHIKLMNPQRSTVWY" # 20 стандартных аминокислот
    valid_combos = [] # список для хранения допустимых комбинаций

    for i in range(len(seq)):
        if i + 1 in exclude_positions:  # проверяем, является ли текущая позиция запрещенной
            valid_combos.append(seq[i])  # если да, то добавляем текущую аминокислоту без изменений
        else:
            for a in alphabet:
                if a != seq[i]:
                    combo = seq[:i] + a + seq[i+1:] # создаем новую комбинацию, заменяя текущую аминокислоту на a
                    valid_combos.append(combo) # добавляем допустимую комбинацию в список

    return valid_combos
