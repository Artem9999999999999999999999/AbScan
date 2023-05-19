from typing import List


def get_combinations_fast(seq: str) -> List[str]:
    alphabet = "ACDEFGHIKLMNPQRSTVWY" 
    valid_combos = [] 

    for i in range(len(seq)):
        for a in alphabet:
            if a != seq[i]:
                combo = seq[:i] + a + seq[i+1:] 
                valid_combos.append(combo)

    return valid_combos


def get_combinations_with_exclude(seq: str, exclude_positions: List[int]) -> List[str]:
    alphabet = "ACDEFGHIKLMNPQRSTVWY" 
    valid_combos = [] 

    for i in range(len(seq)):
        if i + 1 in exclude_positions:  
            valid_combos.append(seq[i])  
        else:
            for a in alphabet:
                if a != seq[i]:
                    combo = seq[:i] + a + seq[i+1:] 
                    valid_combos.append(combo) 

    return valid_combos
