def find_mutation_position(true_seq: str, mutant_seq: str) -> int:
    """
    Finds the position of the mutation between two sequences of equal length.
    """
    if len(true_seq) != len(mutant_seq):
        raise ValueError("Sequences must be of equal length to compare.")

    for i in range(len(true_seq)):
        if true_seq[i] != mutant_seq[i]:
            return i+1

    return -1