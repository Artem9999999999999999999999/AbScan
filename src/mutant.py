import time
import logging
from typing import Optional

from build_dict_from_scalop import get_dict_from_scalop
from build_frame import build_frame
from itertools import product


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def generate_combinations(input_string, alphabet):
    possible_replacements = [alphabet] * len(input_string)
    combinations = product(*possible_replacements)
    return [''.join(c) for c in combinations if c != input_string]


def get_all_mutant(seq: str, chain: str, output_file: str = "mutant_aa_result.txt") -> Optional[None]:
    start_time = time.time()
    true_family = None
    true_cdr = None
    seq_for_scalop: str = build_frame(seq, chain)
    my_dict: dict = get_dict_from_scalop(seq_for_scalop, chain=chain)

    for key, value in my_dict.items():
        true_family = key
        true_cdr = value[0]['cdr_sequence']
        break

    ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

    combinations = generate_combinations(seq, ALPHABET)
    logger.debug('List of combinations: %s', combinations)

    logger.info('Working with SCALOP')
    result: set = set()
    for seq in combinations:
        frame_seq = build_frame(seq, chain)
        mutant_dict = get_dict_from_scalop(frame_seq, chain=chain)

        for mutant_key, mutant_value in mutant_dict.items():
            mutant_cdr = mutant_value[0]['cdr_sequence']
            if mutant_key == true_family and mutant_cdr != true_cdr:
                result.add(mutant_cdr)

    logger.info('Writing output file')
    with open(output_file, 'w') as f:
        for item in result:
            f.write(item + '\n')

    end_time = time.time()
    logger.info(f'Total time: {end_time - start_time}')

    return None

