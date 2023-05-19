import json
import time
import logging
from typing import Optional, List

from positions_mutant_aa import find_mutation_position
from combinations import get_combinations_fast, get_combinations_with_exclude
from aligment import find_best_alignment


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_all_mutant_with_aligment(seq: str, exclude_positions: List[int], output_file: str = "result.json") -> Optional[None]:
    start_time = time.time()
    true_cdr = seq
    true_family = find_best_alignment(seq)

    if exclude_positions:
        combinations = get_combinations_with_exclude(seq, exclude_positions)
    else:
        combinations = get_combinations_fast(seq)
    logger.debug('List of combinations: %s', combinations)

    logger.info('Building mutant dict')
    result: dict = {}
    for combination in combinations:
        mutant_cdr = combination
        mutant_family = find_best_alignment(combination)
        if mutant_family == true_family and mutant_cdr != true_cdr:
            mutation_pos = find_mutation_position(true_cdr, mutant_cdr)
            if mutation_pos not in result:
                result[mutation_pos] = []
            result[mutation_pos].append(mutant_cdr)

    logger.info('Writing output file')
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=4)

    end_time = time.time()
    logger.info(f'Total time: {end_time - start_time}')

    return None

