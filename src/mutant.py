import json
from typing import Optional, List
import time
import os
import logging

from combinations import get_combinations_fast, get_combinations_with_exclude
from positions_mutant_aa import find_mutation_position
from build_dict_from_scalop import get_dict_from_scalop
from build_frame import build_frame
from write_fasta_file import write_fasta_file

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_all_mutant(seq: str, chain: str, exclude_positions: List[int], output_file: str = "mutant_aa_result.json") -> Optional[None]:
    start_time = time.time()
    true_family = None
    true_cdr = None
    seq_for_scalop: str = build_frame(seq, chain)
    my_dict: dict = get_dict_from_scalop(seq_for_scalop, chain=chain)

    for key, value in my_dict.items():
        true_family = key
        true_cdr = value[0]['cdr_sequence']
        break

    if exclude_positions:
        combinations = get_combinations_with_exclude(seq, exclude_positions)
    else:
        combinations = get_combinations_fast(seq)
    logger.debug('List of combinations: %s', combinations)

    list_for_fasta = []

    for seq in combinations:
        frame_seq = build_frame(seq, chain)
        list_for_fasta.append(frame_seq)

    logger.info('Writing fasta file')
    write_fasta_file("sequences.fasta", list_for_fasta)
    mutant_dict = get_dict_from_scalop("sequences.fasta", chain=chain)

    logger.info('Building mutant dict')
    result: dict = {}
    for mutant_key, mutant_value in mutant_dict.items():
        for mutant in mutant_value:
            mutant_seq = mutant['cdr_sequence']
            if mutant_key == true_family and mutant_seq != true_cdr:
                mutation_pos = find_mutation_position(true_cdr, mutant_seq)
                if mutation_pos not in result:
                    result[mutation_pos] = []
                result[mutation_pos].append(mutant_seq)

    logger.info('Writing output file')
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=4)

    end_time = time.time()
    logger.info(f'Total time: {end_time - start_time}')

    os.remove("sequences.fasta")
    return None
