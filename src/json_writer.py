import json
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def data_write_in_json(cdr: str, family: str, file_output: str) -> None:
    with open('../data_abscan/freq_aa_in_canonical_form.json', 'r') as f:
        data = json.load(f)
        
    data_for_true_family = data.get(family, [])
    data_to_write = {cdr: data_for_true_family}

    with open(file_output, 'a') as f:
        json.dump(data_to_write, f, indent=4)
        f.write('\n')
        logger.info(f"Data for {family} has been written to {file_output} file.\n")
