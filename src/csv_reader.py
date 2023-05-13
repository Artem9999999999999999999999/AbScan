import csv
from typing import List, Dict, Union


def read_csv_file(file_path: str) -> List[Dict[str, Union[str, None]]]:
    sequences = []
    with open(file_path, 'r') as csv_file:
        reader = csv.reader(csv_file)
        headers = next(reader, None)  # читаем заголовки колонок
        for row in reader:
            if len(row) == 1:
                sequence = row[0]
                family = None
            else:
                sequence, family = row
            sequences.append({"sequence": sequence, "chain": family})
    return sequences
