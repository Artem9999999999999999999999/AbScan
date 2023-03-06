import urllib.request
import urllib.error
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1
import pandas as pd

# Загрузить данные sabdab в DataFrame
df = pd.read_csv('sabdab_data_pdb_id2.csv')

# Открыть файл для записи результатов
with open('all_sequences.fasta', 'w') as output_file:

    # Пройти по всем строкам DataFrame
    for index, row in df.iterrows():

        # Получить PDB ID и цепочки
        pdb_id = row['pdb']
        chains = f"{row['Hchain']},{row['Lchain']}"

        # Попытаться загрузить структуру
        try:
            # Скачать PDB файл и прочитать его
            pdb_filename = f'{pdb_id}.pdb'
            pdb_url = f'https://files.rcsb.org/download/{pdb_filename}'
            urllib.request.urlretrieve(pdb_url, pdb_filename)
            pdb_parser = PDBParser(QUIET=True)
            structure = pdb_parser.get_structure(pdb_id, pdb_filename)

            # Пройти по всем цепочкам антитела
            for chain_id in chains:
                # Получить последовательность аминокислот цепочки
                try:
                    chain = structure[0][chain_id]
                    seq = seq1("".join([residue.resname for residue in chain]))
                    # Записать последовательность в файл
                    output_file.write(f'>{pdb_id}_{chain_id}\n{seq}\n')
                except KeyError:
                    print(f'Chain {chain_id} not found in PDB file {pdb_id}')
        except urllib.error.HTTPError:
            print(f'PDB ID {pdb_id} not found on RCSB')
