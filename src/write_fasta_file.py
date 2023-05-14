from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def write_fasta_file(filename, sequences):
    records = [SeqRecord(Seq(seq), id=f"sequence{i+1}", description="") for i, seq in enumerate(sequences)]
    SeqIO.write(records, filename, "fasta")
