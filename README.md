# Tool repository for the identification of structurally critical amino acids in CDR1 and CDR2 antibodies


+ [Description of the project](#Description-of-the-project)
+ [Repository structure](#Repository-structure)
+ [AbScan](#AbScan)
+ [Examples](#Examples)
+ [Contacts](#Contacts)

## Description of the project
The goal of the project is to develop a tool for identifying structurally critical amino acids in CDR1 and CDR2 antibodies.

## Repository structure
data/ - folder with source data

data/median_structure_cdr/ - information about all CDR centroids

data/pairwise_align _canonical_form_and_median_structure/ - this folder contains all pairwise alignments of CDR sequences extracted from sabdab to CDR centroids. Each fasta file corresponds to one canonical form and its centroid structure

data/freq_aa_in_canonical_form/ - here is a json file with all amino acid frequencies for different families and positions


/The file "sabdab_summary_all_tsv" from the database "sabdab" contains the latest information for March 6th. The file sabdab_data_pdb_id.csv was obtained by removing all columns except pdb_id and heavy and light chain information from "sabdab_summary_all_tsv", then converted to CSV. Detailed information about the "sabdab_summary_all_tsv" structure can be found at https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/.

get_data_for_research/ - folder with data extraction code: sabdab(scv) => fasta => json | txt

## AbScan

/src - this scripts allows you to get information about the canonical form and frequency of amino acids present in antibodies. For this, files with data on the frequency of amino acids in canonical forms are used (AbScan/data
/median_and_frequency_data/) and the local version of SCALOP, which must be installed on the computer and added to the PATH variable.
```
usage: abscan_run.py [-h] [-s SEQUENCE] [-c CDR] [-n NUMBER] [-f FAMILY] [-m] [-i INPUT] [-o OUTPUT] [-v]

Retrieve the canonical form and/or frequencies of amino acids in an antibody sequence. To work with the amino acid sequence, you must have a local version of SCALOP installed on your device and a PATH
defined. if the chain parameter is not specified, then the canonical family will be determined using alignment. Otherwise, the family will be defined with SCALOP

options:
  -h, --help            show this help message and exit
  -s SEQUENCE, --sequence SEQUENCE
                        Amino acid sequence for which the canonical form and frequencies will be obtained
  -c CDR, --cdr CDR     CDR sequence for which frequency data and its corresponding canonical form will be retrieved
  -n NUMBER, --number NUMBER
                        The number of alternative amino acids to display (default=3)
  -f FAMILY, --family FAMILY
                        Path to csv file containing one sequence column or two sequence and circuit, the possible options are: L1, L2, L3, H1, H2
  -m, --mutant          Search for mutations that do not affect the canonical form
  -i INPUT, --input INPUT
                        File name/path for reading. The file must be in csv format, either with one sequence column, or with two sequence and chaine,the possible options for chain are: L1, L2, L3, H1,
                        H2
  -o OUTPUT, --output OUTPUT
                        File name/path for writing results
  -v, --variabel        modifier, used when the input file contains an antibody chain sequence


```
## Examples

### 1. Example using -s SEQUENCE:
```
request: python3 abscan_run.py -s EVQLVQPGAELRNSGASVKVSCKASGYRFTSYYIDWVRQAPGQGLEWMGRIDPEDGGTKYAQKFQGRVTFTADTSTSTAYVELSSLRSEDTAVYYCARNEWETVVVGDLMYEYEYWGQGTQVTVSSASTKGPSVFPLAPALGCLVKDYFPEPVTVSGVHTFPAVLQSSGLYSLSSVVNVNHK

answer:

Frequencies for CDR 'GYRFTSY': 

Canonical form for your sequence: H1-7-A

Amino acids frequencies:
+----------+------------+-----------+------------------------------+
| Position | Amino Acid | Frequency |         Alternatives         |
+----------+------------+-----------+------------------------------+
|    1     |     G      |   97.45%  |                              |
|          |            |           |  E(0.9%), D(0.7%), R(0.5%)   |
|    2     |     Y      |   22.05%  |                              |
|          |            |           |  F(58.9%), G(9.0%), D(2.1%)  |
|    3     |     R      |   1.94%   |                              |
|          |            |           | T(47.7%), S(17.3%), N(10.6%) |
|    4     |     F      |   66.51%  |                              |
|          |            |           | L(12.5%), I(11.3%), V(6.5%)  |
|    5     |     T      |   21.1%   |                              |
|          |            |           |  S(46.7%), K(5.7%), D(5.0%)  |
|    6     |     S      |   23.96%  |                              |
|          |            |           | D(22.1%), N(16.9%), G(8.5%)  |
|    7     |     Y      |   61.5%   |                              |
|          |            |           |  F(9.5%), S(7.4%), N(5.2%)   |
+----------+------------+-----------+------------------------------+

Frequencies for CDR 'DPEDGG': 

Canonical form for your sequence: H2-6-A

Amino acids frequencies:
+----------+------------+-----------+------------------------------+
| Position | Amino Acid | Frequency |         Alternatives         |
+----------+------------+-----------+------------------------------+
|    1     |     D      |   13.95%  |                              |
|          |            |           | N(28.6%), Y(19.2%), S(12.5%) |
|    2     |     P      |   70.03%  |                              |
|          |            |           |  T(11.0%), A(6.6%), W(3.2%)  |
|    3     |     E      |   4.15%   |                              |
|          |            |           | G(21.9%), Y(13.8%), N(12.4%) |
|    4     |     D      |   14.6%   |                              |
|          |            |           | N(26.5%), S(23.2%), T(13.6%) |
|    5     |     G      |   72.75%  |                              |
|          |            |           |  S(12.0%), D(8.7%), T(1.7%)  |
|    6     |     G      |   15.45%  |                              |
|          |            |           | N(15.9%), D(13.6%), Y(13.1%) |
+----------+------------+-----------+------------------------------+

```
### 2. Example using -c CDR:
```
request: python3 abscan_run.py -c SDRES

answer:

Canonical form for your sequence: L3-5-A

Amino acids frequencies:
+----------+------------+-----------+---------------------------------------------------------------------------+
| Position | Amino Acid | Frequency |                                Alternatives                               |
+----------+------------+-----------+---------------------------------------------------------------------------+
|    1     |     S      |    5.7%   |                                                                           |
|          |            |           |                         Q(88.6%), N(3.3%), W(1.4%)                        |
|    2*    |     D      |           |                                                                           |
|          |            |           |      Q(74.3%), V(10.0%), A(9.0%), H(2.9%), I(1.4%), T(1.4%), C(1.0%)      |
|    3*    |     R      |           |                                                                           |
|          |            |           | Y(49.0%), F(20.0%), L(16.2%), S(6.7%), P(4.8%), M(2.4%), I(0.5%), Q(0.5%) |
|    4     |     E      |   79.0%   |                                                                           |
|          |            |           |                         R(5.7%), Y(5.7%), L(4.3%)                         |
|    5     |     S      |    3.3%   |                                                                           |
|          |            |           |                        F(64.3%), T(18.1%), Y(6.2%)                        |
+----------+------------+-----------+---------------------------------------------------------------------------+

```

### 3. Example using -c CDR and -f FAMILY:
```
request: python3 abscan_run.py -c SDRES -f H2

answer:

Canonical form for your sequence: H2-5-A

Amino acids frequencies:
+----------+------------+-----------+------------------------------+
| Position | Amino Acid | Frequency |         Alternatives         |
+----------+------------+-----------+------------------------------+
|    1     |     S      |   20.19%  |                              |
|          |            |           | Y(27.6%), W(15.1%), T(12.0%) |
|    2     |     D      |   4.57%   |                              |
|          |            |           | S(25.6%), Y(24.7%), W(8.4%)  |
|    3     |     R      |   4.05%   |                              |
|          |            |           | G(31.4%), S(28.0%), D(21.0%) |
|    4     |     E      |   2.53%   |                              |
|          |            |           |  G(83.4%), D(9.4%), S(1.6%)  |
|    5     |     S      |   42.46%  |                              |
|          |            |           | N(12.1%), D(12.0%), T(11.8%) |
+----------+------------+-----------+------------------------------+

```


## Contacts
manasanartem4@gmail.com
