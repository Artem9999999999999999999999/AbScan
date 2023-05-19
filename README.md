# AbScan

+ [Description of the project](#Description-of-the-project)
+ [Installation instructions](#Installation-instructions)
+ [Repository structure](#Repository-structure)
+ [AbScan run](#AbScan-run)
+ [Examples usage](#Examples-usage)
+ [Examples input csv file](#Examples-input-csv-file)
+ [Contacts](#Contacts)

# Tool repository for the identification of structurally critical amino acids in CDR1 and CDR2 antibodies
## Description of the project
The goal of the project is to develop a tool for identifying structurally critical amino acids in CDR1 and CDR2 antibodies.

## Installation instructions

Download the zip archive.

To start, from **/src**, run the command:  **chmod +x AbScan**

Enjoy

## Repository structure
**/data** - folder with source data

**/data/pairwise_align _canonical_form_and_median_structure/** - this folder contains all pairwise alignments of CDR sequences extracted from sabdab to CDR centroids. Each fasta file corresponds to one canonical form and its centroid structure

**data_abscan/freq_aa_in_canonical_form/** - here is a json file with all amino acid frequencies for different families and positions

**data_abscan/median_structure_cdr/** - information about all CDR centroids

The file **"sabdab_summary_all_tsv"** from the database "sabdab" contains the latest information for March 6th. The file sabdab_data_pdb_id.csv was obtained by removing all columns except pdb_id and heavy and light chain information from **"sabdab_summary_all_tsv"**, then converted to CSV. Detailed information about the **"sabdab_summary_all_tsv"** structure can be found at https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/.

**get_data_for_research/** - folder with data extraction code: sabdab(scv) => fasta => json | txt

## AbScan run

**src/AbScan** - this scripts allows you to get information about the canonical form and frequency of amino acids present in antibodies. For this, files with data on the frequency of amino acids in canonical forms are used (**data_abscan/**) and the local version of SCALOP, which must be installed on the computer and added to the PATH variable.

The file AbScan is used to run the script.
```
usage: AbScan [-h] [-c CDR] [-n NUMBER] [-f FAMILY] [-m] [-p [POSITIONS ...]] [-o OUTPUT]

Retrieve the canonical form and/or frequencies of amino acids in an antibody sequence. To work with the amino acid sequence, you must have a local version of SCALOP installed on your device and a PATH
defined. if the chain parameter is not specified, then the canonical family will be determined using alignment. Otherwise, the family will be defined with SCALOP

options:
  -h, --help            show this help message and exit
  -c CDR, --cdr CDR     CDR sequence for which frequency data and its corresponding canonical form will be retrieved
  -n NUMBER, --number NUMBER
                        The number of alternative amino acids to display (default=3)
  -f FAMILY, --family FAMILY
                        Path to csv file containing one sequence column or two sequence and circuit, the possible options are: L1, L2, L3, H1, H2
  -m, --mutant          Search for mutations that do not affect the canonical form
  -p [POSITIONS ...], --positions [POSITIONS ...]
                        List of positions that should not undergo mutations (1-20)
  -o OUTPUT, --output OUTPUT
                        File name/path for writing results. If not specified, the default path is "output_result/result.json"

```
## Examples usage


### 1. Example using -c CDR and -f FAMILY:

This function accepts a cdr sequence, a chain indication and a loop number

##### Recommended for use.

```
request: ./AbScan -c SDRES -f H2

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

### 2. Example using -c SEQUENCE:

This function can accept antibody chains

```
request: ./AbScan -c EVQLVQPGAELRNSGASVKVSCKASGYRFTSYYIDWVRQAPGQGLEWMGRIDPEDGGTKYAQKFQGRVTFTADTSTSTAYVELSSLRSEDTAVYYCARNEWETVVVGDLMYEYEYWGQGTQVTVSSASTKGPSVFPLAPALGCLVKDYFPEPVTVSGVHTFPAVLQSSGLYSLSSVVNVNHK

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
### 3. Example using -c CDR:

In the case where the chain argument is not specified, an alignment algorithm is used whose accuracy is 76% compared to the SCALOP results.

```
request: ./AbScan -c SDRES

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

### 4. Example using -c input_file.csv [-o output_filename.json]

Сsv files are supported. As an output, a json file that contains information about the frequencies of all amino acids for each position of your CDR. To change the default output filename use the -o argument

```
request: ./AbScan -c input_file.csv -o output_filename.json

answer: output_filename.json

```
### 5. Example using -c CDR -f FAMILY -m [-o output_filename.json]

This function determines which substitutions do not change the canonical form. No more than one substitution is made at a time. The output is a json file whose keys are positions and their values are sequences of CDRs of the same canonical form, replaced at that position.


##### Recommended for use.

```
request: ./AbScan -c SDRES -f H2 -m -o output_filename.json

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
INFO:mutant:Writing fasta file
INFO:mutant:Building mutant dict
INFO:mutant:Writing output file
INFO:mutant:Total time: 1.2449893951416016

output_filename.json
```

### 6. Example using -c CDR -m [-o output_filename.json]

the function is similar to the previous one, the difference is the lack of information about the chain. To determine the canonical forms, alignment is used, which is a less accurate method.

```
request: ./AbScan -c SDRES -m -o output_filename.json

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
INFO:mutant_with_aligment:Building mutant dict
INFO:mutant_with_aligment:Writing output file
INFO:mutant_with_aligment:Total time: 8.164558410644531

output_filename.json
```

### 7. Example using -c CDR -m -p POSITIONS [-o output_filename.json]

This function allows you not to mutate amino acids in selected positions

```
request: ./AbScan -c SDRES -f H2 -m -p 2 3 5 output_filename.json

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
INFO:mutant:Writing fasta file
INFO:mutant:Building mutant dict
INFO:mutant:Writing output file
INFO:mutant:Total time: 0.9216213226318359

output_filename.json
```


## Examples input csv file 

Sample files are stored in the **file_examples/** folder

### 1. Example csv file with chain:

```
sequence,chain
SDRES,H2
QAWDSSTAWV,L3
KSSHSVLYSSNNKDFFA,L1

```

### 2. Example csv file without chain:

```
sequence
SDRES
QAWDSSTAWV
KSSHSVLYSSNNKDFFA

```

### 3 Example csv file with sequence:

```
sequence
EVQLVQPGAELRNSGASVKVSCKASGYRFTSYYIDWVRQAPGQGLEWMGRIDPEDGGTKYAQKFQGRVTFTADTSTSTAYVELSSLRSEDTAVYYCARNEWETVVVGDLMYEYEYWGQGTQVTVSSASTKGPSVFPLAPALGCLVKDYFPEPVTVSGVHTFPAVLQSSGLYSLSSVVNVNHK
DIQMTQSPSSLSASLGDRVTITCQASQSISSYLAWYQQKPGQAPNILIYGASRLKTGVPSRFSGSGSGTSFTLTISGLEAEDAGTYYCQQYASVPVTFGQGTKVELKRTVAAPSVFIFPPSVVCLLNNFYPREAKVQWQESVTEQDSKDSTYSLSSTCEVTHQGLSSPVTKSF

```

## Contacts
manasanartem4@gmail.com
