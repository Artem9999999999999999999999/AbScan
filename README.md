# Репозиторий инструмента идентификации структурно-критичных аминокислот в CDR1 и CDR2 антител

## Описание проекта
Цель проекта - разработать инструмент для идентификации структурно-критичных аминокислот в CDR1 и CDR2 антител. 

## Структура репозитория
data/ - папка с исходными данным

data/median_structure_cdr/ - информация о всех CDR центроид

data/pairwise_align _canonical_form_and_median_structure/ - в этой папке лежат все парные выравнивания извлеченных из sabdab последовательностей CDR на CDR центроиды. Каждый fasta- файл соответствует одной канонической форме и ее структуре-центроиде

data/freq_aa_in_canonical_form/ - здесь лежит json-файл со всем частотами аминокислот для различных семейств и позиций


/Файл "sabdab_summary_all_tsv" из базы данных "sabdab" содержит актуальную информацию на 6 марта. Файл sabdab_data_pdb_id.csv был получен путем удаления всех колонок, кроме pdb_id и информации о тяжелой и легкой цепях из "sabdab_summary_all_tsv", далее был конвертирован в CSV. Подробную информацию о структуре "sabdab_summary_all_tsv" можно посмотреть по адресу https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/.

get_data_for_research/ - папка с кодом для извлечения данных: sabdab(scv) => fasta => json | txt

## Получаем частоты

AbScan.py - Данный скрипт позволяет получить информацию о канонической форме и частотах аминокислот, присутствующих в последовательности антитела. Для этого используется файлы с данными о частотах аминокислот в канонических формах и локальная версия SCALOP, которая должна быть установлена на компьютере и добавлена в переменную PATH.
```
usage: AbScan.py [-h] [-s SEQUENCE] [-c CDR] [-n NUMBER] [-f FAMILY]
                          [-m] [-o OUTPUT]

Retrieve the canonical form and/or frequencies of amino acids in an antibody
sequence. To work with the amino acid sequence, you must have a local version
of SCALOP installed on your device and a PATH defined. if the chain parameter
is not specified, then the canonical family will be determined using
alignment. Otherwise, the family will be defined with SCALOP

options:
  -h, --help            show this help message and exit
  -s SEQUENCE, --sequence SEQUENCE
                        amino acid sequence for which the canonical form and
                        frequencies will be obtained
  -c CDR, --cdr CDR     CDR sequence for which frequency data and its
                        corresponding canonical form will be retrieved
  -n NUMBER, --number NUMBER
                        the number of alternative amino acids to display
                        (default=3)
  -f FAMILY, --family FAMILY
                        enter the name of the circuit, the possible options
                        are: L1, L2, L3, H1, H2
  -m, --mutant          Search for mutations that do not affect the canonical
                        form
  -o OUTPUT, --output OUTPUT
                        File name for writing results
```

## Контакты
manasanartem4@gmail.com
