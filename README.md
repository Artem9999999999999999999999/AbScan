# Репозиторий для инструмента идентификации структурно-критичных аминокислот в CDR1 и CDR2 антител

## Описание проекта
Цель проекта - разработать инструмент для идентификации структурно-критичных аминокислот в CDR1 и CDR2 антител. 

## Структура репозитория
data/ - папка с исходными данным
data/median_structure_cdr/ - информация о всех cdr центроид
data/pairwise_align _canonical_form_and_median_structure/ - в этой папке лежат все парные выравнивания извлеченных из sabdab последовательностей CDR на CDR центроиды. Каждый fasta- файл соответствует одной канонической форме и ее структуре-центроиде
data/freq_aa_in_canonical_form/ - здесь лежит json-файл со всем частотами аминокислот для различных семейств и позиций


/Файл "sabdab_summary_all_tsv" из базы данных "sabdab" содержит актуальную информацию на 6 марта. Файл sabdab_data_pdb_id.csv был получен путем удаления всех колонок, кроме pdb_id и информации о тяжелой и легкой цепях из "sabdab_summary_all_tsv", далее был конвертирован в CSV. Подробную информацию о структуре "sabdab_summary_all_tsv" можно посмотреть по адресу https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/.

get_data_for_research/ - папка с кодом для извлечения данных: sabdab(scv) => fasta => json | txt

get_freq.py - скрипт, который на вход получает семейство и последовательность аминокислот cdr, а на выходе дает информацию о всех частотах, если информации о какой-то аминокислоте нет в наших данных, то будут предложена информация о всех известных АК для данной позиции

# Пример работы:

input: 

Введите семейство и последовательность аминокислот CDR через пробел: H1-7-A GYTFTCW

output:

Частоты аминокислот:

В позиции 1 аминокислота G имеет частоту 97.45%

В позиции 2 аминокислота Y имеет частоту 22.05%

В позиции 3 аминокислота T имеет частоту 47.67%

В позиции 4 аминокислота F имеет частоту 66.51%

В позиции 5 аминокислота T имеет частоту 21.1%

В позиции 6 аминокислота C имеет частоту 0.01%

В позиции 7 аминокислота W имеет частоту 0.09%


README.md - файл с описанием проекта

## Контакты
manasanartem4@gmail.com
