# Репозиторий для инструмента идентификации структурно-критичных аминокислот в CDR1 и CDR2 антител

## Описание проекта
Цель проекта - разработать инструмент для идентификации структурно-критичных аминокислот в CDR1 и CDR2 антител. 

## Структура репозитория
data/ - папка с исходными данным

data/median_structure_cdr/ - информация о всех CDR центроид

data/pairwise_align _canonical_form_and_median_structure/ - в этой папке лежат все парные выравнивания извлеченных из sabdab последовательностей CDR на CDR центроиды. Каждый fasta- файл соответствует одной канонической форме и ее структуре-центроиде

data/freq_aa_in_canonical_form/ - здесь лежит json-файл со всем частотами аминокислот для различных семейств и позиций


/Файл "sabdab_summary_all_tsv" из базы данных "sabdab" содержит актуальную информацию на 6 марта. Файл sabdab_data_pdb_id.csv был получен путем удаления всех колонок, кроме pdb_id и информации о тяжелой и легкой цепях из "sabdab_summary_all_tsv", далее был конвертирован в CSV. Подробную информацию о структуре "sabdab_summary_all_tsv" можно посмотреть по адресу https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/.

get_data_for_research/ - папка с кодом для извлечения данных: sabdab(scv) => fasta => json | txt

get_freq.py - Данный скрипт позволяет получить информацию о канонической форме и частотах аминокислот, присутствующих в последовательности антитела. Для этого используется файл с данными о частотах аминокислот в канонических формах.

Также предоставляется возможность получить данные о частотах аминокислот для определенной последовательности CDR (complementarity-determining region) и соответствующей ей канонической формы. Для этого используется локальная версия SCALOP, которая должна быть установлена на компьютере и добавлена в переменную PATH.

После запуска предлагается выбрать один из двух режимов работы - 1 или 2. При выборе первого режима, пользователь должен ввести последовательность аминокислот, для которой необходимо получить каноническую форму и частоты аминокислот. При выборе второго режима, пользователь должен ввести каноническую форму и последовательность CDR, для которой необходимо получить частоты аминокислот.

После выполнения на экран будет выведена информация о канонической форме, частотах аминокислот в последовательности CDR, а также информация о частотах аминокислот для каждой позиции CDR. Если для некоторой позиции CDR информация о частотах отсутствует, программа выведет список аминокислот, которые встречаются в этой позиции, а также их частоты в других канонических формах.



README.md - файл с описанием проекта

## Контакты
manasanartem4@gmail.com
