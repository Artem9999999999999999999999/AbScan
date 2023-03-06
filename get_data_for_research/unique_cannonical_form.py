import json

# Загружаем данные из файла
with open("cannonical_form_of_all_cdr_sequences.json", "r") as f:
    data = json.load(f)

# Создаем словарь для хранения количества CDR для каждой канонической формы
cdr_count = {}

# Проходимся по ключам в словаре data
for canonical_form in data:
    # Получаем количество CDR для текущей канонической формы
    num_cdrs = len(data[canonical_form])
    # Добавляем запись в словарь cdr_count
    cdr_count[canonical_form] = num_cdrs

# Выводим результаты
print("Количество CDR для каждой канонической формы:")
for canonical_form in cdr_count:
    print(f"{canonical_form}: {cdr_count[canonical_form]}")

# Записываем результаты в файл
with open("cdr_count.txt", "w") as f:
    f.write(f"Всего канонических форм: {len(list(data.keys()))}\nКоличество CDR для каждой канонической формы:\n")
    for canonical_form in cdr_count:
        f.write(f"{canonical_form}: {cdr_count[canonical_form]}\n")
