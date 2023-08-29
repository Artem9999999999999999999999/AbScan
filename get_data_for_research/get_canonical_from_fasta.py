import json
import subprocess
import os

output_folder = "/home/artem9/Рабочий стол/NIR/data"
# Задаем аргументы для вызова SCALOP
args = ["SCALOP", "-i", "all_sequences_from_sabdab.fasta", "--numbering_scheme", "chothia", "--cdr_definition", "chothia"]

# Вызываем SCALOP и получаем результаты
try:
    output = subprocess.run(args, capture_output=True, text=True, check=True)
except subprocess.CalledProcessError as e:
    print(f"An error occurred: {e}")

# Создаем пустой словарь
json_data = {}

# Проходимся по строкам в output.stdout
for line in output.stdout.split("\n")[1:]:
    fields = line.strip().split()
    if len(fields) == 5:
        # Извлекаем необходимые данные из fields
        input_id = fields[0]
        cdr_id = fields[1]
        cdr_sequence = fields[2]
        canonical_id = fields[3]
        median_id = fields[4]

        # Создаем новый элемент для словаря
        new_entry = {
            "pdb": input_id,
            "cdr": cdr_id,
            "cdr_sequence": cdr_sequence,
            "median": median_id
        }

        # Добавляем новый элемент в словарь
        if canonical_id not in json_data:
            json_data[canonical_id] = []
        json_data[canonical_id].append(new_entry)

# Записываем результат
with open(os.path.join(output_folder,"cannonical_form_of_all_cdr_sequences_and_mediany6.json"), "w") as f:
    json.dump(json_data, f, indent=4)


