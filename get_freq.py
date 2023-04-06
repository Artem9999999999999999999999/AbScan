import json

with open('freq_aa_in_canonical_form.json') as f:
    data = json.load(f)

family, seq = input("Введите семейство и последовательность аминокислот CDR через пробел: ").split()
seq_len = len(seq)

try:
    freq = {}
    for i, aa in enumerate(seq):
        amino_found = False
        for d in data[family]:
            if str(i + 1) in d:
                if aa in d[str(i + 1)]:
                    freq[i + 1] = d[str(i + 1)][aa]
                    amino_found = True
                else:
                    freq[i + 1] = None
                    amino_found = True
                break
        else:
            freq[i + 1] = None

    print()
    print("Частоты аминокислот:")
    for k, v in freq.items():
        if v is None:
            for d in data[family]:
                if str(k) in d:
                    print(f"Информации о частотах аминоксилоты {seq[k - 1]} в позиции {k} нет, но в этой позиции встречаются следующие аминокислоты: ")
                    for aa_key, freq_value in d[str(k)].items():
                        print(f"{aa_key}: {freq_value}%")
                    break
        else:
            print(f"В позиции {k} аминокислота {seq[k - 1]} имеет частоту {v}%")

    if seq_len > len(data[family][0]):
        print()
        if seq_len == len(data[family][0]) + 1:
            print(f"Для аминокислоты в позиции {len(data[family][0]) + 1} информации нет, скорее всего с вашей CDR что-то не так =(")
        else:
            print(f"Для аминокислот в позициях с {len(data[family][0])+1} по {seq_len} информации нет, скорее всего с вашей CDR что-то не так =(")

except KeyError:
    print("Ошибка: ключ не найден.")

