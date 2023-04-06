import json
import subprocess

with open('freq_aa_in_canonical_form.json') as f:
    data = json.load(f)


def get_freq_from_cdr(family, seq):
    seq_len = len(seq)
    try:
        freq = {}
        for i, aa in enumerate(seq):
            amino_found = False
            for d in data[family]:
                if str(i + 1) in d:
                    if aa in d[str(i + 1)]:
                        freq[i + 1] = d[str(i + 1)][aa]
                    else:
                        freq[i + 1] = None
                    break
            else:
                freq[i + 1] = None

        print()
        print("Amino acids frequencies:")
        for k, v in freq.items():
            if v is None:
                for d in data[family]:
                    if str(k) in d:
                        print(
                            f"There is no information about the frequencies of the amino acid {seq[k - 1]} in position {k}, but the following amino acids occur in this position:")
                        for aa_key, freq_value in d[str(k)].items():
                            print(f"{aa_key}: {freq_value}%")
                        break
            else:
                print(f"At position {k} amino acid {seq[k - 1]}: {v}%")

        if seq_len > len(data[family][0]):
            print()
            if seq_len == len(data[family][0]) + 1:
                print(
                    f"There is no information for the amino acid at the position {len(data[family][0]) + 1}, most likely something is wrong with your CDR =(")
            else:
                print(
                    f"There is no information for amino acids in positions from {len(data[family][0]) + 1} to {seq_len}, most likely something is wrong with your CDR =(")

    except KeyError:
        print("Error: key not found.")


def get_freq_from_scalop(seq_for_scalop):
    args = ["SCALOP", "-i", seq_for_scalop, "--numbering_scheme", "chothia", "--cdr_definition", "chothia"]

    try:
        output = subprocess.run(args, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")

    my_dict = {}

    for line in output.stdout.split("\n")[3:]:
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
                "cdr_sequence": cdr_sequence,
            }

            # Добавляем новый элемент в словарь
            if canonical_id not in my_dict:
                my_dict[canonical_id] = []
            my_dict[canonical_id].append(new_entry)

    for key, value in my_dict.items():
        print()
        print(f"Frequencies for CDR '{value[0]['cdr_sequence']}' of the canonical form {key}: ")
        get_freq_from_cdr(key, value[0]['cdr_sequence'])


if __name__ == "__main__":
    print("Greetings, user. How may I assist you today?")
    print("1. Obtain the canonical form and frequencies of amino acids present in the antibody sequence.")
    print("2. Retrieve the frequency data for a specific CDR sequence and its corresponding canonical form.")
    print("CAUTION! Please ensure that the local version of SCALOP is installed on your device and has been added to the PATH variable.")
    ans = input("Please enter 1 or 2: ")
    if ans == '1':
        seq_for_scalop = str(input("Please enter the amino acid sequence: "))
        get_freq_from_scalop(seq_for_scalop)
    elif ans == '2':
        family, seq = input("Please enter the canonical form and CDR sequence, separated by a space: ").split()
        get_freq_from_cdr(family, seq)
    else:
        print("I'm terribly sorry, I'm not able to comprehend your request at the moment. But I promise to learn and improve soon!")



