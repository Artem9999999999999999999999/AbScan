import argparse
from freq_in_file import get_freq_in_file_for_seq, get_freq_in_file_for_cdr
from freq_in_console import get_freq_from_scalop, get_freq_from_cdr, get_freq_from_scalop_with_chain
from mutant import get_all_mutant


def main():
    parser = argparse.ArgumentParser(
        description='Retrieve the canonical form and/or frequencies of amino acids in an antibody sequence. To work '
                    'with the amino acid sequence, you must have a local version of SCALOP installed on your device '
                    'and a PATH defined. if the chain parameter is not specified, then the canonical family will be '
                    'determined using alignment. Otherwise, the family will be defined with SCALOP')
    parser.add_argument('-s', '--sequence',
                        help='Amino acid sequence for which the canonical form and frequencies will be obtained',
                        type=str,
                        default='d')
    parser.add_argument('-c', '--cdr',
                        help='CDR sequence for which frequency data and its corresponding canonical form will be '
                             'retrieved',
                        type=str)
    parser.add_argument('-n', '--number',
                        help='The number of alternative amino acids to display (default=3)',
                        type=int,
                        default=3)
    parser.add_argument('-f', '--family',
                        help='Path to csv file containing one sequence column or two sequence and circuit, '
                             'the possible options are: L1, L2, L3, H1, H2')
    parser.add_argument('-m', '--mutant',
                        help='Search for mutations that do not affect the canonical form',
                        action='store_true')
    parser.add_argument('-i', '--input',
                        help='File name/path for reading. The file must be in csv format, either with one sequence column, or with two sequence and chaine,'
                             'the possible options for chain are: L1, L2, L3, H1, H2',
                        type=str)
    parser.add_argument('-o', '--output',
                        help='File name/path for writing results',
                        type=str,
                        default='result.txt')
    parser.add_argument('-v', '--variabel',
                        help='modifier, used when the input file contains an antibody chain sequence',
                        action='store_true')

    args = parser.parse_args()

    if not args.sequence and not args.cdr:
        parser.print_help()
        return

    if args.sequence and not args.input:
        seq_for_scalop = args.sequence
        number = args.number
        get_freq_from_scalop(seq_for_scalop, number)

    if args.cdr and not args.family:
        seq = args.cdr
        number = args.number
        get_freq_from_cdr(seq, number_alter=number)

    if args.cdr and args.family:
        seq = args.cdr
        number = args.number
        chain = args.family
        get_freq_from_scalop_with_chain(seq, chain, number_alter=number)

    if args.input and not args.variabel:
        input_filename = args.input
        output_filename = args.output
        get_freq_in_file_for_cdr(input_filename, output_filename[:-3] + 'json')

    if args.input and args.variabel:
        input_filename = args.input
        output_filename = args.output
        get_freq_in_file_for_seq(input_filename, output_filename[:-3] + 'json')

    if args.mutant and args.cdr and args.family:
        seq = args.cdr
        output_filename = args.output
        chain = args.family
        get_all_mutant(seq, chain, output_filename[:])


if __name__ == "__main__":
    main()