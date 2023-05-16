import argparse
from freq_in_file import get_freq_in_file
from freq_in_console import get_freq_from_scalop, get_freq_from_cdr, get_freq_from_scalop_with_chain
from mutant import get_all_mutant
from mutant_with_aligment import get_all_mutant_with_aligment


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
    parser.add_argument('-o', '--output',
                        help='File name/path for writing results',
                        type=str,
                        default='result.json')

    args = parser.parse_args()

    if not args.sequence and not args.cdr:
        parser.print_help()
        return

    if args.sequence[-4:]:
        if args.sequence[-4:] != '.csv':
            seq_for_scalop = args.sequence
            number = args.number
            get_freq_from_scalop(seq_for_scalop, number)

    if args.cdr and not args.family:
        if args.cdr[-4:] != '.csv':
            seq = args.cdr
            number = args.number
            get_freq_from_cdr(seq, number_alter=number)

    if args.cdr and args.family:
        if args.cdr[-4:] != '.csv':
            seq = args.cdr
            number = args.number
            chain = args.family
            get_freq_from_scalop_with_chain(seq, chain, number_alter=number)

    if args.cdr:
        if args.cdr[-4:] == '.csv':
            input_filename = args.cdr
            output_filename = args.output
            get_freq_in_file(input_filename, output_filename)

    if args.mutant and args.cdr and args.family:
        seq = args.cdr
        output_filename = args.output
        chain = args.family
        get_all_mutant(seq, chain, output_filename)

    if args.mutant and args.cdr and not args.family:
        seq = args.cdr
        output_filename = args.output
        get_all_mutant_with_aligment(seq, output_filename)


if __name__ == "__main__":
    main()
