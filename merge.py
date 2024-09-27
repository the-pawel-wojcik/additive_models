#!/usr/bin/env python

import argparse
import json
import sys


def get_args():
    parser = argparse.ArgumentParser(
        description="Take two files prepared to be inputs for xsim's 'better "
        "energies' file and add them togher: energies and 'model' names.")
    parser.add_argument('first',
                        help="file in the xsim's 'better energies' format.")
    parser.add_argument('second',
                        help="file in the xsim's 'better energies' format.")
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    with open(args.first) as first_json:
        first = json.load(first_json)

    with open(args.second) as second_json:
        second = json.load(second_json)

    for state in first:
        found_it = False
        for addition in second:
            if addition['irrep'] != state['irrep']:
                continue

            energy_corr = addition['energy']['transition']
            state['energy']['transition']['eV'] += energy_corr['eV']
            state['energy']['transition']['au'] += energy_corr['au']

            state['model'] += '+' + addition['model']

            found_it = True
            break

        if found_it is False:
            print("Warning! Second file is missing data about "
                  f"{state['irrep']}", file=sys.stderr)

    print(json.dumps(first))

    return 0


if __name__ == "__main__":
    main()
