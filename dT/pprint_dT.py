#!/usr/bin/env python3

import argparse
import json
import sys

ha2eV = 27.211386245988
eV2cm = 8065.543937
ha2cm = ha2eV * eV2cm


def get_args():
    parser = argparse.ArgumentParser(
        description="Prints the higher CC level correction to the EOM"
        " energies, e.g., ΔT/ANO1 := EOM-CCSDT/ANO1 - EOM-CCSD/ANO1.",
        epilog="The input files can be generated with a combination of the"
        " cfour parser (xcfour.py)"
        " and the cfour processor (print_roots.py -c)."
    )

    parser.add_argument(
        'better',
        help="JSON file with ab initio energies at a higher level of theory."
    )

    parser.add_argument(
        'worse',
        help="JSON file with ab initio energies at a lower level of theory."
    )

    parser.add_argument(
        '-u', '--units', type=str, default="eV",
        help="Choose energy units."
    )

    parser.add_argument(
        '-x', '--xsim', default=False, action="store_true",
        help="Print correction in the xsim's 'better energies' format."
    )

    args = parser.parse_args()
    return args


def states_match(foo, bar):
    """
    Compares if the two dictionaries talk about the same state.
    """
    if foo['irrep'] != bar['irrep']:
        return False

    return True


def main():
    args = get_args()
    with open(args.better) as better_json:
        better = json.load(better_json)

    with open(args.worse) as worse_json:
        worse = json.load(worse_json)

    basis = better['basis']
    if worse['basis'] != basis:
        print(
            "Error! Both calculations need to use the same basis set!",
            file=sys.stderr
        )
        return 1

    better_level = better['calclevel']
    worse_level = worse['calclevel']

    worse_len = len(worse_level)
    if worse_len >= len(better_level):
        print(
            f"Error! The better level of theory, {better_level}, is expeted"
            f" to have a longer name than the worse one, {worse_level}.",
            file=sys.stderr
        )
        return 1

    core = better_level[:worse_len]

    if core != worse_level:
        print(f"Warning! Correcting {worse_level} with the seemingly"
              f" incompatibile {better_level}.", file=sys.stderr)

    better_cc_energy = better['cc_energy']
    worse_cc_energy = worse['cc_energy']

    conversion_factors = {
        "au": 1.0,
        "eV": ha2eV,
        "cm": ha2cm,
    }

    if args.units not in conversion_factors:
        print(f"Error! Unknown units choose from: {conversion_factors.keys()}",
              file=sys.stderr)
        return 1

    correction_name = "Δ" + better_level[worse_len:] + "/" + basis
    conversion = conversion_factors[args.units]

    message = ""
    message += f"The {correction_name} correction.\n"
    message += f"Energies in {args.units}.\n\n"
    message += f"{'State':5} {worse_level:>6} {better_level:>6}"
    message += f" {correction_name} Err. est.\n"

    # container for xsim's output
    better_energies = []
    float_fmt = "6.3f"
    for better_state in better['EOM']:
        state_str = str(better_state['irrep']['energy #'])
        state_str += better_state['irrep']['name']
        found_it = False

        for worse_state in worse['EOM']:
            if states_match(better_state, worse_state):
                worse_eom_energy = worse_state['energy'] - worse_cc_energy
                better_eom_energy = better_state['energy'] - better_cc_energy
                eom_energy_correction = better_eom_energy - worse_eom_energy
                error_est = 0.5 * abs(eom_energy_correction)
                message += f"{state_str:5}"
                message += f" {worse_eom_energy*conversion:{float_fmt}}"
                message += f" {better_eom_energy*conversion:{float_fmt}}"
                message += f" {eom_energy_correction*conversion:{float_fmt}}"
                message += f" {error_est*conversion:{float_fmt}}"
                message += "\n"

                state = {
                    "irrep": better_state['irrep'],
                    "model": correction_name,
                    "energy": {
                        "transition": {
                            "eV": eom_energy_correction *
                            conversion_factors['eV'],
                            "au": eom_energy_correction *
                            conversion_factors['au'],
                        },
                    },
                }
                better_energies += [state]

                found_it = True
                break

        if found_it is False:
            print("Warning! No match in the lower level calculation for the"
                  f" state {state_str}.", file=sys.stderr)

    if args.xsim is True:
        print(json.dumps(better_energies))
    else:
        print(message)

    return 0


if __name__ == "__main__":
    main()
