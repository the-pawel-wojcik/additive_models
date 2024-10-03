#!/usr/bin/env python3

import argparse
import json
from turn_cbs_into_xsim_input import prepare_xsim_input
from find_cbs import basis2n

ha2eV = 27.211386245988
eV2cm = 8065.543937
ha2cm = ha2eV * eV2cm


def get_args():
    parser = argparse.ArgumentParser(
        description="Prints out the CBS energies for each EOM state "
        "together with half of the CBS–last-ab-initio difference.")

    parser.add_argument(
        'ab_initio', help="JSON file with ab initio energies for each "
        "basis set.")
    parser.add_argument('cbs', help="JSON file with CBS energies.")
    parser.add_argument(
        '-x', '--xsim',
        help="Print output as a `better energies` input for xsim.",
        action='store_true', default=False)
    parser.add_argument(
        '-s', '--summary',
        help="Print output to standard output.",
        action='store_true', default=False)

    args = parser.parse_args()
    return args


def flip_data_to_xsim_like(dataset):
    """
    Takes as an input a dataset that lists the EOM-CC calculation results
    for various basis sets. The input should be in the same format as the one
    which is accepted by the CBS limit calculations.

    Returns the data converted to look like the CBS limit output, i.e.,
    xsim-like 'better energies' input, except that it is a list of them for
    every basis set.
    ```
    [{
            'basis': str(),
            'cc_energy': float(),
            'EOM': [{
                'irrep': {'energy #': int(), 'name': str()},
                'model': str(),
                'energy': {'transition': {'au': float(), 'eV': float()}},
                },]
    },]
    ```
    """

    basis_data = list()
    for data in dataset:
        if 'cc_energy' not in data:
            raise RuntimeError(
                "Some entries are missing the 'cc_energy' value. Most likely"
                " you used the SCF at a basis where CC is not possible. Remove"
                " the extra point from the ab initio json file."
            )
        cc_energy = data['cc_energy']
        eom_states = list()
        for state in data['EOM']:
            e_total_au = state['energy']
            e_eom_au = e_total_au - cc_energy
            ready_state = {
                'irrep': state['irrep'],
                'model': state['model'],
                'energy': {
                    'transition': {
                        'au': e_eom_au,
                        'eV': e_eom_au * ha2eV,
                    },
                },
            }
            eom_states += [ready_state]
        outpack = {
            'basis': data['basis'],
            'cc_energy': cc_energy,
            'EOM': eom_states,
        }
        basis_data += [outpack]

    return basis_data


def states_match(foo, bar):
    """
    Compares if the two dictionaries talk about the same state.
    """
    if foo['irrep'] != bar['irrep']:
        return False

    if foo['model'] != bar['model']:
        return False

    return True


def main():
    args = get_args()
    with open(args.ab_initio) as ai_json:
        dataset = json.load(ai_json)

    with open(args.cbs) as cbs_json:
        cbs = json.load(cbs_json)

    cbs_final = prepare_xsim_input(cbs)
    basis_data = flip_data_to_xsim_like(dataset)
    basis_data.sort(key=lambda x: basis2n[x['basis']])
    best_ab_initio = basis_data[-1]['EOM']

    better_energies = []
    string_out = ""
    for ai_state in best_ab_initio:
        for cbs_state in cbs_final:
            if not states_match(ai_state, cbs_state):
                continue
            name = str(cbs_state['irrep']['energy #'])
            name += cbs_state['irrep']['name']
            ecbs = cbs_state['energy']['transition']['eV']
            eai = ai_state['energy']['transition']['eV']
            error_est = 0.5 * abs(eai - ecbs)
            string_out += f"{name:4} = {ecbs:6.3f} ± {error_est:5.3f} eV\n"
            state = {
                "irrep": cbs_state['irrep'],
                "eom model": cbs_state['model'],
                "energy": {
                    "transition": {
                        "eV": ecbs,
                        "au": ecbs / ha2eV,
                    },
                },
            }
            better_energies += [state]

    if args.summary is True:
        print(string_out)

    if args.xsim is True:
        print(json.dumps(better_energies))

    return 0


if __name__ == "__main__":
    main()
