#!/usr/bin/env python3

import argparse
import json

ha2eV = 27.211386245988
eV2cm = 8065.543937
ha2cm = ha2eV * eV2cm


def prepare_xsim_input(cbs):
    """
    Prepares a list of EOM states. Each state has the form:
        ``` python
        {'irrep': {'energy #': int(), 'name': str()},
         'model': str(),
         'energy': {'transition': {'eV': float(), 'au': float(), } } }
        ```
    """
    out_pack = []
    for eom_state in cbs['EOM']:
        eom_ex_cbs_au = eom_state['correlation'] - cbs['cc_correlation']
        out_state = {
            'irrep': eom_state['irrep'],
            'model': eom_state['model'],
            'energy': {
                'transition': {
                    'eV': eom_ex_cbs_au * ha2eV,
                    'au': eom_ex_cbs_au,
                }
            }
        }
        out_pack += [out_state]
    return out_pack


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('cbs')
    args = parser.parse_args()

    with open(args.cbs) as cbs_json:
        cbs = json.load(cbs_json)

    out_pack = prepare_xsim_input(cbs)
    print(json.dumps(out_pack))


if __name__ == "__main__":
    main()
