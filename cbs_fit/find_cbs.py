#!/usr/bin/env python3

import argparse
import json
import numpy as np
import fit_correlation as fc
import fit_scf as fscf

ha2eV = 27.211386245988
eV2cm = 8065.543937
ha2cm = ha2eV * eV2cm

basis2n = {
    'ANO0': 1,
    'ANO1': 2,
    'ANO2': 3,
    'PWCVDZ': 2,
    'PWCVTZ': 3,
    'PWCVQZ': 4,
    'PWCV5Z': 5,
    'aug-pwCVDZ': 2,
    'aug-pwCVTZ': 3,
    'aug-pwCVQZ': 4,
    'aug-pwCV5Z': 5,
}


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'ab_initio', help="JSON file with ab initio energies for each "
        "basis set."
    )
    parser.add_argument(
        '-b', '--basis', default="",
        help="String name for the basis set family."
    )
    parser.add_argument(
        "-a", "--use_all", default=False, action="store_true",
        help="Fit using all data points. Default: Fit using minimal number of"
        " points, i.e., only the largest basis sets."
    )
    args = parser.parse_args()
    return args


def add_correlation_energies(dataset):
    for data in dataset:
        scf = data['scf']

        if 'cc_energy' in data:
            cc = data['cc_energy']
            data['cc_correlation'] = cc - scf

        for eom in data['EOM']:
            ee = eom['energy']
            eom['correlation'] = ee - scf


def get_dataset(file_name):
    with open(file_name) as input:
        data = json.load(input)

    data.sort(key=lambda x: basis2n[x['basis']])
    add_correlation_energies(data)

    cc_level = None
    for basis in data:
        if 'calclevel' in basis:
            loc_cc_level = basis['calclevel']
            if cc_level is None:
                cc_level = loc_cc_level
                continue

            if cc_level != loc_cc_level:
                raise RuntimeError("Data set contains varying CC level: "
                                   f"{cc_level}, and {loc_cc_level}")

    return data, cc_level


def fit_scf(dataset, use_best: bool = False, basis_str: str = ""):
    """
    If use_best is set to True, the fit will use only the three points from
    the largest basis sets.
    """
    print("\n\nFitting SCF energy to the exponential model.")
    zetas = [basis2n[data['basis']] for data in dataset]
    scf_energies = [data['scf'] for data in dataset]
    exp_model_scf_fit_parameters = fscf.fit_scf_to_exp_model(
        zetas, scf_energies, use_best=use_best)
    fscf.show_SCF_fitting_result(
        zetas, scf_energies, exp_model_scf_fit_parameters, basis_str)
    scf_cbs = exp_model_scf_fit_parameters[0]
    return scf_cbs


def fit_cc(dataset, name, key, use_best: bool = False, basis_str: str = ""):
    r"""
    Fit using only the last two points if use_best is set.
    Energy assumption:
        E = E _\infty - b / n ** 3
    """
    print(
        "\n\n"
        f"Fitting {name} correlation energy to the 1/n^3 model."
    )
    zetas = [basis2n[data['basis']] for data in dataset]
    correlation_energies = [data[key] for data in dataset]
    cc_parameters = fc.fit_to_cubic_model(
        zetas, correlation_energies, use_best=use_best)

    fc.show_fit_results(zetas, correlation_energies,
                        cc_parameters, name, basis_str)
    correlation_cbs = cc_parameters[0]
    correlation_cbs_error_est = 0.5 * np.abs(
        correlation_cbs - correlation_energies[-1])

    return correlation_cbs, correlation_cbs_error_est


def fit_eom(dataset, use_best):
    """
    For each data in the EOM section find extrapolated energy.
    """
    print(
        "\n\n"
        "Fitting EOM correlation energy to the 1/n^3 model."
    )
    zetas = [basis2n[data['basis']] for data in dataset]

    eom_model = dataset[0]['EOM'][0]['model']
    eom_cbs = []
    for eom_state in dataset[0]['EOM']:
        eom_cbs += [{
            'irrep': eom_state['irrep'],
            'model': eom_model,
        }]

    for desired_eom_state in eom_cbs:
        # For each state in the EOM list run extrapolation
        correlation_energies = list()

        # Collect the desired_eom_state's energies from every basis size
        for one_basis_data in dataset:
            for eom_state in one_basis_data['EOM']:
                if eom_state['irrep'] != desired_eom_state['irrep']:
                    continue

                if eom_state['model'] != eom_model:
                    raise RuntimeError(f"Warning! Varying EOM models:"
                                       f"{eom_model} and {eom_state['model']}")

                correlation_energies += [eom_state['correlation']]

        cc_parameters = fc.fit_to_cubic_model(
            zetas, correlation_energies, use_best=use_best)

        irrep_no = desired_eom_state['irrep']['energy #']
        irrep_name = desired_eom_state['irrep']['name']
        name = eom_model + f" {irrep_no}{irrep_name}"
        fc.show_fit_results(zetas, correlation_energies, cc_parameters, name)
        correlation_cbs = cc_parameters[0]
        correlation_cbs_error_est = 0.5 * abs(
            correlation_cbs - correlation_energies[-1])

        desired_eom_state['correlation'] = correlation_cbs
        desired_eom_state['correlation error est'] = correlation_cbs_error_est

    return eom_cbs


def main():
    args = get_args()
    basis_str = args.basis
    dataset, cc_level = get_dataset(args.ab_initio)

    use_best = not args.use_all
    scf_cbs = fit_scf(dataset, use_best=use_best, basis_str=basis_str)

    # Allow for a case where for the largest basis sets only SCF is available,
    # e.g. aug-pwCVDZ, aug-pwCVTZ, but for QZ only SCF is possible, and it is
    # needed because SCF if fitted to three points at least

    cc_dataset = [data for data in dataset if 'cc_correlation' in data]

    cc_correlation_CBS, cc_corr_error_est = fit_cc(
        cc_dataset, name=cc_level, key='cc_correlation', use_best=use_best,
        basis_str=basis_str)

    cbs_eom = fit_eom(cc_dataset, use_best=use_best)
    cbs = {
        'basis': 'CBS',
        'scf': scf_cbs,
        'calclevel': cc_level,
        'cc_correlation': cc_correlation_CBS,
        'EOM': cbs_eom,
    }
    print(json.dumps(cbs))


if __name__ == "__main__":
    main()
