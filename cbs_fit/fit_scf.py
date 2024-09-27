#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def exp_model(n, cbs_energy, b, c):
    energy = cbs_energy - b * np.exp(- c * n)
    return energy


def initial_guess(data_n, data_e):
    r"""
    The exponential fit takes three parameters. With three data points
    the equations can be solved for the unknown parameters. Solving the set
    of equations with the last three data points generates the intial guess.

    The energy dependence ansatz
    E(n) = E _\infty - b e ^{-c n}
    """
    if len(data_n) != len(data_e):
        print("Initial guess failed: data of different lengths.")
        exit(1)

    if len(data_n) < 3:
        print("Initial guess failed: less than three entries availabe to fit.")
        exit(1)

    e3 = data_e[-1]
    e2 = data_e[-2]
    e1 = data_e[-3]
    n3 = data_n[-1]
    n2 = data_n[-2]
    n1 = data_n[-3]

    # this works one if the ns are in  steps of 1
    if n3-n2 != 1 or n2-n1 != 1:
        print("Error the SCF fit requires data set in zeta steps of 1.")

    c = np.log((e2 - e1) / (e3 - e2))
    b = (e2 - e3) / (np.exp(- n3 * c) - np.exp(- n2 * c))
    cbs_energy = e3 + b * np.exp(- n3 * c)

    return cbs_energy, b, c


def print_exp_model_paramerters(parameters, header):
    print(f"{header} values:")
    print(f"  E cbs = {parameters[0]:.6f}")
    print(f"  B = {parameters[1]:.6f}")
    print(f"  c = {parameters[2]:.6f}")
    return


def get_SCF_CBS_value(data_n, data_e):
    guess = initial_guess(data_n, data_e)
    fit_parameters, covariance = curve_fit(exp_model, data_n, data_e, p0=guess)
    return fit_parameters[0]


def fit_scf_to_exp_model(data_n, data_e, use_best: bool = False):
    r"""
    If use_best is set to True, the exact fit to the last three points will be
    used.

    Fit SCF energies to the CBS limit assuming that the energies follow
    E(n) = E _\infty - b e ^{-c n}
    """
    guess = initial_guess(data_n, data_e)
    msg = "Initial guess (extrapolation of the last three points)"
    print_exp_model_paramerters(guess, msg)

    if use_best is True:
        return guess

    fit_parameters, covariance = curve_fit(exp_model, data_n, data_e, p0=guess)
    print_exp_model_paramerters(fit_parameters, "Fitting result")
    return fit_parameters


def show_SCF_fitting_result(data_n, data_e, fit_parameters,
                            basis_str: str = ""):
    plt.figure(figsize=(4, 3))
    plt.plot(data_n, data_e, 'o', label="ab inito SCF", zorder=10)
    first = 3
    last = 7
    npts = 100
    step = (last - first) / npts
    expanded_n = [first + i * step for i in range(npts)]
    expanded_n = np.array(expanded_n)
    plt.plot(expanded_n, exp_model(expanded_n, *fit_parameters),
             label=r"$E _\infty - b * e ^{- c * n}$")
    plt.axhline(fit_parameters[0],
                label=r"$E _\infty = $" + f"{fit_parameters[0]:.6f}")
    bottom, top = plt.ylim(bottom=fit_parameters[0] - 0.001)
    ticks = [int(i) for i in range(first, last+1)]
    plt.xticks(ticks=ticks, labels=[str(t) for t in ticks])
    plt.legend()
    if basis_str != "":
        plt.xlabel("n, " + basis_str)
    plt.title("SCF convergence")
    plt.tight_layout()
    plt.show()


def main():
    pwCVnZ = [3, 4, 5]
    pwCVnZ_SCF = [-224.34218395, -224.35967856, -224.36413383]

    data_n = np.asarray(pwCVnZ)
    data_e = np.asarray(pwCVnZ_SCF)
    fit_parameters = fit_scf_to_exp_model(data_n, data_e)
    show_SCF_fitting_result(data_n, data_e, fit_parameters)

    return 0


if __name__ == "__main__":
    main()
