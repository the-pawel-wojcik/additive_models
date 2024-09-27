#!/usr/bin/env python3

import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def cube_decay_model(n, cbs_energy, b):
    energy = cbs_energy - b / n**3
    return energy


def initial_guess(data_n, data_e):
    if len(data_n) != len(data_e):
        print("Initial guess failed: data of different lengths.",
              file=sys.stderr)
        exit(1)

    if len(data_n) < 2:
        print("Initial guess failed: less than two entries availabe to fit.",
              file=sys.stderr)
        exit(1)

    e1 = data_e[-1]
    e2 = data_e[-2]
    n1 = data_n[-1]
    n2 = data_n[-2]

    b = (e1 - e2) / (1/n2**3 - 1 / n1**3)
    cbs_energy = e1 + b / n1 ** 3

    return cbs_energy, b


def print_model_paramerters(parameters, header):
    print(f"{header} values:")
    print(f"  E cbs = {parameters[0]:.6f}")
    print(f"  B = {parameters[1]:.6f}")
    return


def get_correlation_CBS(data_n, data_e):
    guess = initial_guess(data_n, data_e)
    fit_parameters, covariance = curve_fit(
        cube_decay_model, data_n, data_e, p0=guess)
    return fit_parameters[0]


def fit_to_cubic_model(data_n, data_e, use_best: bool = False):
    r"""
    Returns (E _infty, b)
    from
        E = E _\infty - b / n ** 3
    """
    guess = initial_guess(data_n, data_e)
    print_model_paramerters(guess, "Initial guess")
    if use_best is True:
        return guess

    fit_parameters, covariance = curve_fit(
        cube_decay_model, data_n, data_e, p0=guess)
    print_model_paramerters(fit_parameters, "Fitting result")
    return fit_parameters


def show_fit_results(data_n, data_e, fit_parameters, title_extra: str = "",
                     basis_str: str = ""):
    plt.figure(figsize=(4, 3))
    plt.plot(data_n, data_e, 'o', label="ab initio correlation", zorder=10)
    npts = 100
    first = 3
    last = 7
    step = (last - first) / npts
    expanded_n = [first + i * step for i in range(npts)]
    expanded_n = np.array(expanded_n)
    plt.plot(expanded_n, cube_decay_model(expanded_n, *fit_parameters),
             label=r"$E _\infty - b / n ^3$")
    plt.axhline(fit_parameters[0],
                label=r"$E _\infty = $" + f"{fit_parameters[0]:.6f} a.u.")
    ticks = [int(i) for i in range(first, last+1)]
    plt.xticks(ticks=ticks, labels=[str(t) for t in ticks])
    plt.legend()
    plt.ylabel("correlation energy, a.u.")
    xlabel = "n"
    if basis_str != "":
        xlabel += ", " + basis_str
    plt.xlabel(xlabel)
    if title_extra != "":
        plt.title(title_extra)
    plt.tight_layout()
    plt.show()


def main():
    pwCVnZ = [3, 4]
    pwCVnZ_CCSDT_correlation = [-0.96318413, -1.02456663]

    data_n = np.asarray(pwCVnZ)
    data_e = np.asarray(pwCVnZ_CCSDT_correlation)
    fit_parameters = fit_to_cubic_model(data_n, data_e)
    show_fit_results(data_n, data_e, fit_parameters)

    return 0


if __name__ == "__main__":
    main()
