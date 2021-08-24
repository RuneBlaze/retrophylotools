"""
Written by Erin Molloy (molloy.erin.k@gmail.com) in Summer 2019.
"""
import argparse
import numpy
import numpy as np
import os
import math
from numpy.random import *


def set_indelible_params(ngens, freqs, rates, dist, args, threshold, sqlen, outf):
    """
    Parameters
    ----------
    ngens : int
            Number of genes
    freqs : list of floats (length 4)
            Alpha values (A, C, G, T) to define a Dirichlet distribution from
            which GTR base frequencies will be drawn.
    rates : list of floats (length 6)
            Alpha values (AC, AG, AT, CG, CT, GT) to define a Dirichlet
            distribution from which GTR transition rates will be drawn.
    alpha : list of floats (length 2)
            Meanlog and sdlog to define a Lognormal distribution from which
            alpha (to define a gamma for site rate heterogeneity) will be drawn.
    sqlen : int
            Sequence length
    outf : string
           output file name

    Returns
    -------
    Nothing, writes an output file
    """
    with open(outf, "w") as f:
        f.write("GENE,SQLN,ALPH,fA,fC,fG,fT,rAC,rAG,rAT,rCG,rCT,rGT\n")
        fmat = numpy.random.dirichlet((freqs[0], freqs[1], freqs[2], freqs[3]), ngens)
        rmat = numpy.random.dirichlet(
            (rates[0], rates[1], rates[2], rates[3], rates[4], rates[5]), ngens
        )
        avec = getattr(numpy.random, dist)(*args, ngens)
        for i in range(avec.shape[0]):
            while avec[i] < threshold:
                avec[i] = getattr(numpy.random, dist)(*args)
        gene = 1
        for j in range(ngens):
            fj = fmat[j, :]
            rj = rmat[j, :]
            f.write(
                "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"
                % (
                    gene,
                    sqlen() if callable(sqlen) else sqlen,
                    avec[j],
                    fj[0],
                    fj[1],
                    fj[2],
                    fj[3],
                    rj[0],
                    rj[1],
                    rj[2],
                    rj[3],
                    rj[4],
                    rj[5],
                )
            )
            gene = gene + 1


def main(args):
    getattr(numpy.random, args.dist)
    set_indelible_params(
        args.ngens,
        args.freqs,
        args.rates,
        args.dist,
        args.args,
        args.threshold,
        eval(args.sqlen.replace("λ", "lambda")),
        args.output,
    )

    os._exit(0)  # CRITICAL ON BLUE WATERS LOGIN NODE


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-n", "--ngens", type=int, required=True, help="Specify number of genes."
    )
    parser.add_argument(
        "-f",
        "--freqs",
        type=float,
        nargs=4,
        required=True,
        help="GTR base frequencies: Specify alpha parameters "
        "(A, C, G, T) for Dirichlet distribution.",
    )
    parser.add_argument(
        "-r",
        "--rates",
        type=float,
        nargs=6,
        required=True,
        help="GTR transition rate matrix: Specify alpha "
        "parameters (AC, AG, AT, CG, CT, GT) for Dirichlet"
        " distribution.",
    )
    #     parser.add_argument("-a", "--alpha", type=float, nargs='+', required=True,
    #                         help="Alpha: specify parameters for distribution")
    parser.add_argument(
        "-d",
        "--dist",
        type=str,
        required=True,
        help="Dist: distribution, " "for examplenumpy.random.[lognormal]",
    )
    parser.add_argument(
        "-a",
        "--args",
        type=float,
        nargs="+",
        required=True,
        help="Args: arguments for the distribution",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=-math.inf,
        help="Treshold: alpha values below will be redrawn",
    )
    parser.add_argument(
        "-l",
        "--sqlen",
        type=str,
        required=True,
        help="Specify sequence length (eval parameter)",
    )
    # for example, for ASTRAL-II
    # (λ lm: (λ sg: (λ: lognormal(lm, sg))))(uniform(5.7,7.3))(uniform(0.0,0.3))
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file")
    main(parser.parse_args())
