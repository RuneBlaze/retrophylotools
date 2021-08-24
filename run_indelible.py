"""
Creates INDELible control file (using parameters from set_indelible_params.py)
and then runs INDELible; see INDELible tutorial pages here:
http://abacus.gene.ucl.ac.uk/software/indelible/tutorial/NUCLEOTIDE.shtml
http://abacus.gene.ucl.ac.uk/software/indelible/tutorial/rates.shtml

Written by EKM (molloy.erin.k@gmail.com) in Spring 2018.
Edited by EKM in Summer 2019.
"""
import argparse
from decimal import Decimal
import dendropy
import numpy
import pandas
import os
import sys


def branch_lengths_2_decimals(str_newick_tree):
    """
    Replaces branch lengths in scientific notation with decimals taken
    github.com/ngcrawford/CloudForest/blob/master/cloudforest/phybase.py

    COPY OF LICENSE FILE:

    Copyright (c) 2011, Nick G. Crawford & Brant C. Faircloth

    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    + Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
    + Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
    + Neither the name of the University of California, Boston University, nor
    the names of its contributors may be used to endorse or promote products
    derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
    """
    colon_s = 0
    comma_back_paren_s = 0
    num = ""
    new_tree = ""
    for count, char in enumerate(str_newick_tree):
        if char == ":":
            colon_s = count
            continue
        if char in (")", ","):
            comma_back_paren_s = 1
            num = "%1.12f" % Decimal(num)
            new_tree += ":" + num
            colon_s = 0
            num = ""
        if colon_s != 0:
            num = num + char
        if colon_s == 0:
            new_tree += char
    new_tree = new_tree.strip("'").strip('"').strip("'")
    return new_tree


def cleanup(tree):
    for n in tree.preorder_node_iter():
        if not n.is_leaf():
            n.label = None


def run_indelible(indelible, params, trees, outdir, tmpdir):
    """
    Parameters
    ----------
    indelible : INDELible executable (may need to include full path)
    params : pandas dataframe
             Each row has INDELible simulation parameters
                + GENE: Gene ID
                + fT, fC, fA, fG: GTR base frequencies
                + rAC, rAG, rAT, CG, CT, GT: GTR transition rate matrix
                + ALPH: Gamma distribution (across sites rate heterogeneity)
                + SQLN: Sequence length
    trees : list of strings
            Each string is a newick string for a rooted gene tree; trees[i]
            corresponds to Gene ID i (from params)
    outdir : string
             Output directory (may need to include full path)
    tmpdir : string
             Temporary directory (may need to include full path)

    Returns
    -------
    Nothing, writes output files
    """
    os.makedirs(tmpdir)
    os.chdir(tmpdir)

    ipath, iexec = indelible.rsplit("/", 1)
    os.system("cp " + indelible + " .")

    pad = len(str(len(trees)))
    for i, r in params.iterrows():
        gene = str(int(r["GENE"])).zfill(pad)

        tree = dendropy.Tree.get(
            string=trees[i],
            schema="newick",
            rooting="force-rooted",
            preserve_underscores=True,
        )
        cleanup(tree)
        tree = tree.as_string(schema="newick").replace("'", "")
        tree = branch_lengths_2_decimals(tree)[5:]

        with open("control.txt", "w") as f:
            f.write("[TYPE] NUCLEOTIDE 1\n")
            f.write("[MODEL] modelname\n")
            # [MODEL] GTR CT AT GT AC CG [divide by AG, i.e., AG = 1]
            f.write(
                "[submodel] GTR %f %f %f %f %f\n"
                % (
                    r["rCT"] / r["rAG"],
                    r["rAT"] / r["rAG"],
                    r["rGT"] / r["rAG"],
                    r["rAC"] / r["rAG"],
                    r["rCG"] / r["rAG"],
                )
            )
            # [statefreq] T C A G
            f.write("[statefreq] %f %f %f %f\n" % (r["fT"], r["fC"], r["fA"], r["fG"]))
            f.write("[rates] 0 %f 0 \n" % r["ALPH"])
            f.write("[TREE] treename  " + tree + "\n")
            f.write("[PARTITIONS] partitionname\n")
            f.write("[treename modelname %d]\n" % r["SQLN"])
            f.write("[EVOLVE] partitionname 1 %s\n" % gene)
        os.system("./" + iexec)
        os.rename(gene + "_TRUE.phy", gene + ".phy")
        os.remove(gene + ".fas")

    os.remove(iexec)
    os.remove("trees.txt")
    os.remove("control.txt")
    os.remove("LOG.txt")
    os.system("mv " + tmpdir + "/*.phy " + outdir)
    os.system("rm -rf " + tmpdir)


def main(args):
    with open(args.trees, "r") as f:
        trees = [l for l in f]

    params = pandas.read_csv(args.params)
    params = params[(params["GENE"] >= args.start) & (params["GENE"] <= args.end)]

    tmpdir = args.output + "/" + "tmp-" + str(args.start) + "-" + str(args.end)

    if not os.path.exists(tmpdir):
        run_indelible(args.indelible, params, trees, args.output, tmpdir)

    os._exit(0)  # CRITICAL ON BLUE WATERS LOGIN NODE


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-x",
        "--indelible",
        type=str,
        help="Path to INDELible executable",
        required=True,
    )
    parser.add_argument("-s", "--start", type=int, help="Start index", required=True)
    parser.add_argument("-e", "--end", type=int, help="End index", required=True)
    parser.add_argument(
        "-p",
        "--params",
        type=str,
        help="Parameter list file (one row per tree)",
        required=True,
    )
    parser.add_argument("-t", "--trees", type=str, help="Tree list file", required=True)
    parser.add_argument(
        "-o", "--output", type=str, help="Output directory", required=True
    )

    main(parser.parse_args())
