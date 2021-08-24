import argparse
import dendropy
import os
import sys


def read_g2s_map(ifil):
    """
    Parameters
    ----------
    ifil : string
           name of input gene to species label map file (ASTRAL-multi)

    Returns
    -------
    g2sm : python dictionary
           maps gene labels to species labels
    """
    g2sm = {}
    with open(ifil, "r") as f:
        for line in f.readlines():
            [species, genes] = line.split(":")
            genes = genes.split(",")
            genes[-1] = genes[-1].replace("\n", "")
            for gene in genes:
                g2sm[gene] = species
    return g2sm


def prepare_for_stag(ifil, mfil, odir):
    """

    Parameters
    ----------
    ifil : string
           name of input file (one newick string per line)
    odir :
    """
    g2sm = read_g2s_map(mfil)

    # Make output directory
    os.mkdir(odir)
    os.mkdir(odir + "/GeneTrees")

    species = set()
    with open(ifil, "r") as f:
        for l, line in enumerate(f.readlines()):
            temp = "".join(line.split())
            taxa = dendropy.TaxonNamespace()
            tree = dendropy.Tree.get(
                data=temp,
                schema="newick",
                rooting="force-unrooted",
                taxon_namespace=taxa,
            )

            for node in tree.postorder_node_iter():
                if node.is_leaf():
                    species.add(g2sm[node.taxon.label])
                else:
                    # Remove internal node label
                    node.label = None

            # Write gene tree
            with open(odir + "/GeneTrees/" + str(l) + ".txt", "w") as f:
                f.write(tree.as_string(schema="newick")[5:].replace("'", ""))

    # Write gene to species map
    with open(odir + "/SpeciesMap.txt", "w") as f:
        for s in species:
            f.write(s + "_* " + s + "\n")


def main(args):
    base = args.input.rsplit(".", 1)
    prepare_for_stag(args.input, args.map, base[0])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, help="Input file", required=True)
    parser.add_argument("-a", "--map", type=str, help="Input file", required=True)

    main(parser.parse_args())
