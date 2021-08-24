import argparse
import dendropy
import os
import sys


def prepare_for_stag(ifil, odir):
    """

    Parameters
    ----------
    ifil : string
           name of input file (one newick string per line)
    odir :
    """
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
                    species.add(node.taxon.label.split()[0])
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
    prepare_for_stag(args.input, base[0])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, help="Input file", required=True)

    main(parser.parse_args())
