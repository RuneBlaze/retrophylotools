import argparse
import dendropy
from dendropy.calculate.treecompare import false_positives_and_negatives
from itertools import izip


def compare_trees(tr1, tr2):
    """
    Compares two trees

    Parameters
    ----------
    tr1 : dendropy tree object
            First tree (typically the model tree)
    tr2 : dendropy tree object
            Second tree (typically the estimated tree)

    Returns
    -------
    nl : int
         Size of the shared leaf set, i.e., the number of leaves in both trees
    ei1 : int
          Number of internal edges in first tree (after restricting it to the shared leaf set)
    ei2 : int
          Number of internal edges in second tree (after restricting it to the shared leaf set)
    fn : int
         Number of edges in the first tree that are not in the second tree
    fp : int
         Number of edges in the second tree that are not in the first tree
    rf : float
         Normalized Robinson-Foulds (RF) distance between the first and second trees

    Example
    -------
    If tree 1 corresponds to "(((A,B,C),D),E);" and tree 2 corresponds to "((((A,B),C),D),E);",
    then the output is "5 1 2 0 1 0.25". In this example,
      + first and second trees share 5 leaves (A, B, C, D, E).
      + first tree has one internal edge "A,B,C|D,E"
      + second tree has two internal edges "A,B|C,D,E" and "A,B,C|D,E"
      + one edges in the first tree that are missing from the second tree
      + no edge "A,B|C,D,E" in the second tree that is missing in the first tree
      + normalized RF distance is (FP+FN)/(2*NL-6) = (1+0)/(2*5-6) = 0.25
    """

    # Unroot the two trees!
    tr1.is_rooted = False
    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    tr2.is_rooted = False
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    # Restrict the two trees to the same leaf set if necessary!
    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])

    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)

        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)

        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)

    # Compare trees!
    tr1.update_bipartitions()
    tr2.update_bipartitions()

    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))

    if nl < 4:
        return (nl, ei1, ei2, "NA", "NA", "NA")

    [fn, fp] = false_positives_and_negatives(tr1, tr2)
    rf = (fn + fp) / (2.0 * nl - 6.0)

    return (nl, ei1, ei2, fn, fp, rf)


def main(args):
    if args.prefix is None:
        p = ""
    else:
        p = str(args.prefix + ",")

    with open(args.output, "aw") as fo, open(args.index, "r") as fi, open(
        args.treelist1, "r"
    ) as f1, open(args.treelist2, "r") as f2:

        for li, l1, l2 in izip(fi, f1, f2):
            i = int(li)

            taxa = dendropy.TaxonNamespace()

            tre1 = dendropy.Tree.get(
                string=l1,
                schema="newick",
                rooting="force-unrooted",
                taxon_namespace=taxa,
            )

            tre2 = dendropy.Tree.get(
                string=l2,
                schema="newick",
                rooting="force-unrooted",
                taxon_namespace=taxa,
            )

            [nl, ei1, ei2, fn, fp, rf] = compare_trees(tre1, tre2)
            if rf == "NA":
                fo.write("%s%d,%d,%d,%d,%s,%s,%s\n" % (p, i, nl, ei1, ei2, fn, fp, rf))
            else:
                fo.write(
                    "%s%d,%d,%d,%d,%d,%d,%1.6f\n" % (p, i, nl, ei1, ei2, fn, fp, rf)
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-t1", "--treelist1", type=str, help="Input tree list 1 file", required=True
    )
    parser.add_argument(
        "-t2", "--treelist2", type=str, help="Input tree list 2 file", required=True
    )
    parser.add_argument(
        "-i", "--index", type=str, help="Index list file", required=True
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        help="Append prefix to each row of CSV",
        required=False,
    )
    parser.add_argument(
        "-o", "--output", type=str, help="Output CSV file", required=True
    )

    main(parser.parse_args())
