import argparse
import dendropy


def create_species_counter(stax, ltax):
    counter = {}
    for x in stax:
        counter[x] = 0

    for x in ltax:
        if x[:5] != "Lost-":
            y = x.split()[0]
            counter[y] = counter[y] + 1

    return counter


def main(args):
    stax = dendropy.TaxonNamespace()
    stre = dendropy.Tree.get(
        path=args.stree, schema="newick", rooting="force-unrooted", taxon_namespace=stax
    )
    stax = sorted([x.label for x in stax])

    with open(args.output, "aw") as fo:
        # Write CSV HEADER
        if args.prefix is None:
            fo.write("GENE,")
        else:
            prefix = str(args.prefix + ",")
            if args.column is None:
                fo.write("PREFIX,GENE,")
            else:
                fo.write(args.column + ",GENE,")
        for x in stax:
            fo.write("NCPY_" + x + ",")
        fo.write("GTRE_NLEA,GTRE_NTAX\n")

        # Compare trees
        with open(args.ltreelist, "r") as fi:
            for l, line in enumerate(fi.readlines()):
                ltax = dendropy.TaxonNamespace()
                ltre = dendropy.Tree.get(
                    string=line,
                    schema="newick",
                    rooting="force-unrooted",
                    taxon_namespace=ltax,
                )
                ltax = [x.label for x in ltax]

                counter = create_species_counter(stax, ltax)

                fo.write("%s%d," % (prefix, l + 1))
                nl = 0
                nx = 0
                for x in stax:
                    cr = counter[x]
                    fo.write("%d," % cr)
                    nl = nl + cr
                    if cr > 0:
                        nx = nx + 1
                fo.write("%d,%d\n" % (nl, nx))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-s", "--stree", type=str, help="Input species tree file", required=True
    )
    parser.add_argument(
        "-l", "--ltreelist", type=str, help="Input gene tree list file", required=True
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        help="Append prefix to each row of CSV",
        required=False,
    )
    parser.add_argument(
        "-c", "--column", type=str, help="Column labels for prefix ", required=False
    )
    parser.add_argument(
        "-o", "--output", type=str, help="Output CSV file", required=True
    )

    main(parser.parse_args())
