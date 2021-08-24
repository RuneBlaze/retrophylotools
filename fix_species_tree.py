"""
Copyright 2019 Erin K. Molloy

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
 
3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
"""
import argparse
import dendropy
import os
import sys


def is_binary(tree):
    """
    Checks if a tree is binary

    Parameters
    ----------
    tree : dendropy tree object

    Returns
    -------
    True if the tree is binary and False otherwise
    """
    nodes = [n for n in tree.preorder_node_iter()]
    for n in nodes[1:]:
        if not n.is_leaf():
            children = n.child_nodes()
            if len(children) != 2:
                return False
    return True


def is_rooted(tree):
    """
    Checks if a tree is rooted

    Parameters
    ----------
    tree : dendropy tree object

    Returns
    -------
    True if the tree is rooted and False otherwise
    """
    n = tree.seed_node
    children = n.child_nodes()
    if len(children) != 2:
        return False
    return True


def scale_branch_lengths(tree, constant):
    """
    Rescales all branch lengths in tree by a constant factor

    Parameters
    ----------
    tree : dendropy tree object
    constant : float

    Returns
    -------
    Nothing
    """
    constant = float(constant)
    for e in tree.preorder_edge_iter():
        if e.length is not None:
            e.length = e.length * constant


def force_ultrametric(tree):
    """
    Forces a tree to have ultrametric branch lengths

    Parameters
    ----------
    tree : dendropy tree object

    Returns
    -------
    Nothing

    Example
    -------
    Suppose that (A:y,(B:x,C:x):z); is *not* ultrametric, that is,

        y + H(A) != z + H(A,B)

    where H(A) is the height node A (which equals 0) and
          H(B,C) is the height of node (B,C) (which equals x).

    Then, taking y + z = l, we adjust the branch lengths as follows:

        y' = (l - H(A) + H(B,C))/2
        z' = l - y'

    to force the resulting tree (A:y',(B:x,C:x):z'); to be ultrametric.
    """
    if not is_binary(tree):
        sys.exit("Tree is not binary!")
    if not is_rooted(tree):
        sys.exit("Tree is not binary!")

    for n in tree.postorder_node_iter():
        if n.is_leaf():
            n.height = 0.0
        else:
            children = n.child_nodes()
            c1 = children[0]
            c2 = children[1]

            h1 = c1.edge.length + c1.height
            if h1 != c2.edge.length + c2.height:
                l = c1.edge.length + c2.edge.length

                e1 = (l + c2.height - c1.height) / 2.0
                e2 = l - e1

                c1.edge.length = e1
                c2.edge.length = e2

                h1 = c1.edge.length + c1.height
                if h1 != c2.edge.length + c2.height:
                    sys.exit("Unable to force tree to be ultrametric!")

            n.height = h1
        last = n.height

    sys.stdout.write("Species Tree Height: %f\n" % last)


def main(args):
    tree = dendropy.Tree.get(path=args.input, schema="newick")
    scale_branch_lengths(tree, args.factor)
    force_ultrametric(tree)
    tree.write(path=args.output, schema="newick")

    os._exit(0)  # CRITICAL ON BLUE WATERS LOGIN NODE


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fixes species tree.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input file")
    parser.add_argument(
        "-f", "--factor", type=str, required=True, help="Branch length scale factor"
    )
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file")
    main(parser.parse_args())
