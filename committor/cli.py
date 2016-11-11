#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK                                -*- mode: python -*-

import argparse
try:
    import argcomplete
except ImportError:
    argcomplete = None
import os
import sys

import committor


def auto_complete(p, ext):
    """Enable auto-completion of file name extension."""
    if not argcomplete:
        return

    p.completer = argcomplete.completers.FilesCompleter(allowednames=ext)


def main():
    parser = argparse.ArgumentParser(description='committor function solver')

    p = parser.add_argument('-g', '--grad-file', metavar='GRAD-FILE', type=str,
                            required=True, help='path to file containing '
                            'the gradient')
    auto_complete(p, 'grad')

    p = parser.add_argument('-p', '--pmf-file', metavar='PMF-FILE', type=str,
                            required=False, help='path to file containing '
                            'the potential of mean force')
    auto_complete(p, 'pmf')

    p = parser.add_argument('-e', '--est-file', metavar='EST-FILE', type=str,
                            required=False, help='path to file containing '
                            'the estimated gradient of the potential of mean '
                            'force')
    auto_complete(p, 'est')

    p = parser.add_argument('-b', '--boundaries', metavar='CFG-FILE', type=str,
                            required=True, help='definition of the reactant and '
                            'the product subdomains')
    auto_complete(p, 'cfg')

    p = parser.add_argument('-o', '--output', metavar='PVD-FILE', type=str,
                            required=True, help='path to file where the '
                            'resulting committor function will be saved')
    auto_complete(p, 'pvd')

    parser.add_argument('-n', '--num-cells', metavar='NUM-CELLS', type=int,
                        required=True, help='number of cells in each '
                        'dimension')

    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress information')

    if argcomplete:
        argcomplete.autocomplete(parser)

    args = parser.parse_args(sys.argv[1:])

    if args.pmf_file is None:
        args.pmf_file = os.path.extsep.join([args.grad_file, 'pmf'])
    if args.est_file is None:
        args.est_file = os.path.extsep.join([args.grad_file, 'est'])

    committor.committor(args.grad_file, args.pmf_file, args.est_file,
                        args.output, args.num_cells,
                        args.boundaries, args.verbose)


if __name__ == '__main__':
    main()
