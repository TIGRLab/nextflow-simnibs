'''
Apply orientation transformation on a set of coordinates
'''

import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Apply orientation transform on a set of coordinates")

    parser.add_argument("coords",
                        type=str,
                        help="Comma separated ANTS 3D coordinates file")
    parser.add_argument("matrix", type=str, help="Path to orientation matrix")
    parser.add_argument("output", type=str, help="Output file")
    parser.add_argument("--invert", required=False, action="store_true")
    args = parser.parse_args()

    m = np.genfromtxt(args.matrix, sep=",", dtype=float)
    coords = np.genfromtxt(args.coords, sep=",", skip_header=1, dtype=float)

    affine = np.zeros((4, 4), dtype=float)

    if args.invert:
        affine[:3, :3] = np.lingalg.pinv(m)
    else:
        affine[:3, :3] = m

    res = m @ coords[:, np.newaxis]
    np.savetxt(args.output, res, header="x,y,z,t")


if __name__ == "__main__":
    main()
