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

    m = np.genfromtxt(args.matrix, delimiter=",", dtype=float)
    coords = np.genfromtxt(args.coords,
                           delimiter=",",
                           skip_header=1,
                           dtype=float)

    if coords.shape[0] == 3:
        coords = np.append(coords, 0)

    affine = np.zeros((4, 4), dtype=float)

    if args.invert:
        affine[:3, :3] = np.linalg.pinv(m)
    else:
        affine[:3, :3] = m

    res = affine @ coords[:, np.newaxis]
    np.savetxt(args.output,
               res.reshape((1, 4)),
               header="x,y,z,t",
               comments="",
               fmt="%.10f",
               delimiter=",")


if __name__ == "__main__":
    main()
