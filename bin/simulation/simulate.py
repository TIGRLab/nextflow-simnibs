'''
Run SimNIBS simulation
'''

import argparse
import numpy as np
from simnibs.simulation import sim_struct


def main():
    parser = argparse.ArgumentParser(
        description="Run SimNIBS simulation on mesh.")
    parser.add_argument("mesh", type=str, help="Head model .msh")
    parser.add_argument("matsimnibs", type=str, help="Matsimnibs .npy matrix")
    parser.add_argument("coil", type=str, help="Coil definition file")
    parser.add_argument("--gifti",
                        required=False,
                        action="store_true",
                        help="Output GIFTI files")
    parser.add_argument("--m2m-path",
                        required=False,
                        type=str,
                        help="Path to m2m folder, required if using --gifti")

    args = parser.parse_args()

    s = sim_struct.SESSION()
    s.fnamehead = args.mesh
    s.pathfem = "."
    if args.gifti and not args.m2m_path:
        raise AttributeError(
            "GIFTI output required but --m2m-path not supplied")
    else:
        s.map_to_fsavg = args.gifti
        s.map_to_surf = args.gifti
        s.subpath = args.m2m_path

    tmslist = s.add_tmslist()
    tmslist.fnamecoil = args.coil

    pos = tmslist.add_position()
    pos.matsimnibs = np.load(args.matsimnibs)

    s.run()


if __name__ == '__main__':
    main()
