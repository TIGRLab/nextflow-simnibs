import argparse
from collections import namedtuple
import simnibs.msh.mesh_io as mesh_io
import numpy as np
from numpy import arcsin, arctan2, deg2rad, cos, sin

Surf = namedtuple("Surf", ["mesh", "coords", "triangles"])
Line = namedtuple("Line", ["start", "end"])


def make_surf(msh, tag):
    m = msh.crop_mesh(tags=tag)
    return Surf(m, m.nodes.node_coord, m.elm.node_number_list[:, :-1] - 1)


def gen_qc(head, brain, outfile, lines=None):

    import meshplot as mp
    mp.website()

    p = mp.plot(head.coords, head.triangles, c=np.array([0.7, 0.7, 0.7]))

    p._Viewer__objects[0]['material'].transparent = True
    p._Viewer__objects[0]['material'].opacity = 0.5
    p._Viewer__objects[0]['material'].metalness = 0.1
    p.add_mesh(brain.coords, brain.triangles)

    if lines is not None:
        starts = np.array([l.start for l in lines])
        ends = np.array([l.end for l in lines])
        p.add_lines(starts, ends, shading={"line_width": 20})

    p.save(outfile)


def get_matsimnibs(centre, n, twist):

    # Compute euler angles
    alpha = arctan2(-n[1], n[2])
    beta = arcsin(n[0])
    gamma = deg2rad(twist)

    # Define rotation matrices for XYZ variant of
    # tait-bryan rotation matrices
    R_ap = np.array([[1, 0, 0], [0, cos(alpha), -sin(alpha)],
                     [0, sin(alpha), cos(alpha)]])

    R_lr = np.array([[cos(beta), 0, sin(beta)], [0, 1, 0],
                     [-sin(beta), 0, cos(beta)]])

    R_tw = np.array([[cos(gamma), -sin(gamma), 0], [sin(gamma),
                                                    cos(gamma), 0], [0, 0, 1]])

    # Get full rotation matrix of coil onto scalp plane
    R = R_ap @ R_lr @ R_tw

    # Invert normal and X to maintain right-handedness
    R[:3, 1] = -R[:3, 1]
    R[:3, 2] = -R[:3, 2]

    # Construct matsimnibs
    msn = np.zeros((4, 4))
    msn[:3, :3] = R
    msn[:3, 3] = centre

    return msn


def main():
    parser = argparse.ArgumentParser(
        description="Place coil on scalp "
        "closest to given coordinate and orient coil. "
        "Returns SimNIBS matsimnibs matrix")

    parser.add_argument("coord", type=str, help="3D coordinate")
    parser.add_argument("mesh", type=str, help="SimNIBS head mesh .msh")
    parser.add_argument("twist",
                        type=float,
                        help="Twist angle of coil. "
                        "Corresponds to BrainSight twist angle")
    parser.add_argument("outfile",
                        type=str,
                        help="matsimnibs .npy output file name")
    parser.add_argument("--qc-file",
                        type=str,
                        help="Generate QC image to verify coil placement")

    args = parser.parse_args()

    msh = mesh_io.read_msh(args.mesh)
    coord = np.load(args.coord)

    head = make_surf(msh, 1005)
    brain = make_surf(msh, 1002)
    head_coord, head_ind = head.mesh.nodes.find_closest_node(coord,
                                                             return_index=True)
    head_ind = head_ind - 1

    # Get centre at coil
    coil_z = head.mesh.nodes_normals().value[head_ind, :]
    msn = get_matsimnibs(head_coord, coil_z, args.twist)
    np.save(args.outfile, msn)

    if args.qc_file:
        project_line = Line(coord, head_coord)
        coil_y = Line(head_coord, head_coord + msn[:3, 1].flatten() * 10)
        gen_qc(head, brain, args.qc_file, lines=[project_line, coil_y])


if __name__ == '__main__':
    main()
