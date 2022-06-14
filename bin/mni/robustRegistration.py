import argparse
from pathlib import Path

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
import nipype.interfaces.io as nio
from nipype.interfaces.image import Reorient
from nipype.interfaces.utility import Select
from niworkflows.anat.ants import init_brain_extraction_wf
from niworkflows.interfaces.norm import SpatialNormalization


def _directory(string):
    directory_path = Path(string)
    if directory_path.is_dir():
        return str(directory_path.absolute())
    else:
        raise NotADirectoryError


def _file(string):
    file_path = Path(string)
    if file_path.parent.exists():
        return str(file_path.absolute())
    else:
        raise OSError(f"Parent directory {file_path.parent} does not exist!")


def init_skullstrip_wf():
    """
    Perform Niworkflows brain extraction on a single T1 image
    """

    wf = pe.Workflow(name='skullstrip_t1')
    inputnode = pe.Node(niu.IdentityInterface(fields=['t1_img']),
                        name='inputnode')

    be_wf = init_brain_extraction_wf(normalization_quality='precise')

    select_t1 = pe.Node(Select(index=[0]), name='select_t1')
    select_mask = pe.Node(Select(index=[0]), name='select_mask')

    reorient_t1 = pe.Node(Reorient(orientation='RAS'), name='reorient_t1')
    reorient_mask = pe.Node(Reorient(orientation='RAS'), name='reorient_mask')

    outputnode = pe.Node(niu.IdentityInterface(fields=["stripped_t1", "mask"]),
                         name="outputnode")

    wf.connect([
        (inputnode, be_wf, [('t1_img', 'inputnode.in_files')]),
        (be_wf, select_t1, [("outputnode.out_file", "inlist")]),
        (be_wf, select_mask, [("outputnode.out_mask", "inlist")]),
        (select_t1, reorient_t1, [("out", "in_file")]),
        (select_mask, reorient_mask, [("out", "in_file")]),
        (reorient_t1, outputnode, [('out_file', 'stripped_t1')]),
        (reorient_mask, outputnode, [('out_file', 'mask')]),
    ])
    return wf


def main():

    parser = argparse.ArgumentParser(
        description="Perform Robust skullstrip and"
        " spatial normalization of T1 images to a reference image")

    parser.add_argument('t1', type=_file, help='Path to T1 file')
    parser.add_argument('reference',
                        type=_file,
                        help='Path to reference image')
    parser.add_argument('out_dir',
                        type=_directory,
                        help="Path to output directory")
    parser.add_argument('out_prefix',
                        type=str,
                        help="Prefix to use for all outputs")
    parser.add_argument('--nthreads',
                        type=int,
                        help='Number of threads to use')
    parser.add_argument('--work',
                        type=_directory,
                        help="Nipype working directory")

    args = parser.parse_args()

    wf = pe.Workflow(name="robust_registration")
    wf.base_dir = args.work

    skullstrip_wf = init_skullstrip_wf()
    skullstrip_wf.inputs.inputnode.t1_img = args.t1

    spatial_norm = pe.Node(SpatialNormalization(reference_image=args.reference,
                                                flavor='precise',
                                                num_threads=args.nthreads),
                           name='spatial_norm',
                           num_threads=args.nthreads)
    datasink = pe.Node(nio.DataSink(), name='data_sink')
    datasink.inputs.base_directory = args.out_dir
    datasink.inputs.substitutions = [
        ("ants_t1_to_mni", f"{args.out_prefix}_"), ("_Warped", "Warped"),
        ("09_relabel_wm_mask_xform_ras", f"{args.out_prefix}_mask"),
        (f"{args.out_prefix}_corrected_xform_masked_ras",
         f"{args.out_prefix}_brain")
    ]

    wf.connect([(skullstrip_wf, spatial_norm, [('outputnode.stripped_t1',
                                                'moving_image')]),
                (spatial_norm, datasink, [('warped_image', '@warped_image'),
                                          ('composite_transform',
                                           '@composite_transform'),
                                          ('inverse_composite_transform',
                                           '@inverse_composite_transform')]),
                (skullstrip_wf, datasink, [('outputnode.stripped_t1',
                                            '@skullstripped'),
                                           ('outputnode.mask', '@mask')])])

    wf.run()


if __name__ == '__main__':
    main()
