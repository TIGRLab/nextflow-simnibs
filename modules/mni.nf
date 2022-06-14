/*
* Transformations into MNI space
*/


process convert_to_nifti {
    /*
    * Convert Freesurfer MGZ file to NIFTI
    *
    *
    * Arguments:
    *   sub (str): Subject ID
    *   t1 (path): Path to Freesurfer MGZ file
    *
    * Outputs:
    *   
    *   nifti (queue): (sub, nifti) Converted NIFTI file
    */

    label 'freesurfer'

    input:
    tuple val(sub), path(mgz)

    output:
    tuple val(sub), path("${sub}.nii.gz"), emit: nifti

    shell:
    '''
    mri_convert !{mgz} !{sub}.nii.gz
    '''

}

process antsRegistration {

    /*
    * Run antsBrainExtraction and antsRegistration-based Nipype workflow
    *
    * Arguments:
    *   sub (str): Subject ID
    *   t1 (path): Path to 3D moving image
    *   mni (path): Path to fixed 3D image (registration target)
    *
    * Outputs:
    *   warped (queue): (sub, warped) Warped moving image
    *   warp (queue): (sub, warp) Forward warp field
    *   inverseWarp (queue): (sub, warp) Inverse warp field
    *   brain (queue): (sub, brain) Skullstrip brain output
    *   mask (queue): (sub, mask) Skullstrip mask output
    */

    label 'niworkflows'
    label 'bin'

    input:
    tuple val(sub), path(t1), path(mni)


    output:
    tuple val(sub), path("${sub}_Warped.nii.gz"), emit: warped
    tuple val(sub), path("${sub}_Composite.h5"), emit: warp
    tuple val(sub), path("${sub}_InverseComposite.h5"), emit: inverseWarp
    tuple val(sub), path("${sub}_brain.nii.gz"), emit: brain
    tuple val(sub), path("${sub}_mask.nii.gz"), emit: mask


    shell:
    '''
    mkdir work
    python /scripts/mni/robustRegistration.py \
        !{t1} !{mni} \
        "." "!{sub}" \
        --nthreads !{task.cpus} \
        --work ./work
    '''

}

process antsRegistrationQC{
/* Generate QC image for ANTs registration
*
* Arguments:
*   subject (String): Subject identifier key
*   moving (Path): moving image
*   fixed (Path): fixed image
*
* Output:
*   qcImage (Queue): [subject,
*                     Path qcImage] : Path to QC image]
*/

    label 'niviz'

    input:
    tuple val(subject), path(moving), path(fixed)

    output:
    tuple val(subject), path("${subject}_qc-registration.svg"), emit: qcImage

    shell:
    '''
    niviz single \
        registration \
        --set bg_nii=!{fixed} \
        --set fg_nii=!{moving} \
        !{subject}_qc-registration.svg
    '''
}

process _prepareCoordsForWarp{

    input:
    tuple val(x), val(y), val(z)

    output:
    path("lps_coords.csv"), emit: coordsCsv

    shell:
    '''
    echo "x,y,z,t" > lps_coords.csv
    echo "!{x}, !{y}, !{z}, 0" >> lps_coords.csv
    '''
}

process _antsApplyWarpToCoordinates{

    /*
    * Perform antsApplyWarpToCoordinates
    * Arguments:
    *   subject (String): Subject key
    *   warpFile (Path): Path to warp file
    *   coordinates (Path): Path to ants compatible coordinates.csv
    *
    * Outputs:
    *   warpedCoordinates: (subject, Path warped) coordinates warped into target space
    */

    label 'ants'

    input:
    tuple val(subject), path(warpFile), path(coordinates)

    output:
    tuple val(subject), path("${subject}_warpedCoordinates.csv"), emit: warpedCoordinates

    shell:
    '''
    antsApplyTransformsToPoints \
        -d 3 \
        -p 1 \
        -i !{coordinates} \
        -o !{subject}_warpedCoordinates.csv \
        -t !{warpFile}
    '''
}

process antsToNumpy{
    /*
    * Transform ANTS style CSV into numpy array
    */

    label 'numpy'

    input:
    tuple val(subject), path(antsCoords)

    output:
    tuple val(subject), path("${subject}_coords.npy"), emit: coords

    shell:
    '''
    #!/usr/bin/env python

    import numpy as np

    a = np.loadtxt("!{antsCoords}", skiprows=1, delimiter=",")

    # LPS --> RAS
    a[:2] = -a[:2]

    # Remove t
    np.save("!{subject}_coords.npy", a[:3])
    '''
}

workflow antsApplyWarpToCoordinates{
/*
* ANTS apply transforms to coordinates workflow
*
* Arguments:
*   coordinates (Channel): [subject, [x,y,z] coordinates]
*   warps (Channel): [subject, warpFile]
*
* Outputs:
*   warpedCoordinates (Channel): [subject, Path: warpedCoordinates]
*/

    take:
        coordinates
        warps

    main:
        lps_coordinates = coordinates.map { x, y, z -> [-x, -y, z] } 
        lps_coordinates | view
        _prepareCoordsForWarp(lps_coordinates)
        _antsApplyWarpToCoordinates(warps.combine(_prepareCoordsForWarp.out.coordsCsv))

    emit:
        warpedCoordinates = _antsApplyWarpToCoordinates.out.warpedCoordinates

}


workflow registerFreesurferToMNI {
    /*
    * Register Freesurfer outputs to MNI standard space
    *
    * Arguments:
    *   freesurfer (channel): (sub, fs_dir) Freesurfer directories
    *   mni (str): MNI template (brain)
    *
    * Outputs:
    *   warped (channel): (sub, warped) Freesurfer T1 warped into MNI
    *   warp (channel): (sub, warp) ANTS .h5 warp file
    *   inverseWarp (channel): (sub, inverseWarp) ANTS .h5 inverse warp file
    */
    take:
        freesurfer
        mni


    main:

        convert_to_nifti(freesurfer.map { s, fsdir -> [
            s, "${fsdir}/mri/T1.mgz"
        ]})

        antsRegistration(convert_to_nifti.out.nifti.combine(mni))

        antsRegistrationQC(
            antsRegistration.out.warped.combine(mni)
        )

    emit:
        warped = antsRegistration.out.warped
        warp = antsRegistration.out.warp
        inverseWarp = antsRegistration.out.inverseWarp
        qcImage = antsRegistrationQC.out.qcImage
        brain = antsRegistration.out.brain
        mask = antsRegistration.out.mask
}
