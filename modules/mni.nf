/*
* Transformations into MNI space
*/


process antsRegistration {

    /*
    * ANTS settings from
    * https://github.com/nipreps/niworkflows/blob/master/niworkflows/data/t1w-mni_registration_precise_000.json
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
    */

    label 'ants'

    input:
    tuple val(sub), path(t1), path(mni)


    output:
    tuple val(sub), path("${sub}_Warped.nii.gz"), emit: warped
    tuple val(sub), path("${sub}Composite.h5"), emit: warp
    tuple val(sub), path("${sub}InverseComposite.h5"), emit: inverseWarp


    shell:
    '''
    antsRegistration --collapse-output-transforms 1 --dimensionality 3 \
        --initialize-transforms-per-stage 0 --interpolation LanczosWindowedSinc \
        --output [ !{sub}, !{sub}_Warped.nii.gz ] \
        --transform Rigid[ 0.05 ] \
        --metric \
            Mattes[ !{mni}, !{t1}, 1, 56, Regular, 0.25 ] \
        --convergence [ 100x100, 1e-06, 20 ] \
        --smoothing-sigmas 2.0x1.0vox --shrink-factors 2x1 \
        --use-estimate-learning-rate-once 1 --use-histogram-matching 1 \
        --transform Affine[ 0.08 ] \
        --metric \
            Mattes[ !{mni}, !{t1}, 1, 56, Regular, 0.25 ] \
        --convergence [ 100x100, 1e-06, 20 ] \
        --smoothing-sigmas 1.0x0.0vox --shrink-factors 2x1 \
        --use-estimate-learning-rate-once 1 --use-histogram-matching 1 \
        --transform SyN[ 0.1, 3.0, 0.0 ] \
        --metric CC[ !{mni}, !{t1}, 1, 4, None, 1 ] \
        --convergence [ 100x70x50x20, 1e-06, 10 ] \
        --smoothing-sigmas 3.0x2.0x1.0x0.0vox --shrink-factors 8x4x2x1 \
        --use-estimate-learning-rate-once 1 --use-histogram-matching 1 \
        --winsorize-image-intensities [ 0.005, 0.995 ]  --write-composite-transform 1 \
        -v
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
        antsRegistration(
            freesurfer
                .map { s, f -> [s, "${f}/mri/brainmask.mgz"] }
                .combine(mni)
        )

        antsRegistrationQC(
            antsRegistration.out.warped.combine(mni)
        )

    emit:
        warped = antsRegistration.out.warped
        warp = antsRegistration.out.warp
        inverseWarp = antsRegistration.out.inverseWarp
        qcImage = antsRegistrationQC.out.qcImage
}
