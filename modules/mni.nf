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

process _antsWarpInfoMatrix{
    /*
    * Transform input coordinates to match coordinate convention of warp
    * Arguments:
    *   subject (String): Subject key
    *   warpFile (Path): Warpfile
    *
    * Outputs:
    *   subject (String): Subject key
    *   warpFile (Path): Warpfile
    *   transform (Path): Coordinate transform used in warpfile
    */

    input:
    tuple val(subject), path(warpFile)

    output:
    tuple val(subject), path('transform.csv'), emit: transform

    shell:
    '''
    antsTransformInfo !{warpFile} | grep -m 1 -A 3 "Direction" | sed 1d | \
        tr ' ' ',' > transform.csv
    '''
}

process _prepareCoordsForWarp{

    label 'numpy'
    label 'bin'
    /*
    * Transform coordinates with a rotation matrix
    * Arguments:
    *   subject (String): Subject key
    *   matrix (Path): Path to rotation matrix
    *   x (Float): X coordinates
    *   y (Float): Y coordinates
    *   z (Float): Z coordinates
    * Outputs:
    *   coords: (subject, Path coords): transformed coordinates
    */

    input:
    tuple val(subject), path(matrix), val(x), val(y), val(z)

    output:
    tuple val(subject), path("${subject}_fixed_coords.csv"), emit: coords

    shell:
    '''
    #!/bin/bash
    printf "x,y,z,t\n!{x},!{y},!{z},0" > coords.txt
    python /scripts/mni/transformCoords.py \
        coords.txt \
        !{matrix} \
        !{subject}_fixed_coords.csv
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

process _untransformCoordinates{
    /*
    * Perform inverse transform of warp
    * Arguments:
    *   subject (String): subject key
    *   matrix (Path): Path to transformation matrix file
    *   coordinates (Path): Path to ants compatible coordinates file
    *
    * Outputs:
    *   fixedCoordinates: (subject, Path fixed) warped coordinates with corrected orientation
    */

    label 'numpy'
    label 'bin'

    input:
    tuple val(subject), path(matrix), path(coordinates)

    output:
    tuple val(subject), path("${subject}_warpedFixedCoordinates.csv"), emit: fixedCoordinates

    shell:
    '''
    #!/bin/bash
    python /scripts/mni/transformCoords.py \
        !{coordinates} \
        !{matrix} \
        !{subject}_warpedFixedCoordinates.csv \
        --invert
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
    np.save("!{subject}_coords.npy", a[:3])
    '''
}

workflow antsApplyWarpToCoordinates{
/*
* ANTS apply transforms to coordinates workflow
*
* Arguments:
*   coordinates (Channel): [subject, HashMap [x,y,z] coordinates]
*   warps (Channel): [subject, warpFile]
*
* Outputs:
*   warpedCoordinates (Channel): [subject, HashMap [x,y,z] warpedCoordinates]
*/

    take:
        coordinates
        warps

    main:

        // Get warp orientation transform and apply to coordinates
        _antsWarpInfoMatrix(warps)

        _prepareCoordsForWarp(
            _antsWarpInfoMatrix.out.transform
                .combine(coordinates)
                .map { it.flatten() }
        )

        // Apply warp
        _antsApplyWarpToCoordinates(
            warps.join(_prepareCoordsForWarp.out.coords)
        )

        // Revert orientation transform
        _untransformCoordinates(
            _antsWarpInfoMatrix.out.transform
                .join(_antsApplyWarpToCoordinates.out.warpedCoordinates)
        )


    emit:
        warpedCoordinates = _untransformCoordinates.out.fixedCoordinates
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

    emit:
        warped = antsRegistration.out.warped
        warp = antsRegistration.out.warp
        inverseWarp = antsRegistration.out.inverseWarp
}
