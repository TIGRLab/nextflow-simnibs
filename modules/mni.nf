/*
* Transformations into MNI space
*/


process antsRegistration {

    /*
    * ANTS settings from 
    * https://github.com/nipreps/niworkflows/blob/master/niworkflows/data/t1w-mni_registration_precise_000.json
    */

    input:
    tuple val(sub), path(t1), path(mni)


    output:
    tuple val(sub), path("${sub}_Warped.nii.gz"), emit: warped
    tuple val(sub), path("${sub}_Composite.h5"), emit: warp
    tuple val(sub), path("${sub}_InverseComposite.h5"), emit: inverseWarp


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
        --winsorize-image-intensities [ 0.005, 0.995 ]  --write-composite-transform 1
    '''
    
}


workflow fs_to_mni {
    /*
    * Register Freesurfer outputs to MNI standard space
    *
    * Arguments:
    *   freesurfer (channel): (sub, fs_dir) Freesurfer directories
    *   mni (str): MNI template (brain)
    *
    * Outputs:
        warped (channel): (sub, warped) Freesurfer T1 warped into MNI
        warp (channel): (sub, warp) ANTS .h5 warp file
        inverseWarp (channel): (sub, inverseWarp) ANTS .h5 inverse warp file
    */
    take:
        freesurfer
        mni


    main:
        antsRegistration(
            freesurfer
                .map { s, f -> [s, "${f}/mri/brainmask.gz"] }
                .join(mni)
        )

    output:
        warped = antsRegistration.out.warped
        warp = antsRegistration.out.warp
        inverseWarp = antsRegistration.out.inverseWarp
}
