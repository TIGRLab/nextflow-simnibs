nextflow.preview.dsl = 2


// Example process
process fs_to_gifti{

    // Label allows us to define multi-process settings
    // I.e use a Freesurfer container for all processes with
    // label 'freesurfer'
    label 'freesurfer'

    /*
    * Convert Freesurfer file to GIFTI
    *
    * Inputs:
    *   subject_ID (str)
    *   hemi (str): [lh, rh]: Freesurfer Hemisphere
    *   fs_white (path): Path to freesurfer white matter surface
    *   freesurfer_file: Path to surface file to be converted
    *
    * Output: (subject_ID: str, hemi: [lh, rh], gifti_file: path)
    */
    input:
    tuple val(sub), val(hemi), path(fs_white), path(freesurfer_file)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.shape.gii"), emit: gifti

    '''
    mris_convert -c !{freesurfer_file} !{fs_white} !{sub}.!{hemi}.shape.gii
    '''
}


// To define
process freesurfer_prep{
    // We'll use the connectome workbench container for all processes with
    // label 'workbench'
    label 'workbench'
}

process metric_resample{
}

process create_dense_scalar{
}

// Make a process for each step in the simnibs_2_cifti_surface script

workflow simnibs2cifti{
    /*
    * Maps SimNIBS
    * Arguments:
    *   simfiles (channel): (subject_id, hemi: [lh,rh], simulation_file) - simnibs simulation files
    *   fs_dir (channel): (subject_id, fs_dir) - simnibs freesurfer directories
    *   cifti_dir (channel): (subject_id, ciftify_dir) - ciftify directories
    *   atlas_dir (value): path - single value path to standard mesh atlases
    *
    * Returns:
    *   (channel): (subject_id, field_dscalar) - simulation file resampled to CIFTI space
    */
    take:
        simfiles
        fs_dir
        cifti_dir
        atlas_dir

    main:

        // join works like a table join. It'll link up channels using a key (subject_ID)
        // map is like python/R's `map` function where you can apply a function
        // to each item in a channel
        // each item in the joined channel is a tuple of
        // (subject, hemisphere, simulation_file, freesurfer_directory)
        i_fs_to_gifti = simfiles.join(fs_dir)
                                .map{sub, hemi, simfile, fs ->
                                    [
                                       sub, hemi, simfile,
                                       "${fs}/surf/${hemi}.white"
                                    ]
                                }

        fs_to_gifti(i_fs_to_gifti)

        // fs_to_gifti.out.gifti <- method to access output variable of process
        // you can pipe (|) variables into functions like `view`, which will display
        // the channel
        fs_to_gifti.out.gifti | view


    // The final output to return from this workflow
    output:
        sim_dscalar
}
