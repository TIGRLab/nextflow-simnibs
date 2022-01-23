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

    shell:
    '''
    mris_convert -c !{freesurfer_file} !{fs_white} !{sub}.!{hemi}.shape.gii
    '''
}


// To define
process freesurfer_prep{
    // We'll use the connectome workbench container for all processes with
    // label 'workbench'
    label 'workbench'
    /* Prepare surfaces for resampling
    *
    Inputs:
    *   subject_ID (str)
    *   hemi (str): [lh, rh]: Freesurfer Hemisphere
    * fs_white (path): Path to freesurfer white matter surface
    * fs_pial (path): Path to freesurfer pial surface
    * fs_sphere (path): path to freesurfer current sphere
    * resample_sphere (path): path to new HCP sphere (32k or 164k) to resample to
    * hcp_hemi (str): [L, R]: HCP Hemisphere
    * Surface vertices (str): [32k, 164k]: total number of greyordinate vertices per hemi

    * Outputs: (subject_ID: str, hemi: [lh, rh], gifti_current_midthick_file: path )
               (subject_ID: str, Hemi: [L, R], hcp_vertice: [32k,164k], gifti_new_midthick_file: path )
               (subject_ID: str, hemi: [lh, rh], gifti_sphere_file: path)  

    */
    input:
    tuple val(sub), val(hemi), path(fs_white), path(fs_pial), path(fs_sphere), path(resample_sphere)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.midthickness.surf.gii"), emit: current_midthickness
    tuple val(sub), val(hemi), path("${sub}.*.midthickness.*_fs_LR.surf.gii"), emit: new_midthickness
    tuple val(sub), val(hemi), path("${sub}.${hemi}.sphere.reg.surf.gii", emit: current_gifti_sphere)

    shell:
    '''
    hcp_hemi=$(echo !{hemi} | cut -c 1 | tr "[:lower:]" "[:upper:]"
    hcp_vertices=$(echo !{resample_sphere} | grep -oE "[0-9]+k")

    wb_shortcuts -freesurfer-resample-prep !{fs_white} !{fs_pial} !{fs_sphere} !{resample_sphere} \
                                           !{sub}.!{hemi}.midthickness.surf.gii \
                                           !{sub}.!{hcp_hemi}.midthickness.!{hcp_vertices}_fs_LR.surf.gii \
                                           !{sub}.!{hemi}.sphere.reg.surf.gii
    '''
}    

process metric_resample{

    label 'workbench'

    input:
    tuple val(sub), val(hemi), path(gifti_file), path(current_gifti_sphere_file),\
     path(resample_sphere), path(current_midthickness_file), path(new_midthickness_file)

    output:
    tupe val(sub), val(hemi), path("${sub}.E.norm.*.*_fs_LR.func.gii"), emit:e_norm_resampled

    shell:
    '''
    hcp_hemi=$(echo !{hemi} | cut -c 1 | tr "[:lower:]" "[:upper:]"
    hcp_vertices=$(echo !{resample_sphere} | grep -oE "[0-9]+k")

    wb_command -metric-resample !{gifti_file} !{current_sphere_file} !{resample_sphere} \
    ADAP_BARY_AREA !{sub}.E.norm.!{hcp_hemi}.!{hcp_vertices}_fs.LR.func.gii \
    -area-surf !{current_midthickness_file} !{new_midthickness_file}
    '''  
}

process create_dense_scalar{
    label 'workbench'

    input:
    tuple val(sub), path(resampled_L_file), path(resampled_R_file)

    output:
    tuple val(sub), path("${sub}.norm.E.${hcp_vertices}_fs_LR.dscalar.nii"), emit:sim_dscalar

    shell:
    '''
    wb_command -cifti-create-dense-scalar !{sub}.norm.E.${hcp_vertices}_fs_LR.dscalar.nii \
    -left-metric !{resampled_L_file} -right-metric !{resample_R_file}
    '''
}

// Make a process for each step in the simnibs_2_cifti_surface script

workflow simnibs2cifti{
    /*
    * Maps SimNIBS
    * Arguments:
    *   simfiles (channel): (subject_id, simulation_file) - simnibs simulation files
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
	
	split_hemi = fs_dir.multimap { sub, fs -> left: [sub, lh, ${fs}/surf/lh.white ],
						  right: [sub, rh, ${fs}/surf/rh.white ]
				}
	split_hemi.right.mix(split_hemi.left)

        i_fs_to_gifti = simfiles.join(split_hemi)
                                .map{sub, simfile, fs ->
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

        freesurfer_prep(
                    fs_dir.map{sub, hemi, fs -> [
                                                    sub, hemi,
                                                    "${fs}/surf/${hemi}.white",
                                                    "${fs}/surf/${hemi}.pial",
                                                    "{fs}/surf/${hemi}.sphere"
                                                 ]}.combine(atlas_dir)
        )                        

        // fs_to_gifti.out.gifti 
        // freesurfer_prep.out.current_gifti_sphere
        // freesurfer_prep.out.current_midthickness
        // freesurfer_prep.out.new_midthickness

        i_metric_resample = fs_to_gifti.out.gifti.join(freesurfer_prep.out.current_gifti_sphere, by:[0,1])
                                                 .combine(atlas_dir).map { sub, hemi, gifti, sphere, atlas -> [ sub, hemi, gifti, sphere, "${atlas}/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii", $atlas/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii]}
                                                 .join(freesurfer_prep.out.current_midthickness)
                                                 .join(freesurfer_prep.out.new_midthickness)

        metric_resample(i_metric_resample)
	i_metric_resample = metric_resample.out.e_norm_resampled.groupTuple(by: 0)
        create_dense_scalar(metric_resample.out.e_norm_resampled)

    emit:
	sim_dscalar = create_dense_scalar.out.sim_dscalar    


    // The final output to return from this workflow
}
