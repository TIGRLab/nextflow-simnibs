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
    tuple val(sub), val(hemi), path(gifti_file), path(current_sphere_file),\
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
    tuple val(sub), path("${sub}.norm.E.32k_fs_LR.dscalar.nii"), emit:sim_dscalar

    shell:
    '''
    wb_command -cifti-create-dense-scalar !{sub}.norm.E.32k_fs_LR.dscalar.nii \
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
    //    cifti_dir
        atlas_dir

    main:

        // join works like a table join. It'll link up channels using a key (subject_ID)
        // map is like python/R's `map` function where you can apply a function
        // to each item in a channel
        // each item in the joined channel is a tuple of
        // (subject, hemisphere, freesurfer_directory, simulation_file)
	
    

	    split_hemi = fs_dir.multiMap { sub, fs -> 
                                left: [sub, "lh", "${fs}/surf/lh.white" ]
						        right: [sub, "rh", "${fs}/surf/rh.white" ]
				}
	
        sim_inputs = simfiles.transpose().map { s, sim -> [s, (sim =~ ~/[l,r]h/)[0], sim]}

        i_fs_to_gifti = split_hemi.left.mix(split_hemi.right).join(sim_inputs,by:[0,1])

        fs_to_gifti(i_fs_to_gifti)

        // fs_to_gifti.out.gifti <- method to access output variable of process
        // you can pipe (|) variables into functions like `view`, which will display
        // the channel

        //fs_to_gifti.out.gifti | view

        fs_files = fs_dir.multiMap { sub, fs -> 
                left: [ sub, "lh", "${fs}/surf/lh.white", "${fs}/surf/lh.pial", "${fs}/surf/lh.sphere" ]
                right: [ sub, "rh", "${fs}/surf/rh.white", "${fs}/surf/rh.pial", "${fs}/surf/rh.shpere" ]
                }

        fs_input = fs_files.left.mix(fs_files.right)

        atlas_files = atlas_dir.multiMap { i -> 
                left: ["lh", "${i}/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii" ]
                right: ["rh", "${i}/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii"]
                }

        atlas_input = atlas_files.left.mix(atlas_files.right)

        i_freesurfer_prep = fs_input.map { sub, hemi, white, pial, sphere -> 
                                        [hemi, sub, white, pial, sphere]}.join(atlas_input,by:[0]).map 
                                        { hemi, sub, white, pial, sphere, atlas -> 
                                        [sub, hemi, white, pial, sphere, atlas]}
        
        freesurfer_prep(i_freesurfer_prep)
                               
        // fs_to_gifti.out.gifti  [sub, hemi, shape.gii]
        // freesurfer_prep.out.current_gifti_sphere [sub, hemi, current_sphere)
        // freesurfer_prep.out.current_midthickness [sub, hemi, current_midthickness]
        // freesurfer_prep.out.new_midthickness [sub, hemi, new_fs_32k_midthickness]

        i_metric_resample = fs_to_gifti.out.gifti.join(freesurfer_prep.out.current_gifti_sphere, by:[0,1]) //[sub,hemi,gifti,current_sphere]
                                                 .map { sub, hemi, gifti, current_sphere -> [hemi, sub, gifti, current_sphere]}.join(atlas_input,by:[0])
                                                 .map { hemi, sub, gifti, current_sphere, atlas -> [ sub, hemi, gifti, current_sphere, atlas]} //[sub,hemi,gifti,current_sphere,atlas]
                                                 .join(freesurfer_prep.out.current_midthickness,by:[0,1]) // [sub, hemi, gifti, current_sphere, atlas, current_midthickness]
                                                 .join(freesurfer_prep.out.new_midthickness, by:[0,1]) // [sub, hemi, gifti, current_sphere, atlas, current_midthickness, new_midthickness]

        metric_resample(i_metric_resample)
	    i_create_dense_scalar = metric_resample.out.e_norm_resampled.groupTuple(by: [0], sort:true)
                                                                    .map { it.flatten() }.map { s, hL, hR, rL, rR -> [s, rL, rR]}
        create_dense_scalar(i_create_dense_scalar)

    emit:
	    sim_dscalar = create_dense_scalar.out.sim_dscalar    


    // The final output to return from this workflow
}
