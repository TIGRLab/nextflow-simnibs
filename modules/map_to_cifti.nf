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


process freesurfer_prep{
    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(fs_white), path(fs_pial), path(fs_sphere), path(resample_sphere)

    output:
    tuple val(sub), val(hemi), path("${sub}.${hemi}.midthickness.surf.gii"), emit: current_midthickness
    tuple val(sub), val(hemi), path("${sub}.${hemi}.midthickness.32k_fs_LR.surf.gii"), emit: new_midthickness
    tuple val(sub), val(hemi), path("${hemi}.sphere.reg.surf.gii"), emit: current_gifti_sphere
    
    shell:
    '''
    wb_shortcuts -freesurfer-resample-prep !{fs_white} !{fs_pial} !{fs_sphere} !{resample_sphere} \
    !{sub}.!{hemi}.midthickness.surf.gii \
    !{sub}.!{hemi}.midthickness.32k_fs_LR.surf.gii \
    !{hemi}.sphere.reg.surf.gii
    ''' 
}

process metric_resample{
    label 'connectome'

    input:
    tuple val(sub), val(hemi), path(gifti_file), path(current_sphere_file), path(resample_sphere), path(current_midthickness_file), path(new_midthickness_file)

    output:
    tuple val(sub), val(hemi), path("${sub}.E.norm.${hemi}.32k_fs_LR.func.gii"), emit: new_enorm_resampled
    
    shell:
    '''
    wb_command -metric-resample !{gifti_file} !{current_sphere_file} !{resample_sphere} \
    ADAP_BARY_AREA !{sub}.E.norm.!{hemi}.32k_fs_LR.func.gii \
    -area-surfs !{current_midthickness_file} !{new_midthickness_file}
    ''' 
}

process create_dense_scalar{
    label 'connectome'

    input:
    tuple val(sub), path(resampled_L_file), path(resampled_R_file)

    output:
    tuple val(sub), path("${sub}.norm.E.32k_fs_LR.dscalar.nii"), emit:sim_dscalar

    shell:
    '''
    wb_command -cifti-create-dense-scalar !{sub}.norm.E.32k_fs_LR.dscalar.nii \
    -left-metric !{resampled_L_file} -right-metric !{resampled_R_file}
    '''
}


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
        atlas_dir

    main:

	    split_hemi = fs_dir.multiMap { sub, fs -> 
                                left: [sub, "lh", "${fs}/surf/lh.white" ]
						        right: [sub, "rh", "${fs}/surf/rh.white" ]
				}
	
        sim_inputs = simfiles.transpose().map { s, sim -> [s, (sim =~ ~/[l,r]h/)[0], sim]}

        i_fs_to_gifti = split_hemi.left.mix(split_hemi.right).join(sim_inputs,by:[0,1])
        fs_to_gifti(i_fs_to_gifti)

        fs_files = fs_dir.multiMap { sub, fs -> 
                left: [ sub, "lh", "${fs}/surf/lh.white", "${fs}/surf/lh.pial", "${fs}/surf/lh.sphere.reg" ]
                right: [ sub, "rh", "${fs}/surf/rh.white", "${fs}/surf/rh.pial", "${fs}/surf/rh.sphere.reg" ]
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

        i_metric_resample = fs_to_gifti.out.gifti.join(freesurfer_prep.out.current_gifti_sphere, by:[0,1])
                                                 .map { sub, hemi, gifti, current_sphere -> [hemi, sub, gifti, current_sphere]}.join(atlas_input,by:[0])
                                                 .map { hemi, sub, gifti, current_sphere, atlas -> [ sub, hemi, gifti, current_sphere, atlas]}
                                                 .join(freesurfer_prep.out.current_midthickness,by:[0,1])
                                                 .join(freesurfer_prep.out.new_midthickness, by:[0,1])

        metric_resample(i_metric_resample)
    
	i_create_dense_scalar = metric_resample.out.new_enorm_resampled.map{s,h,f -> [s,f]}.groupTuple(by: 0, sort: {it.baseName}).map { it.flatten() }

        create_dense_scalar(i_create_dense_scalar)

    emit:
        sim_dscalar = create_dense_scalar.out.sim_dscalar    


}
