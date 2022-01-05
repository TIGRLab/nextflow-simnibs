#!/bin/bash
#SBATCH --partition=low-moby
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=2G
#SBATCH --time=30:00:00
#SBATCH --export=ALL
#SBATCH --job-name="fs_2_32k_LR"
#SBATCH --output=/projects/ttan/UBC-TMS/code/fs_resampling_2_32k_LR_modified_%j.txt
#SBATCH --array=1-2

module load connectome-workbench/1.4.1
module load freesurfer/6.0.1

# Specific inputs
study="UBC-TMS"
coil_type="MRI-B91"
sublist=/projects/ttan/${study}/sublist.txt
fs_indir=/projects/ttan/${study}/simnibs/mri2mesh
resamp_indir=/projects/ttan/${study}/standard_mesh_atlases/resample_fsaverage

(
index() {
   head -n $SLURM_ARRAY_TASK_ID $sublist \
   | tail -n 1
}
sub_outdir=/projects/ttan/${study}/simnibs/fsaverage_LR32k/`index`
mkdir ${sub_outdir}

# Convert Norm E scalar in freesurface surface to gifti format
# Path to simulation outputs
simul_indir=/projects/ttan/${study}/simnibs/simulation_outputs

fs2gifti() {
    mris_convert -c ${simul_indir}/`index`/subject_overlays/${idx,,}h.`index`_TMS_1-0001_${coil_type}_*scalar.central.E.norm \
                    ${fs_indir}/fs_`index`/surf/${idx,,}h.white \
                    ${simul_indir}/`index`/subject_overlays/${idx,,}h.`index`_TMS_1-0001_${coil_type}_scalar_E_norm.func.gii

    echo "Running: " \
        mris_convert -c ${simul_indir}/`index`/subject_overlays/${idx,,}h.`index`_TMS_1-0001_${coil_type}_*scalar.central.E.norm \
        ${fs_indir}/fs_`index`/surf/${idx,,}h.white \
        ${simul_indir}/`index`/subject_overlays/${idx,,}h.`index`_TMS_1-0001_${coil_type}_scalar_E_norm.func.gii
}

# Freesurfer native individual data to fs_LR

# Step 1 prepare left and right surface for resampling process
fsresample() {
    wb_shortcuts -freesurfer-resample-prep  ${fs_indir}/fs_`index`/surf/${idx,,}h.white \
                                            ${fs_indir}/fs_`index`/surf/${idx,,}h.pial \
                                            ${fs_indir}/fs_`index`/surf/${idx,,}h.sphere.reg \
                                            ${resamp_indir}/fs_LR-deformed_to-fsaverage.${idx}.sphere.32k_fs_LR.surf.gii \
                                            ${sub_outdir}/`index`.${idx,,}h.midthickness.surf.gii \
                                            ${sub_outdir}/`index`.${idx}.mdthickness.32k_fs_LR.surf.gii ${sub_outdir}/`index`.${idx,,}h.sphere.reg.surf.gii

    echo "Running: " \
          wb_shortcuts -freesurfer-resample-prep  ${fs_indir}/fs_`index`/surf/${idx,,}h.white \
                                            ${fs_indir}/fs_`index`/surf/${idx,,}h.pial \
                                            ${fs_indir}/fs_`index`/surf/${idx,,}h.sphere.reg \
                                            ${resamp_indir}/fs_LR-deformed_to-fsaverage.${idx}.sphere.32k_fs_LR.surf.gii \
                                            ${sub_outdir}/`index`.${idx,,}h.midthickness.surf.gii \
                                            ${sub_outdir}/`index`.${idx}.mdthickness.32k_fs_LR.surf.gii ${sub_outdir}/`index`.${idx,,}h.sphere.reg.surf.gii
}

# Step 2 resample gifti surface to 32k
giftiresamp() {
    gifti_in=${simul_indir}/`index`/subject_overlays/${idx,,}h.`index`_TMS_1-0001_${coil_type}_scalar_E_norm.func.gii
    wb_command -metric-resample ${gifti_in} \
                                ${sub_outdir}/`index`.${idx,,}h.sphere.reg.surf.gii \
                                ${resamp_indir}/fs_LR-deformed_to-fsaverage.${idx}.sphere.32k_fs_LR.surf.gii \
                                ADAP_BARY_AREA ${sub_outdir}/`index`.E.norm.${idx}.32k_fs_LR.func.gii \
                                -area-surfs ${sub_outdir}/`index`.${idx,,}h.midthickness.surf.gii \
                                ${sub_outdir}/`index`.${idx}.mdthickness.32k_fs_LR.surf.gii

    echo "Running: " \
        wb_command -metric-resample ${gifti_in} \
                                ${sub_outdir}/`index`.${idx,,}h.sphere.reg.surf.gii \
                                ${resamp_indir}/fs_LR-deformed_to-fsaverage.${idx}.sphere.32k_fs_LR.surf.gii \
                                ADAP_BARY_AREA ${sub_outdir}/`index`.E.norm.${idx}.32k_fs_LR.func.gii \
                                -area-surfs ${sub_outdir}/`index`.${idx,,}h.midthickness.surf.gii \
                                ${sub_outdir}/`index`.${idx}.mdthickness.32k_fs_LR.surf.gii
}

# Run step by step
for idx in L R;
do
index &&
    fs2gifti &&
    fsresample &&
    giftiresamp
done

# Create dense scalar from L and R hemisphere fs_LR32k

wb_command -cifti-create-dense-scalar ${sub_outdir}/`index`.norm.E.32k_fs_LR.dscalar.nii \
                                        -left-metric ${sub_outdir}/`index`.E.norm.L.32k_fs_LR.func.gii \
                                        -right-metric ${sub_outdir}/`index`.E.norm.R.32k_fs_LR.func.gii

echo "Running: " wb_command -cifti-create-dense-scalar ${sub_outdir}/`index`.norm.E.32k_fs_LR.dscalar.nii \
                                        -left-metric ${sub_outdir}/`index`.E.norm.L.32k_fs_LR.func.gii \
                                        -right-metric ${sub_outdir}/`index`.E.norm.R.32k_fs_LR.func.gii
)
