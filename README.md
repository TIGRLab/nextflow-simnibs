# SimNIBS to Cifti space Nextflow pipeline

Perform surface mapping from SimNIBS native surface outputs to HCP fs_LR32k space.

## Requirements:

- Install [Nextflow](https://www.nextflow.io/)
- Have singularity containers listed in `parameters.json` available
	- Note that the connectome workbench container must also have freesurfer available
	- You may use the same container more than once if it contains the appropriate software
- **If not in Kimel Lab** - add a `profile` in `config/base.nf.config` for your cluster setup. If you are not using a cluster you can specify `-profile local` to run locally on your computer. The profiles available by default are: ['kimel', 'scc', 'local']

## Pipelines
- `pipeline/mni_transform.nf` - perform SimNIBS Freesurfer to MNI *brain* transformation
- `pipeline/simulate_mni.nf` - perform SimNIBS simulation using MNI coordinates. Automates registration, coordinate mapping, coil placement, and cifti transformation
- `pipeline/simnibs_dscalar.nf` - convert SimNIBS outputs to CIFTI format. `subject_overlays/` must be available for each subject. You can generate these using `map_to_surf` when running a SimNIBS simulation

---

**NOTE**: For all `pipeline` scripts you *must* run as follows:

```
nextflow run pipeline/<pipeline>.nf -c <repo_dir>/config/base.nf.config \
	-params-file <repo_dir>/parameters.json [other args...]
```

- `-c` - provides base configuration, in addition to the supported profiles
- `-params-file` - provides list of singularity containers to use


---


## Example usages


**Run simulation using MNI coordinates and obtain CIFTI files**
```
nextflow run simulate_mni.nf -c config/base.nf.config -params-file parameters.json \
	--out_dir path/to/out --mri2mesh_dir path/to/mri2mesh/outputs --create_cifti \
	--mni_coordinates 30,43,23 --twist 155 --coil /path/to/coil \
	--create_cifti -profile local
```


**Help**

```
nextflow run [mni_transform.nf|simnibs_dscalar.nf|simulate_mni.nf] -c config/base.nf.config --help
```
