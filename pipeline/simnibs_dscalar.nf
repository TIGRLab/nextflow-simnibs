nextflow.enable.dsl = 2

include {simnibs2cifti} from "../modules/map_to_cifti.nf" params(params)
include { validateArgs; getArgumentParser} from "../lib/args"

parser = getArgumentParser(
	title: "Simnibs To Cifti Mapping",
	description: "Mapping simnibs space to ciftify surface on Simulated E-Field",
	scriptName: "${workflow.scriptName}".toString(),
	note: """\
	Any parameter can be specificed in Nextflow config file with the following syntax:

	```
	params.<parameter> = <value>

	// i.e
	params.mri2mesh_dir = '/path/to/mri2mesh'
	"""
)

parser.addArgument("--mri2mesh_dir",
                   "Path to mri2mesh output dir",
                   params.mri2mesh_dir.toString(),
                   "MRI2MESH")

parser.addArgument("--simulation_dir",
				   "Path to simnibs simulation output dir",
				   params.simulation_dir.toString(),
				   "SIMULATION")

parser.addArgument("--HCP_atlas_dir",
				   "Path to HCP atlas dir",
				   params.HCP_atlas_dir.toString(),
				   "HCPATLAS"

parser.addArgument("--out_dir",
                   "Path to output directory",
                   params.out_dir.toString(),
                   "OUT_DIR"

parser.addArgument("--fs_img",
				   "Path to Freesurfer singularity image",
				   params.fs_img.toString(),
				   "FS_IMG")

parser.addArgument("--wb_img",
				   "Path to Connectome Workbench singularity image",
				   params.fs_img.toString(),
				   "WB_IMG")

parser.addOptional("--subjects", "Path to text file with list of subjects to run")

missingArgs = parser.isMissingRequired()

if (params.help) {
	print(parser.makeDoc())
	System.exit(0)
}

if (missingArgs) {
    log.error("Missing parameters!")
    missingArgs.each{ log.error("Missing ${it}") }
    print(parser.makeDoc())
    System.exit(0)
}

log.info("Input mri2mesh directory: $params.mri2mesh_dir")
log.info("Input simulation directory: $params.simulation_dir")
log.info("Using ANTS image file: $params.ants_img")
log.info("Using FS image file: $params.fs_img")
log.info("Using wb image file: $params.wb_img")

if (params.subjects) {
    log.info("Subject list file provided: $params.subjects")
}

// Extract subject directories to run
input_dirs = new File(params.mri2mesh_dir).list()
simul_dirs = new File(params.simulation_dir).list()

fs_channel = channel.fromPath("$params.mri2mesh_dir/fs_sub-*", type: 'dir')
						.map{i -> [i.getBaseName(), i]}

simul_channel = channel.fromPath("$params.simul_dir/sub-*", type: 'dir')
						.map{i -> ["fs_" + i.getBaseName(), ["${i}/subject_overlays/*E.norm"]]} //use branch operator

atlas_channel = Channel.value("$params.HCP_atlas_dir", type:"dir")

input_channel = fs_channel.join(simul_channel).combine(atlas_channel)

if (params.subjects){
    subjects_channel = Channel.fromPath(params.subjects)
                            .splitText(){it.strip()}
                            .map{ s -> "fs_${s}"}

    input_channel = input_channel.join(subjects_channel)

}

process publish{
    publishDir path: "$params.out_dir", \
               mode: 'move', \
               overwrite: true
               }