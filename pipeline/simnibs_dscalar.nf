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

parser.addArgument("--atlas_dir",
				   "Path to HCP atlas dir",
				   params.atlas_dir.toString(),
				   "ATLAS_DIR")

parser.addArgument("--out_dir",
                   "Path to output directory",
                   params.out_dir.toString(),
                   "OUT_DIR")

parser.addArgument("--freesurfer_img",
				   "Path to Freesurfer singularity image",
				   params.freesurfer_img.toString(),
				   "FS_IMG")

parser.addArgument("--connectome_img",
				   "Path to Connectome Workbench singularity image",
				   params.connectome_img.toString(),
				   "CONNECTOME_IMG")

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
log.info("Using FS image file: $params.freesurfer_img")
log.info("Using connectome-workbench image file: $params.connectome_img")
log.info("Input HCP atlas directory: $params.atlas_dir")
if (params.subjects) {
    log.info("Subject list file provided: $params.subjects")
}

fs_input = channel.fromPath("$params.mri2mesh_dir/fs_sub-*", type: 'dir')
						.map{i -> [i.getBaseName() - ~/^fs_/, i]}

simul_input = channel.fromPath("$params.simulation_dir/sub-*", type: 'dir')
						.map{i -> [i.getBaseName(), 
									new FileNameByRegexFinder().getFileNames("${i}",".subject_overlays/.*E.norm\$")]} //use branch operator

atlas_input = Channel.fromPath("$params.atlas_dir", type:"dir")

if (params.subjects){
    subjects_channel = Channel.fromPath(params.subjects)
                            .splitText(){it.strip()}

    fs_input = fs_input.join(subjects_channel,by:[0])
	simul_input = simul_input.join(subjects_channel,by:[0])

}

workflow {
	main:
		simnibs2cifti(simul_input, fs_input, atlas_input)
}