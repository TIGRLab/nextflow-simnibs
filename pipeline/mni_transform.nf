nextflow.enable.dsl = 2

include {registerFreesurferToMNI} from "../modules/mni.nf" params(params)
include { validateArgs; getArgumentParser } from "../lib/args"

parser = getArgumentParser(
    title: "Freesurfer MNI Transform",
    description: "Perform antsRegistration MNI transformation on Freesurfer T1",
    scriptName: "${workflow.scriptName}".toString(),
    note: """\
    Any parameter can be specified in a Nextflow config file with the following syntax:

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

parser.addArgument("--out_dir",
                   "Path to output directory",
                   params.out_dir.toString(),
                   "OUT_DIR")

parser.addArgument("--ants_img",
                   "Path to ANTS singularity image",
                   params.ants_img.toString(),
                   "ANTS_IMG")

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
    System.exit(1)
}


log.info("Input mri2mesh directory: $params.mri2mesh_dir")
log.info("Using ANTS image file: $params.ants_img")

if (params.subjects) {
    log.info("Subject list file provided: $params.subjects")
}

// Extract subject directories to run
input_channel = Channel.fromPath("$params.mri2mesh_dir/fs_sub-*", type: 'dir')
                        .map{i -> [i.getBaseName(), i]}

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

    input:
    tuple val(sub), path(warped), path(warp), path(inverseWarp)

    output:
    tuple val(sub), path(warped), path(warp), path(inverseWarp)

    shell:
    '''
    echo "Moving files into !{params.out_dir}
    '''
}

workflow {
    main:
        registerFreesurferToMNI(
            input_channel,
            Channel.of(params.mni_standard)
        )
}
