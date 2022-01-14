nextflow.enable.dsl = 2

include {fs_to_mni} from "../modules/mni.nf" params(params)
include { validateArgs; getUsage } from "../lib/args"

bindings = [
    "subjects": "${params.subjects}",
    "mni_standard": "${params.mni_standard}",
    "mri2mesh_dir": "${params.mri2mesh_dir}",
    "ants_simg": "${params.ants_simg}"
]

usage = getUsage(
    "${workflow.scriptFile.getParent()}/mni_transform.usage.txt",
    bindings)
 
req_param = [
    "--mri2mesh_dir": params.bids,
    "--mni_standard": params.mni_standard,
    "--out_dir": params.out_dir,
    "--ants_simg": params.ants_simg
]
missing_args = validateArgs(req_param)

if (missing_args) {
    log.error("Missing parameters!")
    missing_args.each{ log.error("Missing ${it.key}") }
    print(usage)
    System.exit(0)
}

if (params.help) {
    print(usage)
    System.exit(0)
}

log.info("Input mri2mesh directory: $params.mri2mesh_dir")
log.info("Using ANTS image file: $params.ants_simg")

if (params.subjects) {
    log.info("Subject list file provided: $params.subjects")
}

// Extract subject directories to run
all_dirs = file(params.bids).list()
input_dirs = new File(params.mri2mesh_dir).list()

input_channel = Channel.fromPath("$params.mri2mesh_dir/fs_sub-*", type: 'dir')
                        .map{i -> i.getBaseName()}

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
        fs_to_mni(
            input_channel.map { f -> [ f - ~/^fs_/, f ] },
            params.mni_standard
        )
}
