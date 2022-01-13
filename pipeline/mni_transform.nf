nextflow.enable.dsl = 2

include {fs_to_mni} from "../modules/mni.nf" params(params)

usage = file("${workflow.scriptFile.getParent()}/mni_transform.usage.txt")
bindings = [
    "subjects": "${params.subjects}",
    "mni_standard": "${params.mni_standard}",
    "mri2mesh_dir": "${params.mri2mesh_dir}",
    "ants_simg": "${params.ants_simg}"
]

engine = new groovy.txt.SimpleTemplateEngine()
toprint = engine.createTemplate(usage.text).make(bindings)
printhelp = params.help

req_param = [
    "--mri2mesh_dir": "${params.bids}",
    "--mni_standard": "${params.mni_standard}",
    "--out_dir": "${params.out_dir}"
]

req_config_param = [
    "--ants_simg": "${params.ants_simg}",
]

missing_req = req_param.grep{ (it.value == null || it.value == "") }
missing_req_config = req_config_param.grep{ (it.value == null || it.value == "") }

if (missing_req) {
    log.error("Missing required command-line arguments!")
    missing_req.each{ log.error("Missing ${it.key}") }
    printhelp = true
}

if (missing_req) {
    log.error("Config file missing required parameters!")
    missing_req_config.each{ log.error("Missing ${it.key}") }
    printhelp = true
}

if (printhelp) {
    print(toprin)
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
output_dirs = new File(params.out_dir).list()

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
