nextflow.preview.dsl = 2

usage = file("${workflow.scriptFile.getParent()}/simulate_mni.usage.txt")
printhelp = params.help

// Import required modules
include {fs_to_mni} from "../modules/mni.nf" params(params)
include {run_simulation} from "../modules/simnibs.nf" params(params)
include {simnibs2cifti} from "../modules/map_to_cifti.nf" params(param)

workflow {
}
