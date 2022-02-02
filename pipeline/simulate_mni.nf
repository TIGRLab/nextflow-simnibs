nextflow.preview.dsl = 2

include {registerFreesurferToMNI; antsApplyWarpToCoordinates; antsToNumpy } from "../modules/mni.nf"
include { simnibs2cifti } from "../modules/map_to_cifti.nf"
include { getArgumentParser } from "../lib/args"
include { runSimulate } from "../modules/simulate.nf"

parser = getArgumentParser(
    title: "SimNIBS MNI Simulations",
    description: "Perform MNI coordinate transformation and simulation on SimNIBS output",
    scriptName: "${workflow.scriptName}".toString(),
)

parser.addArgument("--mri2mesh_dir",
    "Path to mri2mesh directory",
    params.mri2mesh_dir.toString(),
    "MRI2MESH")

parser.addArgument("--ants_img",
    "Path to ANTS Singularity image",
    params.ants_img.toString(),
    "ANTS_IMG")

parser.addArgument("--mni_coordinates",
    "MNI coordinates to simulate at",
    params.mni_coordinates.toString(),
    "X,Y,Z")

parser.addArgument("--twist",
    "Coil Twist Angle (BrainSight convention)",
    params.twist.toString(),
    "TWIST_ANGLE")

parser.addArgument("--coil",
    "Coil definition file (.nii.gz)",
    params.coil.toString(),
    "COIL_NII_GZ")

parser.addOptional("--subjects", "Path to text file with list of subjects to run")
parser.addOptional("--warps_file",
    "Path to text file containing (subject, warp, inverseWarp) tuples",
    "WARPS_FILES")
parser.addOptional("--mni_standard",
    "Path to MNI standard registration target, required if --warp_files not provided",
    "MNI_STANDARD")

parser.addOptional("--create_cifti", "Generate CIFTI outputs")


missingArgs = parser.isMissingRequired()
if (params.help){
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

if (params.warps_file){
    log.info("Using warps map: $params.warps_file")
} else if (!params.mni_standard) {
    log.error("MNI standard registration target not provided!")
    log.error("--mni_standard required if --warp_files not used")
    System.exit(1)
}

if (params.subjects){
    log.info("Using subjects file: $params.subjects")
}

if (params.create_cifti){
    log.info("Will output CIFTI files")
}

fs_input = Channel.fromPath("$params.mri2mesh_dir/fs_sub-*", type: 'dir')
            .map { i -> [i.getBaseName() - ~/^fs_/, i] }
m2m_input = Channel.fromPath("$params.mri2mesh_dir/m2m_sub-*", type: 'dir')
            .map { i -> [i.getBaseName() - ~/^m2m_/, i] }
mesh_input = Channel.fromPath("$params.mri2mesh_dir/*.msh")
            .map { i -> [i.getBaseName() - ~/.msh$/, i] }

// Filter subjects
if (params.subjects){
    subjects = Channel.fromPath(params.subjects)
            .splitText() { it.strip() }
    fs_input = fs_input.join(subjects)
    m2m_input = m2m_input.join(subjects)
    mesh_input = mesh_input.join(subjects)
} else {
    subjects = fs_input.map { s,f -> s }
}

// Check to make sure dataset is complete
fs_input.join(m2m_input, failOnMismatch: true)
        .join(mesh_input, failOnMismatch: true)

workflow getOrCreateWarps{
/*
* Obtain Freesurfer to MNI warps
* Outputs:
*   warps (channel): [subject, warp, inverseWarp] hash map
*/

    main:

        if (params.warps_file){
            warps = Channel.fromPath(params.warps_file)
                    .splitCsv(header: ['subject', 'warp', 'inverseWarp'])

            // Make sure we have a warp for each subject
            subjects.join(warps, failOnMismatch: true)
        } else {
            log.info("No warps provided, running registration")
            registerFreesurferToMNI(fs_input, Channel.of(params.mni_standard))
            warps = registerFreesurferToMNI.out.warp
                        .join(registerFreesurferToMNI.out.inverseWarp)
                        .map { s, w, iw ->
                            [
                                subject: s,
                                warp: w,
                                inverseWarp: iw
                            ]
                        }
        }

    emit:
        warps = warps
}

workflow createCifti {
    take:
        simulation_files
        fs_dir

    main:
        simnibs2cifti(
            simulation_files,
            fs_dir,
            Channel.of(params.atlas_dir)
        )

        //add a publish command

    emit:
        dscalar = simnibs2cifti.out.sim_dscalar
}

workflow {

    main:
        getOrCreateWarps()

        coordinates = Channel.of(params.mni_coordinates).map { it.split(",") }
        antsApplyWarpToCoordinates(
            coordinates,
            getOrCreateWarps.out.warps.map { [it.subject, it.warp] }
        )

        antsToNumpy(antsApplyWarpToCoordinates.out.warpedCoordinates)

        runSimulate(
            mesh_input,
            antsToNumpy.out.coords,
            fs_input,
            m2m_input,
            Channel.of(params.twist),
            Channel.fromPath(params.coil)
           )

        if (params.create_cifti){
            createCifti(
                runSimulate.out.rightSurf
                    .mix(runSimulate.out.leftSurf)
                    .map{ s, surfs -> [s, surfs.norm] },
                fs_input
            )
        }
}
