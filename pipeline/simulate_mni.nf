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

parser.addArgument("--out_dir",
    "Path to output directory",
    params.out_dir.toString(),
    "OUTDIR")

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
parser.addOptional("--didt_file",
    "Path to CSV file containing (subject, DIDT) entries",
    "DIDT_FILE")

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

process publishRegistration{
/*
* Publish registration outputs
*
* Argument:
*   warped: T1 warped into MNI
*   warp: Forward warp file
*   inverseWarp: Inverse warp file
*   qcSvg: Quality control SVG file
*/
    publishDir path: "${params.out_dir}/mniSimulation/${sub}/mniRegistration", \
                mode: 'copy', \
                overwrite: true

    input:
    tuple val(sub), path(warped), path(warp), path(inverseWarp), path(qcSvg)

    output:
    tuple val(sub), path(warped), path(warp), path(inverseWarp), path(qcSvg)

    shell:
    '''
    echo "Moving outputs into !{params.out_dir}/mniSimulation/!{sub}/mniRegistration"
    '''
}

process publishSimulations{
/*
* Publish outputs into `params.out_dir`
*
* Arguments:
*   simFile: SimNIBS .msh file
*   simGeo: SimNIBS coil placement .geo file
*   leftSurf: SimNIBS native subject overlay lh
*   rightSurf: SimNIBS native subject overlay rh
*   leftFsavgSurf: SimNIBS fsavg overlay lh
*   rightFsavgSurf: SimNIBS fsavg overlay rh
*   qcHtml: Coil placement QC HTML file
*/

    publishDir path: "${params.out_dir}/mniSimulation/${sub}/simulations", \
                mode: 'copy', \
                overwrite: true

    input:
    tuple val(sub), path(simFile), path(simGeo),\
    path(leftSurf), path(rightSurf),\
    path(leftFsavgSurf), path(rightFsavgSurf),\
    path(qcHtml)

    output:
    tuple val(sub), path(simFile), path(simGeo),\
    path(leftSurf), path(rightSurf),\
    path(leftFsavgSurf), path(rightFsavgSurf),\
    path(qcHtml)

    shell:
    '''
    echo "Publishing to !{params.out_dir}/mniSimulation/!{sub}/simulations"
    '''

}

process publishCifti{
/* Publish CIFTI dscalar outputs
* Arguments:
*   dscalar: CIFTI dscalar file
*/
    publishDir path: "${params.out_dir}/mniSimulation/${sub}/", \
                mode: 'copy', \
                overwrite: true

    input:
    tuple val(sub), path(cifti)

    output:
    tuple val(sub), path(cifti)

    shell:
    '''
    echo "Publishing to !{params.out_dir}/mniSimulation/!{sub}"
    '''
}

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
                        .join(registerFreesurferToMNI.out.warped)
                        .join(registerFreesurferToMNI.out.qcImage)
                        .map { s, w, iw, warped, qc ->
                            [
                                subject: s,
                                warp: w,
                                inverseWarp: iw,
                                warped: warped, 
                                qcImage: qc
                            ]
                        }
            publishRegistration(
                warps.map { w -> [w.subject, w.warped, w.warp, w.inverseWarp, w.qcImage] }
            )
                            
        }

    emit:
        warps = warps
}

workflow getOrCreateDosage{
    main:
        if (params.didt_file){
            dosages = Channel.fromPath(params.didt_file)
                        .splitCsv(header: ['subject', 'didt'])
            subjects.join(dosages, failOnMismatch: true)
        } else {
            dosages = subjects.combine([params.default_dose])
        }

    emit:
        dosages = dosages
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

        publishCifti(simnibs2cifti.out.sim_dscalar)

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

        getOrCreateDosage()

        runSimulate(
            mesh_input,
            antsToNumpy.out.coords,
            fs_input,
            m2m_input,
            Channel.of(params.twist),
            Channel.fromPath(params.coil),
            getOrCreateDosage.out.dosages
           )

        publishSimulations(
            runSimulate.out.simMsh
                .join(runSimulate.out.simGeo)
                .join(runSimulate.out.leftSurf.map{ sub, surf -> [sub, surf.norm]})
                .join(runSimulate.out.rightSurf.map{ sub, surf -> [sub, surf.norm]})
                .join(runSimulate.out.leftFsavgSurf.map{ sub, surf -> [sub, surf.norm]})
                .join(runSimulate.out.rightFsavgSurf.map{ sub, surf -> [sub, surf.norm]})
                .join(runSimulate.out.qcFile)
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
