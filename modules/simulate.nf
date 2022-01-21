/*
* Simulation modules
*/


process placeCoil{
/*
* Place coil on scalp surface given a coordinate
*/
    label 'simnibs'
    label 'bin'

    input:
    tuple val(sub), path(mesh), path(coordinate), val(twist)

    output:
    tuple val(sub), path("${sub}_coilPos.npy"), emit: matsimnibs
    tuple val(sub), path("${sub}_coilQC.html"), emit: qcFile, optional: true

    shell:
    '''
    python /scripts/placeCoil.py \
        !{coordinate} \
        !{mesh} \
        !{twist} \
        !{sub}_coilPos.npy

    # TODO: After adding container with meshplot
    #--qc-file !{sub}_coilQC.html
    '''
}


workflow runSimulate{

    take:
        mesh
        coordinates
        twist

    main:

        placeCoil(
            mesh.join(coordinates)
                .combine(twist)
        )
}
