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
    python /scripts/simulation/placeCoil.py \
        !{coordinate} \
        !{mesh} \
        !{twist} \
        !{sub}_coilPos.npy \
        --qc-file !{sub}_coilQC.html
    '''
}

process simulate{
    label 'simnibs'
    label 'bin'

    input:
    tuple val(sub), path(mesh), path(matsimnibs), path(coil), path(m2m_path), path(fs_path)

    output:
    tuple val(sub), path("${sub}_TMS*.msh"), emit: simMsh
    tuple val(sub), path("${sub}_TMS*.geo"), emit: simGeo
    tuple val(sub), path("fsavg_overlays/lh.*"), emit: leftSim
    tuple val(sub), path("fsavg_overlays/rh.*"), emit: rightSim

    shell:
    '''
    python /scripts/simulation/simulate.py \
        !{mesh} \
        !{matsimnibs} \
        !{coil} \
        --gifti --m2m-path !{m2m_path}
    '''
}


workflow runSimulate{

    take:
        mesh
        coordinates
        fs_path
        m2m_path
        twist
        coil

    main:

        placeCoil(
            mesh.join(coordinates)
                .combine(twist)
        )

        simulate(
            mesh.join(placeCoil.out.matsimnibs)
                .combine(coil)
                .join(m2m_path)
                .join(fs_path)
        )






    emit:
        qcFile = placeCoil.out.qcFile
        simMsh = simulate.out.simMsh
        simGeo = simulate.out.simGeo
        leftGifti = simulate.out
                        .leftSim.map {sub, fields -> [
                            sub,
                            fields.collectEntries{
                                [it.getName().split("\\.")[ -1 ], it]
                            }
                        ]}
        rightGifti = simulate.out
                        .rightSim.map {sub, fields -> [
                            sub,
                            fields.collectEntries{
                                [it.getName().split("\\.")[ -1 ], it]
                            }
                        ]}

}
