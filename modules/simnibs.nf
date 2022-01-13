/*
* SimNIBS simulation operations
*/

process generate_matsimnibs{
/*
* Generate a matsimnibs matrix from a coordinates and a twist angle
* Twist angle is defined in the BrainSight convention
*
* Arguments;
*   sub (str): Subject ID
*   msh (path): Path to .msh file
*   x (int): x coordinate
*   y (int): y coordinate
*   z (int): z coordinate
*   theta (int): Brainsight twist angle
*
* Output:
*   matsimnibs (path): Matsimnibs file
*/

    label 'simnibs'

    input:
    tuple val(sub), path(msh),\
    val(x), val(y), val(z), val(theta)

    output:
    tuple val(sub), path(matsimnibs), emit: matsimnibs

    shell:
    '''
    '''
}

process simulate {
/*
* Run SimNIBS simulation
*/
}

workflow run_simulation {

/*
* Run SimNIBS simulation using input coordinates
*
* Arguments:
*   msh (channel): [subject, msh] Subject .msh files
*   m2m (channel): [subject, m2m] Subject SimNIBS m2m directories
*   orientation (channel): [subject, x,y,z,theta] Coil orientations to run simulation on
*
* Outputs:
*   sim_msh (channel): [subject, sim_msh] Subject simulation .msh files
*   sim_fs (channel): [subject, hemi, sim_fs] Subject fsaverage simulation files
*/

    take:
        msh
        m2m
        orientation


    main:
        

    emit:
        sim_msh
        sim_fs
        
}
