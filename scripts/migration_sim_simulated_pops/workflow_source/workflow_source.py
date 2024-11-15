

#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys, re, shutil, itertools
import subprocess as sp
from workflow_templates import *

def migration_simulation_simsfs(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Needs to be run twice (from the config folder: ../configurations/collembola/Entomobrya_nicoleti_simulated )
    
    Will run fastsimcoal on files with 2dSFS, across differing migration divides in a model where there was one migration rate before the devide, and another one after.
    You can set the limits for the migration divide range in the config file.
    It will estimate recent Ne, ancestral Ne, split time, recent gene flow and ancestral gene flow.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    # EDIT CONFIGS
    CONFIG = yaml.safe_load(open(config_file)) #opens yml file for reading
    ACCOUNT: str = CONFIG['account'] # gets the setting after account as string
    SPECIES_NAME: str = CONFIG['species_name']
    #OUTPUT_DIR: str = CONFIG['output_directory_path']
    WORK_DIR: str = CONFIG['working_directory_path']
    TAXONOMY: str = CONFIG['taxonomic_group']
    DESCRIBER_LIST: list = CONFIG['simulation_desriber_list']
    SFS_FILES_LIST: list = CONFIG['SFS_file_list']
    MIGRATION_DIVIDE_INTERVAL: str = CONFIG['generations_interval']
    MAX_MIG_DIVIDE: int = CONFIG['max_generations_migration_divide']
    MIN_MIG_DIVIDE: int = CONFIG['min_generations_migration_divide']

    


    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    ) # what does this do?
    


    if not len(DESCRIBER_LIST) == len(SFS_FILES_LIST):
        print("list of files and describers are not equally long. Please adjust so they match.")

    ############################
    ### MAKE DIRECTORIES
    ############################
    
    new_wd=f'{WORK_DIR}/migration_sim/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/fsc/'
    new_wd_noSingletons=f'{WORK_DIR}/migration_sim_noSingletons/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/fsc/'
    #new_out_mig=f'{OUTPUT_DIR}/migration_sim/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/'
    if not os.path.isdir(new_wd):
        os.makedirs(new_wd)

    if not os.path.isdir(f'{new_wd}/simulated_SFS/'):
        os.makedirs(f'{new_wd}/simulated_SFS/')

    if not os.path.isdir(f'{new_wd}/simulated_SFS/'):
        os.makedirs(f'{new_wd}/simulated_SFS/')

    if not os.path.isdir(new_wd_noSingletons):
        os.makedirs(new_wd_noSingletons)

    if not os.path.isdir(f'{new_wd_noSingletons}/simulated_SFS/'):
        os.makedirs(f'{new_wd_noSingletons}/simulated_SFS/')

    if not os.path.isdir(f'{new_wd_noSingletons}/simulated_SFS/'):
        os.makedirs(f'{new_wd_noSingletons}/simulated_SFS/')
    
    for element in range(0, len(DESCRIBER_LIST)):
        if not os.path.isdir(f'{new_wd}/simulated_SFS/{DESCRIBER_LIST[element]}'):
            os.makedirs(f'{new_wd}/simulated_SFS/{DESCRIBER_LIST[element]}')
        
        if not os.path.isdir(f'{new_wd_noSingletons}/simulated_SFS/{DESCRIBER_LIST[element]}'):
            os.makedirs(f'{new_wd_noSingletons}/simulated_SFS/{DESCRIBER_LIST[element]}')
        





    ##############################################
    ### RUN 2dSFS in fastsímcoal (incl. setup) ###
    ##############################################

    for number in range(0, len(SFS_FILES_LIST)):
        #print( number )
        # get correct files and folder name
        describer = DESCRIBER_LIST[number]
        print(describer)
        file = SFS_FILES_LIST[number]
        #print(file)
        # make name to copy the sfs to
        sfs_file_new = f'{new_wd}/simulated_SFS/{describer}/{os.path.basename(file)}'
        print(sfs_file_new)
        # copy SFS into analysis base-folder
        shutil.copy2(file, sfs_file_new)


    ######################################################
    ### split SFS files, if more than one obs per file ###
    ######################################################
        # template which reads if first sentence starts with 1 or higher number of observations, and makes yml list with paths 
        python_make_yml(SFS_file=sfs_file_new, describer=f'{describer}')

        # if file 
        if os.path.exists(sfs_file_new.replace(".obs", ".yaml")):
        #if os.path.exists(check_SFS_for_multi_obs_target.outputs['new_sfs_list_file']):
            print("The needed yaml file exists. Proceeding.")
            

            # get list of coming output files
            yml_list = yaml.safe_load(open(sfs_file_new.replace(".obs", ".yaml"))) #opens yml file for reading
            replicate_sfs_list: list = yml_list['replicate_sfs_list'] # gets the setting after account as string
            #print(replicate_sfs_list)


            # Split .obs file with multiple obs
            split_file(input_filename=sfs_file_new, 
                        path_out=os.path.dirname(os.path.dirname(replicate_sfs_list[0])), 
                        file_name_base=f'pop1_pop2_sfs2d_rep_X_{describer}_jointMAFpop1_0.obs')
                # add X where the repition index should be



                # run with singletons
            # mapping over all migration_divides in the range specified
            migration_range_list = list(range(MIN_MIG_DIVIDE, MAX_MIG_DIVIDE, MIGRATION_DIVIDE_INTERVAL))
            #print(len(migration_range_list))
            migration_range_dict_list = [{'migration_divide': d } for d in range(MIN_MIG_DIVIDE, MAX_MIG_DIVIDE, MIGRATION_DIVIDE_INTERVAL)] 

            # make dict list for map input:
            combinations_dict_list = [{'SFS_file': a, 'migration_divide': b} for a, b in itertools.product(replicate_sfs_list, migration_range_list)]
            print("Number of FSC runs in line: " + str(len(combinations_dict_list)))


            run_simulated_sfs_map = gwf.map(
                name=create_run_fsc,
                template_func = setup_run_FSC_target,
                inputs = combinations_dict_list,
                extra = {'name_pops': f'pop1_pop2'})

            print("")

 
        else:   
            print("ERROR / OBS: The yml file created by the workflow does not exist yet. Rerun workflow, when first templates have been run (where the yaml is being created).")

        # make new list of SFS files
                
    



    return gwf
