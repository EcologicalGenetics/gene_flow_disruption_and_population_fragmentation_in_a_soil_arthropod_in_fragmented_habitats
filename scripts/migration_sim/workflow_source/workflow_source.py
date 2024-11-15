#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys, re
import subprocess as sp
from workflow_templates import *

def migration_simulation(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Makes 2dSFS for each population pair, and simulates present and past migration rates at specified intervals, using FastSimCoal.
    
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
    OUTPUT_DIR: str = CONFIG['output_directory_path']
    WORK_DIR: str = CONFIG['working_directory_path']
    #VCF_BASE_FOLDER: str = CONFIG['vcf_base_folder']
    TAXONOMY: str = CONFIG['taxonomic_group']
    #BED_PATH: str = CONFIG['bed_files_path']
    #GENOME_PATH: str = CONFIG['reference_genome_path']
    #ANNOTATION_GTF: str = CONFIG['annotation_gtf']
    #VCF_FILES: list = CONFIG['vcf_lists']
    EXCLUDE_POPS: list = CONFIG['pops_exclude_list']
    INCLUDE_POPS: list = CONFIG['pops_INCLUDE_list']
    #COLLECTION_SITES_FILE: str = CONFIG['collection_sites_file']
    MIGRATION_DIVIDE_INTERVAL: str = CONFIG['generations_interval']
    MAX_MIG_DIVIDE: int = CONFIG['max_generations_migration_divide']
    MIN_MIG_DIVIDE: int = CONFIG['min_generations_migration_divide']
    VCF_FILE: str = CONFIG['vcf_file']
    


    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    ) # what does this do?
    

    ### MAKE DIRECTORIES
    #print("making directories")
    # make directories
    new_wd=f'{WORK_DIR}/migration_sim/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/fsc/'
    new_wd_noSingletons=f'{WORK_DIR}/migration_sim_noSingletons/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/fsc/'
    new_out_mig=f'{OUTPUT_DIR}/migration_sim/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/'
    if not os.path.isdir(new_wd):
        os.makedirs(new_wd)
    if not os.path.isdir(f'{new_wd}/2dSFS/adata_prep'):
        os.makedirs(f'{new_wd}/2dSFS/adata_prep')
    if not os.path.isdir(new_out_mig):
        os.makedirs(new_out_mig)
    if not os.path.isdir(new_wd_noSingletons):
        os.makedirs(new_wd_noSingletons)

    # defines new files/directories base on config
    #bed_genes = glob.glob(f'{BED_PATH}/*_genomic.genes.bed')[0]
    #bed_intergenes = glob.glob(f'{BED_PATH}/*_genomic.intergenic.bed')[0]
    #bed_repeats = glob.glob(f'{BED_PATH}/*_genomic.repeats.bed')[0]
    #bed_neutral = f'{new_wd}/{os.path.basename(bed_genes).replace("genes","neutral")}'
    # f'{bed_repeats.replace("repeats","neutral")}'
    

    #####################################
    ### MAKE data ready for 2dSFS BED ###
    #####################################

    
	# Parse individual IDs and their counters
    awk_parse_pops_from_vcf(VCF_FILE, f'{new_wd}/2dSFS/adata_prep/{species_abbreviation(SPECIES_NAME)}_parsed_pops.tsv')
    print(f'Individual IDs and counters saved to {new_wd}/2dSFS/adata_prep/{species_abbreviation(SPECIES_NAME)}_parsed_pops.tsv')
    # Generate all possible pairs where counter2 >= counter1
    awk_create_pairs(f'{new_wd}/2dSFS/adata_prep/{species_abbreviation(SPECIES_NAME)}_parsed_pops.tsv', f'{new_wd}/2dSFS/adata_prep/{species_abbreviation(SPECIES_NAME)}_pop_pairs.tsv')
    print(f'Individual pairs saved to {new_wd}/2dSFS/adata_prep/{species_abbreviation(SPECIES_NAME)}_pop_pairs.tsv')


    # Prepare data for making SFS (slow)
    make_ready_for_2dSFS_pairs = gwf.target_from_template(
            name='prepare_data_sfs',
            template=prepare_data(
                vcf_file=VCF_FILE,
                pop_pairs_file=f'{new_wd}/2dSFS/{species_abbreviation(SPECIES_NAME)}_pop_pairs.tsv',
                output_prefix=species_abbreviation(SPECIES_NAME),
                working_dir=f'{new_wd}/2dSFS/')
        )

    #################
    ### RUN 2dSFS ###
    #################
        # funtion used to create dictionaries, based on files created in this source flow.
    #input_dict_list = create_input_dict_2dSFS(pair_file = f'{new_wd}/2dSFS/adata_prep/{species_abbreviation(SPECIES_NAME)}_pop_pairs.tsv', exclude_list = [])
    input_dict_list = create_input_dict_2dSFS(pair_file = f'{new_wd}/2dSFS/adata_prep/{species_abbreviation(SPECIES_NAME)}_pop_pairs.tsv', include_list = INCLUDE_POPS, exclude_list = EXCLUDE_POPS)
    #input_dict_list_noSingl = create_input_dict_2dSFS_noSingl(pair_file = f'{new_wd}/2dSFS/adata_prep/{species_abbreviation(SPECIES_NAME)}_pop_pairs.tsv', include_list = INCLUDE_POPS, exclude_list = EXCLUDE_POPS)
    #new_wd_noSingletons
    #input_dict_list_noSingletons = create_input_dict_2dSFS(pair_file = f'{new_wd}/2dSFS/adata_prep/{species_abbreviation(SPECIES_NAME)}_pop_pairs.tsv', include_list = INCLUDE_POPS, exclude_list = EXCLUDE_POPS)
    #print(EXCLUDE_POPS)
    print(input_dict_list[0])
    
    neutral_vcf_files_runtemplate_map = gwf.map(
        #name=[d['name'] for d in input_dict_list],
        template_func = pair_2DSFS_map_target,
        inputs = input_dict_list,
        extra = {'out2_file': make_ready_for_2dSFS_pairs.outputs['outfile_dat'], 
                 'working_directory': f'{new_wd}/2dSFS/'})
    
    neutral_vcf_files_runtemplate_map_noSingletons = gwf.map(
        #name=create_run_name_fsc_pair,
        template_func = pair_2DSFS_map_target_noSingletons,
        inputs = input_dict_list,
        extra = {'out2_file': make_ready_for_2dSFS_pairs.outputs['outfile_dat'], 
                 'working_directory': f'{new_wd_noSingletons}/2dSFS/',
                 'working_directory_alldat': f'{new_wd}/2dSFS/'})
    

    #############################################################
    ### implement a FSC run with e.g. 100 as migration divide ###
    #############################################################
        # then afterwards I can consider whether one script should run all pop comparisons for one migration divide (as loop)? or
            # if I should loop around migrations divides for each comparison?
            # or each pair-migration divide in each run - probably easiest?

    #print([{'outfile_path': d['outfile_path']} for d in neutral_vcf_files_runtemplate_map.outputs])
    #SFS_paths=[{'outfile_path': d['outfile_path']} for d in neutral_vcf_files_runtemplate_map.outputs] # as list of dics to input in map.
        #Change key
    #for k in SFS_paths:
     #   k['SFS_file'] = k.pop('outfile_path')
    
    SFS_paths = [d['outfile_path'] for d in neutral_vcf_files_runtemplate_map.outputs] # as list of paths to combine with another list? into dicts
    SFS_paths_noSingletons = [d['sfs_newrun_outfile_path'] for d in neutral_vcf_files_runtemplate_map_noSingletons.outputs] # as list of paths to combine with another list? into dicts
    names_list = [d['name'] for d in input_dict_list] # as list of paths to combine with another list? into dicts
    input_dict_list_FSC = [{'SFS_file': f, 'name_pops': n} for f, n in zip(SFS_paths, names_list)]
    input_dict_list_FSC_noSingletons = [{'SFS_file': f, 'name_pops': n} for f, n in zip(SFS_paths_noSingletons, names_list)]

    print(input_dict_list_FSC[0])
    print(input_dict_list_FSC_noSingletons[0])
    

    # set up folder structure, copy .obs into it and change name
    for migration_divide in range(MIN_MIG_DIVIDE, MAX_MIG_DIVIDE, MIGRATION_DIVIDE_INTERVAL):
        #print(migration_divide)
        #migration_divide=600
        #input_dict_list_FSC = [{'SFS_file': f, 'name_pops': n} for f, n in zip(SFS_paths, names_list)]
        setup_run_FastSimCoal = gwf.map(
            name = create_run_name_fsc,
            template_func = setup_run_FSC_map_target ,
            #inputs = input_dict_list_FSC[0:150],
            inputs = input_dict_list_FSC,
            extra = {'migration_divide': migration_divide})

        # run with no singletons:
        setup_run_FastSimCoal = gwf.map(
            name = create_run_name_fsc_noSingl,
            template_func = setup_run_FSC_map_target_nosingletons ,
            #inputs = input_dict_list_FSC_noSingletons[0:150],
            inputs = input_dict_list_FSC_noSingletons,
            extra = {'migration_divide': migration_divide})
    


    return gwf
