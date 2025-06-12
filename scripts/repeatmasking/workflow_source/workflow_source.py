#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def repeat_masking_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Uses RepeatModeler, RepeatMasker and data from available databases to identify and soft-mask repetitive region in a genome assembly.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    GENOME_ASSEMBLY: str = config['genome_assembly']
    REPEAT_DATABASE: str = config['repeat_database']
    WORKING_DIR: str = config['working_directory_path']
    OUTPUT_DIR: str = config['output_directory_path']
    ALT_NAME: str = config['alternative_folder_name']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    top_dir = f'{WORKING_DIR}/{ALT_NAME.replace(" ", "_")}/repeatmasking' if ALT_NAME else f'{WORKING_DIR}/{SPECIES_NAME.replace(" ", "_")}/repeatmasking'
    repeatmodeler_dir = f'{top_dir}/RepeatModeler'
    repeatmodeler_db_dir = f'{repeatmodeler_dir}/RM_DB_{SPECIES_NAME.replace(" ", "_")}'
    repeatmasker_dir = f'{top_dir}/RepeatMasker'
    repeatmasker_repeatmodeler_run_dir = f'{repeatmasker_dir}/RepMod_run'
    repeatmasker_repbase_run_dir = f'{repeatmasker_dir}/RepBase_arthropoda_run'

    os.makedirs(name=top_dir, exist_ok=True)
    os.makedirs(name=repeatmodeler_dir, exist_ok=True)
    os.makedirs(name=repeatmodeler_db_dir, exist_ok=True)
    os.makedirs(name=repeatmasker_dir, exist_ok=True)
    os.makedirs(name=repeatmasker_repeatmodeler_run_dir, exist_ok=True)
    os.makedirs(name=repeatmasker_repbase_run_dir, exist_ok=True)

    build_database = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_build_database',
        template=build_repeatmodeller_database(
            genome_assembly_file=GENOME_ASSEMBLY,
            working_directory=repeatmodeler_db_dir,
            species_name=SPECIES_NAME
        )
    )
    
    repmod = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_RepeatModeler',
        template=repeatmodeler(
            database=build_database.outputs['db_files'],
            working_directory=repeatmodeler_dir,
            species_name=SPECIES_NAME
        )
    )

    repeatmasker_database = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_RepeatMasker_repbase_arthropoda',
        template=repeatmasker(
            genome_assembly_file=GENOME_ASSEMBLY,
            library_file=REPEAT_DATABASE,
            working_directory=repeatmasker_repbase_run_dir,
        )
    )

    repeatmasker_repmod = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_RepeatMasker_repmod',
        template=repeatmasker(
            genome_assembly_file=repeatmasker_database.outputs['repmaskout'][1],
            library_file=repmod.outputs['all'],
            working_directory=repeatmasker_repeatmodeler_run_dir,
        )
    )

    combine = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_combine_RepeatMasker_output',
        template=combine_repeatmasker_runs(
            repeatmasker_run1=repeatmasker_database.outputs['repmaskout'],
            repeatmasker_run2=repeatmasker_repmod.outputs['repmaskout'],
            library_file=REPEAT_DATABASE,
            working_directory=repeatmasker_dir,
            species_name=SPECIES_NAME
        )
    )

    softmasking = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_softmasking',
        template=mask_assembly(
            genome_assembly_file=GENOME_ASSEMBLY,
            annotation_file=combine.outputs['gff'],
            working_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    return gwf