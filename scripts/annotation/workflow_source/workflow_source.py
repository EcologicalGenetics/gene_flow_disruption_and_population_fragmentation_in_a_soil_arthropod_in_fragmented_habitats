#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def braker3_annotation_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Annotates genome assembly using BRAKER3.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    CONFIG = yaml.safe_load(open(config_file))
    ACCOUNT: str = CONFIG['account']
    SPECIES_NAME: str = CONFIG['species_name']
    OUTPUT_DIR: str = CONFIG['output_directory_path']
    GENOME_ASSEMBLY: str = CONFIG['genome_assembly_file']
    BRAKER_SETTING: int = CONFIG['braker_setting']
    RNA_READS: list = CONFIG['rna_sequence_files']
    PROTEIN_DB: str = CONFIG['protein_database_file']
    REFERENCE: str = CONFIG['referene_genome_file']
    GTF: str = CONFIG['gtf_genome_annotation_file']
    ALT_NAME: str = CONFIG['alternative_folder_name']
    BUSCO: str = CONFIG['busco_dataset_name']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    top_dir = f'{OUTPUT_DIR}/{ALT_NAME.replace(" ", "_")}/annotation' if ALT_NAME else f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}/annotation'
    # os.makedirs(top_dir, exist_ok=True)

    if BRAKER_SETTING == 1 or BRAKER_SETTING == 3:
        index = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_STAR_index',
            template=star_index(
                genome_assembly_file=GENOME_ASSEMBLY,
                output_directory=top_dir
            )
        )
        
        rna_alignment = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_STAR_alignment',
            template=star_alignment(
                rna_sequence_files=RNA_READS,
                star_index_directory=f'{top_dir}/indices',
                output_directory=top_dir,
                species_name=SPECIES_NAME
            )
        )

    # Running BRAKER using only RNA evidence
    if BRAKER_SETTING == 1:
        braker = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_braker1',
            template=braker1(
                genome_assembly_file=GENOME_ASSEMBLY,
                rna_alignment_bam=rna_alignment.outputs['bam'],
                output_directory=top_dir,
                species_name=ALT_NAME if ALT_NAME else SPECIES_NAME
            )
        )

    # Running BRAKER using only protein evidence
    elif BRAKER_SETTING == 2:
        if not PROTEIN_DB:
            make_database = gwf.target_from_template(
                name=f'Create_protein_database',
                template=make_protein_db(
                    reference_genome_file=REFERENCE,
                    gtf_annotation_file=GTF,
                    output_directory=OUTPUT_DIR
                )
            )
            PROTEIN_DB = make_database.outputs['proteins']

        braker = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_braker2',
            template=braker2(
                genome_assembly_file=GENOME_ASSEMBLY,
                protein_database_file=PROTEIN_DB,
                output_directory=top_dir,
                species_name=ALT_NAME if ALT_NAME else SPECIES_NAME
            )
        )

    # Running BRAKER using both RNA and protein evidence
    elif BRAKER_SETTING == 3:
        if not PROTEIN_DB:
            make_database = gwf.target_from_template(
                name=f'Create_protein_database',
                template=make_protein_db(
                    reference_genome_file=REFERENCE,
                    gtf_annotation_file=GTF,
                    output_directory=OUTPUT_DIR
                )
            )
            PROTEIN_DB = make_database.outputs['proteins']

        braker = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_braker3',
            template=braker3(
                genome_assembly_file=GENOME_ASSEMBLY,
                rna_alignment_bam=rna_alignment.outputs['bam'],
                protein_database_file=PROTEIN_DB,
                output_directory=top_dir,
                species_name=ALT_NAME if ALT_NAME else SPECIES_NAME
            )
        )
    # elif DATABASE == 2:
        
    #     braker_db_exists_not = gwf.target_from_template(
    #         name=f'{species_abbreviation(SPECIES_NAME)}_braker3',
    #         template=braker3(
    #             genome_assembly_file=GENOME_ASSEMBLY,
    #             rna_alignment_bam=rna_alignment.outputs['bam'],
    #             protein_database_file=make_database.outputs['proteins'],
    #             output_directory=top_dir,
    #             species_name=SPECIES_NAME
    #         )
    #     )

    busco = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_busco_protein',
        template=busco_protein(
            genome_assembly_file=GENOME_ASSEMBLY,
            genome_annotation_gtf=braker.outputs['gtf'],
            output_directory=top_dir,
            busco_dataset=BUSCO
        )
    )

    return gwf