#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def yahs_hic_scaffolding_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Scaffolds draft genome assembly using Hi-C data
    
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
    HIC_READS: list = CONFIG['hic_sequence_files']
    DRAFT_GENOME: str = CONFIG['draft_genome_file']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    top_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}'
    os.makedirs(top_dir, exist_ok=True)

    index = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_index',
        template=index_reference(
            reference_genome_file=DRAFT_GENOME,
            output_directory=top_dir
        )
    )
    
    hic_align = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_HiC_alignment',
        template=hic_alignment_to_draft_assembly(
            hic_seqeuence_files=HIC_READS,
            draft_genome_file=DRAFT_GENOME,
            reference_indices=index.outputs['bwa'],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    mark_duplicates = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_mark_duplicates',
        template=mark_duplicates_picard(
            alignment_bam_file=hic_align.outputs['bam'],
            output_directory=top_dir
        )
    )

    hic_scaffolding = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_HiC_scaffolding',
        template=hic_scaffolding_yahs(
            draft_genome_file=DRAFT_GENOME,
            hic_to_draft_bam_file=mark_duplicates.outputs['markdup'],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    conversion = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_juicer_conversion',
        template=yahs_conversion_manual_curation(
            hic_bin_file=hic_scaffolding.outputs['bin'],
            scaffolds_final_agp_file=hic_scaffolding.outputs['final_agp'],
            draft_assembly_fai_index_file=index.outputs['fai'],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    matrix = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_matrix',
        template=contact_matrix_manual_curation(
            JBAT_text_file=conversion.outputs['jbat'][0],
            JBAT_log_file=conversion.outputs['jbat'][4],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    return gwf

def juicer_hic_scaffolding_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Scaffolds draft genome assembly using Hi-C data with the :script:`juicer` pipeline.
    
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
    HIC_READS: dict = CONFIG['hic_sequence_files']
    HIC_READ1: list = HIC_READS['read1']
    HIC_READ2: list = HIC_READS['read2']
    DRAFT_GENOME: str = CONFIG['draft_genome_file']
    SETTINGS_3DDNA: dict = CONFIG['3ddna']
    INPUT_SIZE: int = SETTINGS_3DDNA['input_size']
    EDIT_ROUDNS: int = SETTINGS_3DDNA['number_of_edit_rounds']
    REVIEW_FILE: str = CONFIG['reviewed_assembly_file']
    POST_EDITING: dict = CONFIG['post_editing']
    CHROM_NUM: int = POST_EDITING['number_of_chromosomes_to_keep']
    BUNDLE: int = POST_EDITING['bundle_debris'] if POST_EDITING['bundle_debris'] else 0
    INSERTION_SIZE: int = POST_EDITING['insertion_size'] if POST_EDITING['insertion_size'] else 1000
    BUSCO: str = CONFIG['busco_dataset_name']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    top_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}/HiC_scaffolding'

    setup = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_juicer_setup',
        template=setup_for_juicer(
            hic_read1_files=HIC_READ1,
            hic_read2_files=HIC_READ2,
            draft_genome_assembly_file=DRAFT_GENOME,
            output_directory=top_dir
        )
    )
    
    juicer = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_juicer_pipeline_start',
        template=juicer_pipeline(
            draft_genome_assembly_file=setup.outputs['reference'],
            chrom_sizes_file=setup.outputs['sizes'],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )

    assembly = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_3d_dna',
        template=assembly_3ddna(
            draft_assembly_file=DRAFT_GENOME,
            merged_nodups_file=juicer.outputs['nodups'],
            output_directory=top_dir,
            edit_rounds=EDIT_ROUDNS,
            input_size=INPUT_SIZE
        )
    )

    if REVIEW_FILE:
        finalize = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_finalize_assembly',
            template=finalize_3ddna(
                reviewed_assembly_file=REVIEW_FILE,
                draft_assembly_file=DRAFT_GENOME,
                merged_nodups_file=juicer.outputs['nodups'],
                number_of_chromosomes=CHROM_NUM,
                final_hic_file=assembly.outputs['final_hic']
            )
        )

        if BUNDLE:
            bundle_debris = gwf.target_from_template(
                name=f'bundle_debris',
                template=bundle_sequences(
                    assembly_fasta_file=finalize.outputs['fasta'],
                    n_chromosomes_to_keep=CHROM_NUM,
                    n_insertion_size=INSERTION_SIZE,
                )
            )

            busco_bundle = gwf.target_from_template(
                name=f'BUSCO_assembly_bundle',
                template=busco_genome(
                    genome_assembly_file=bundle_debris.outputs['bundled'],
                    busco_dataset=BUSCO
                )
            )

        busco_no_debris = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_BUSCO_assembly_no_debris',
            template=busco_genome(
                genome_assembly_file=finalize.outputs['final_fasta'],
                busco_dataset=BUSCO
            )
        )

        busco_with_debris = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_BUSCO_assembly_w_debris',
            template=busco_genome(
                genome_assembly_file=finalize.outputs['fasta'],
                busco_dataset=BUSCO
            )
        )

    return gwf