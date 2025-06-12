#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def freebayes_vcf_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Create :format:`VCF` using :script:`freebayes`.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    SAMPLE_LIST: list = config['sample_list']
    REFERENCE_GENOME: str = config['reference_genome_path']
    OUTPUT_DIR: str = config['output_directory_path']
    WORK_DIR: str = config['working_directory_path']
    PARTITION_SIZE: int = config['partition_size']
    PLOIDY: int = config['sample_ploidy']
    BESTN: int = config['best_n_alleles']
    ALT_FRACTION: float | int = config['min_alternate_fraction']
    ALT_COUNT: int = config['min_alternate_count']

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    if os.path.exists(f'reference_partitions.{PARTITION_SIZE}bp.txt') and os.path.exists('reference_sequences.txt'):
        # If files exists reads data directly from files
        # Loads reference genome partitioning
        with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'r') as infile:
            partitions = [{'num': entry.split(sep='\t')[0].strip(), 'region': entry.split(sep='\t')[1].strip(), 'start': entry.split(sep='\t')[2].strip(), 'end': entry.split(sep='\t')[3].strip()} for entry in infile]
            npadding = len(str(sum(1 for line in partitions)))
        # Loads list of contigs in reference genome
        with open('reference_sequences.txt', 'r') as infile:
            contigs = [{'contig': entry.split(sep='\t')[0].strip()} for entry in infile]
    else:
        # If files don't exist, generate data and write files
        sequences = parse_fasta(REFERENCE_GENOME)
        with open(f'reference_sequences.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
        # Partitions reference genome
        npadding = padding_calculator(parse_fasta=sequences, size=PARTITION_SIZE)
        partitions = partition_chrom(parse_fasta=sequences, size=PARTITION_SIZE, npad=npadding)
        with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
        # Creates list of contigs in reference genome
        contigs = [{'contig': contig['sequence_name']} for contig in sequences]

    top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/freebayes'
    output_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}/freebayes'

    freebayes_parts = gwf.map(
        template_func=freebayes_chrom,
        inputs=partitions,
        name=name_freebayes_chrom,
        extra={'reference_genome_file': REFERENCE_GENOME,
               'bam_file_list': SAMPLE_LIST,
               'output_directory': top_dir,
               'species_name': SPECIES_NAME,
               'ploidy': PLOIDY,
               'best_n_alleles': BESTN,
               'min_alternate_fraction': ALT_FRACTION,
               'min_alternate_count': ALT_COUNT}
        )
    
    concat_freebayes = gwf.target_from_template(
            name='concatenate_freebayes',
            template=concat_vcf(
                files=collect(freebayes_parts.outputs, ['vcf'])['vcfs'],
                output_name=f'{species_abbreviation(SPECIES_NAME)}.freebayes_n{BESTN}_p{PLOIDY}_minaltfrc{ALT_FRACTION}_minaltcnt{ALT_COUNT}',
                output_directory=output_dir,
                compress=True
            )
        )

    return gwf

def freebayes_population_set_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Create :format:`VCF` file for each sample in configuration.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    CONFIG = yaml.safe_load(open(config_file))
    ACCOUNT: str = CONFIG['account']
    SPECIES_NAME: str = CONFIG['species_name']
    REFERENCE_GENOME: str = CONFIG['reference_genome_path']
    WORK_DIR: str = CONFIG['working_directory_path']
    OUTPUT_DIR: str = CONFIG['output_directory_path']
    PARTITION_SIZE: int = CONFIG['partition_size']
    FREEBAYES_SETTINGS: dict = CONFIG['freebayes_settings']
    FREEBAYES_PLOIDY: int = FREEBAYES_SETTINGS['sample_ploidy']
    FREEBAYES_BESTN: int = FREEBAYES_SETTINGS['best_n_alleles']
    FREEBAYES_MINALTFRC: float | int = FREEBAYES_SETTINGS['min_alternate_fraction']
    FREEBAYES_MINALTCNT: int = FREEBAYES_SETTINGS['min_alternate_count']
    FILTERING: dict = CONFIG['filtering']
    FILTERING_MINDP: int = FILTERING['minimum_depth']
    FILTERING_MAXDP: int = FILTERING['maximum_depth']
    SAMPLE_LIST: list = CONFIG['sample_list']

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    if os.path.exists(f'reference_partitions.{PARTITION_SIZE}bp.txt') and os.path.exists('reference_sequences.txt'):
        # If files exists reads data directly from files
        # Loads reference genome partitioning
        with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'r') as infile:
            partitions = [{'num': entry.split(sep='\t')[0].strip(), 'region': entry.split(sep='\t')[1].strip(), 'start': entry.split(sep='\t')[2].strip(), 'end': entry.split(sep='\t')[3].strip()} for entry in infile]
            npadding = len(str(sum(1 for line in partitions)))
        # Loads list of contigs in reference genome
        with open('reference_sequences.txt', 'r') as infile:
            contigs = [{'contig': entry.split(sep='\t')[0].strip()} for entry in infile]
    else:
        # If files don't exist, generate data and write files
        sequences = parse_fasta(REFERENCE_GENOME)
        with open(f'reference_sequences.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
        # Partitions reference genome
        npadding = padding_calculator(parse_fasta=sequences, size=PARTITION_SIZE)
        partitions = partition_chrom(parse_fasta=sequences, size=PARTITION_SIZE, npad=npadding)
        with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
        # Creates list of contigs in reference genome
        contigs = [{'contig': contig['sequence_name']} for contig in sequences]

    top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/vcf'

    for sample in SAMPLE_LIST:
        SAMPLE_NAME: str = sample['sample_name']
        SAMPLE_GROUP: str = sample['sample_group'].lower()
        BAM: str = sample['bam_file']

        freebayes_parts = gwf.map(
            name=name_freebayes_chrom_single_sample,
            template_func=freebayes_chrom_single_sample,
            inputs=partitions,
            extra={'reference_genome_file': REFERENCE_GENOME,
                   'bam_file': BAM,
                   'output_directory': top_dir,
                   'sample_group': SAMPLE_GROUP,
                   'sample_name': SAMPLE_NAME,
                   'ploidy': FREEBAYES_PLOIDY,
                   'best_n_alleles': FREEBAYES_BESTN,
                   'min_alternate_fraction': FREEBAYES_MINALTFRC,
                   'min_alternate_count': FREEBAYES_MINALTCNT}
        )

        concat_freebayes = gwf.target_from_template(
                name=f'concatenate_freebayes_vcf_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
                template=concat_vcf(
                    files=collect(freebayes_parts.outputs, ['vcf'])['vcfs'],
                    output_name=f'{SAMPLE_NAME}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}',
                    output_directory=f'{top_dir}/raw_vcf/{SAMPLE_GROUP}/{SAMPLE_NAME}',
                    compress=True
                )
        )

        filtering = gwf.target_from_template(
            name=f'filter_vcf_{SAMPLE_GROUP}_{SAMPLE_NAME.replace("-", "_")}',
            template=filter_vcf(
                vcf_file=concat_freebayes.outputs['concat_file'],
                output_directory=f'{OUTPUT_DIR}' if OUTPUT_DIR else f'{top_dir}/filtered_vcf',
                sample_group=SAMPLE_GROUP,
                sample_name=SAMPLE_NAME,
                min_depth=FILTERING_MINDP,
                max_depth=FILTERING_MAXDP
            )
        )
    
    return gwf
