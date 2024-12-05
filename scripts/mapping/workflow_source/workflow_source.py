#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def mapping_resequencing_data_population_genetics_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: Align resequencing data to reference genome and do basic filtering.
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(config_file))
	ACCOUNT: str = CONFIG['account']
	TAXONOMY: str = CONFIG['taxonomic_group'].lower()
	SPECIES_NAME: str = CONFIG['species_name']
	REFERENCE_GENOME: str = CONFIG['reference_genome_path']
	WORK_DIR: str = CONFIG['working_directory_path'][:len(CONFIG['working_directory_path']) - 1] if CONFIG['working_directory_path'].endswith('/') else CONFIG['working_directory_path']
	OUTPUT_DIR: str = CONFIG['output_directory_path'][:len(CONFIG['output_directory_path']) - 1] if CONFIG['output_directory_path'].endswith('/') else CONFIG['output_directory_path']
	SAMPLE_LISTS: list = CONFIG['sample_lists']
	ADAPTERREMOVAL_SETTINGS: dict = CONFIG['adapterremoval_settings']
	AR_MIN_QUAL: int | None = ADAPTERREMOVAL_SETTINGS['min_quality'] if ADAPTERREMOVAL_SETTINGS['min_quality'] else 25
	AR_MIN_LENGTH: int | None = ADAPTERREMOVAL_SETTINGS['min_length'] if ADAPTERREMOVAL_SETTINGS['min_length'] else 20
	AR_SEQUENCE1: str | None = ADAPTERREMOVAL_SETTINGS['adaptersequence1'] if ADAPTERREMOVAL_SETTINGS['adaptersequence1'] else 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA'
	AR_SEQUENCE2: str | None = ADAPTERREMOVAL_SETTINGS['adaptersequence2'] if ADAPTERREMOVAL_SETTINGS['adaptersequence2'] else 'AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG'
	SAMTOOLS_FILTER_SETTINGS: dict = CONFIG['filter_settings']
	STF_EXCLUDE: int | None = SAMTOOLS_FILTER_SETTINGS['flags_excluded'] if SAMTOOLS_FILTER_SETTINGS['flags_excluded'] else 3844
	STF_REQUIRED: int | None = SAMTOOLS_FILTER_SETTINGS['flags_required'] if SAMTOOLS_FILTER_SETTINGS['flags_required'] else None
	STF_MIN_MQ: int | None = SAMTOOLS_FILTER_SETTINGS['min_mq'] if SAMTOOLS_FILTER_SETTINGS['min_mq'] else 20
	
	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/{TAXONOMY.replace(" ", "_")}/mapping'
	top_out = f'{OUTPUT_DIR}/alignments/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}'
	
	index = gwf.target_from_template(
		name=f'{SPECIES_NAME.replace(" ", "_")}_index_reference',
		template=index_reference_genome(
			reference_genome_file=REFERENCE_GENOME,
			output_directory=top_dir)
	)

	between_group_qualimap = []

	for GROUP in SAMPLE_LISTS:
		if not GROUP['sample_folder_list']:
			continue
		GROUP_NAME: str = GROUP['group_name'].lower()
		SAMPLE_FOLDER_LIST: list = GROUP['sample_folder_list']
		sample_list = [get_sample_data(path, data_type=2) for path in SAMPLE_FOLDER_LIST]
		within_group_qualimap = []
		for sample in sample_list:	
			adapterremoval = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_adapterremoval',
				template=adapterremoval_pairedend(
					sample_name=sample['sample_name'],
					read1_files=sample['read1_files'],
					read2_files=sample['read2_files'],
					output_directory=f'{top_dir}/{GROUP_NAME}',
					min_quality=AR_MIN_QUAL,
					min_length=AR_MIN_LENGTH,
					adapter1=AR_SEQUENCE1,
					adapter2=AR_SEQUENCE2
				)
			)

			align_paired = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_paired_alignment',
				template=alignment_pairedend(
					read1_file=adapterremoval.outputs['pair1'],
					read2_file=adapterremoval.outputs['pair2'],
					reference_genome_file=index.outputs['symlink'],
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}'
				)
			)

			align_collapsed = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_collapsed_alignment',
				template=alignment_collapsed(
					collapsed_read_files=adapterremoval.outputs['collapsed'],
					reference_genome_file=index.outputs['symlink'],
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}'
				)
			)

			merge = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_merge_alignments',
				template=merge_alignments(
					alignment_files=[align_paired.outputs['alignment'],
									 align_collapsed.outputs['alignment']],
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}'
				)
			)

			mark_duplicates = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_markdup',
				template=mark_duplicates_samtools(
					alignment_file=merge.outputs['merged'],
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}'
				)
			)

			unmapped = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_extract_unmapped_reads',
				template=extract_unmapped_reads(
					alignment_file=mark_duplicates.outputs['markdup'],
					sample_name=sample['sample_name'],
					output_directory=f'{top_out}/{GROUP_NAME}' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/unmapped'
				)
			)

			pre_filter_stats = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_pre_filter_stats',
				template=samtools_stats(
					alignment_file=mark_duplicates.outputs['markdup'],
					output_directory=f'{top_out}/{GROUP_NAME}/{sample['sample_name']}/qc/pre_filtering_stats' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/pre_filtering_stats/{sample["sample_name"]}'
				)
			)

			filter_alignment = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_filtering',
				template=samtools_filter(
					alignment_file=mark_duplicates.outputs['markdup'],
					sample_name=sample['sample_name'],
					output_directory=f'{top_out}/{GROUP_NAME}' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/filtered_alignment',
					flags_excluded=STF_EXCLUDE,
					flags_required=STF_REQUIRED,
					min_mq=STF_MIN_MQ
				)
			)

			post_filter_stats = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_post_filter_stats',
				template=samtools_stats(
					alignment_file=filter_alignment.outputs['filtered'],
					output_directory=f'{top_out}/{GROUP_NAME}/{sample['sample_name']}/qc/pre_filtering_stats' if OUTPUT_DIR else f'{top_dir}/{GROUP_NAME}/post_filtering_stats/{sample["sample_name"]}'
				)
			)

			qualimap = gwf.target_from_template(
				name=f'{sample["sample_name"].replace("-", "_")}_qualimap',
				template=qc_qualimap(
					alignment_file=filter_alignment.outputs['filtered'],
					output_directory=f'{top_out}/{GROUP_NAME}/{sample['sample_name']}/qc/bamqc' if OUTPUT_DIR else f'{top_dir}/bamqc/individual/{GROUP_NAME}/{sample["sample_name"]}'
				)
			)

			within_group_qualimap.append([sample["sample_name"], os.path.dirname(qualimap.outputs["raw"]), sample["sample_name"]])
			between_group_qualimap.append([sample["sample_name"], os.path.dirname(qualimap.outputs["raw"]), GROUP_NAME])

		within_multi_qualimap = gwf.target_from_template(
			name=f'{GROUP_NAME}_multi_qualimap',
			template=qualimap_multi(
				dataset=within_group_qualimap,
				output_directory=f'{top_out}/{GROUP_NAME}/group_bamqc' if OUTPUT_DIR else f'{top_dir}/bamqc/group/{GROUP_NAME}'
			)
		)
	
	between_multi_qualimap = gwf.target_from_template(
		name=f'all_groups_multi_qualimap',
		template=qualimap_multi(
			dataset=between_group_qualimap,
			output_directory=f'{top_out}/species_bamqc' if OUTPUT_DIR else f'{top_dir}/bamqc/all'
		)
	)

	return gwf