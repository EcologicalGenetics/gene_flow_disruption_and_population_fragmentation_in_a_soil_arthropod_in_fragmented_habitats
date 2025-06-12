#!/bin/env python3
from gwf import AnonymousTarget
import os, glob
from datetime import date

########################## Functions ##########################

def species_abbreviation(species_name: str) -> str:
	"""Creates species abbreviation from species name.

	:param str species_name:
		Species name written as *genus* *species*"""
	genus, species = species_name.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

########################## Braker3 ##########################

def star_index(genome_assembly_file: str, output_directory: str):
	"""
	Template: Indexes genome assembly prior to :script:`STAR` alignment of RNA data.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': genome_assembly_file}
	outputs = {'indices': [f'{output_directory}/indices/chrLength.txt',
						   f'{output_directory}/indices/chrName.txt',
						   f'{output_directory}/indices/chrNameLength.txt',
						   f'{output_directory}/indices/chrStart.txt',
						   f'{output_directory}/indices/Genome',
						   f'{output_directory}/indices/genomeParameters.txt',
						   f'{output_directory}/indices/SA',
						   f'{output_directory}/indices/SAindex'],
				'log': f'{output_directory}/indices/Log.out'}
	options = {
		'cores': 30,
		'memory': '40g',
		'walltime': '02:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate annotation
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/indices ] || mkdir -p {output_directory}/indices
	[ -d {output_directory}/tmp ] && rm -rf {output_directory}/tmp
	[ "$(ls -A {output_directory}/indices)" ] || rm -rf {output_directory}/indices/*
	cd {output_directory}
	
	salength=$(awk 'BEGIN{{genome_length = 0}} {{if ($0 ~ /[^>]/) {{genome_length += length($0) - 1}}}} END{{alt = int(((log(genome_length) / log(2)) / 2) - 1); if (alt < 14) {{print alt}} else {{print 14}}}}' {genome_assembly_file})

	STAR \
		--runThreadN {options['cores']} \
		--runMode genomeGenerate \
		--genomeDir {output_directory}/indices \
		--genomeFastaFiles {genome_assembly_file} \
		--outTmpDir {output_directory}/tmp \
		--genomeSAindexNbases "$salength"
	
	mv {output_directory}/Log.out {output_directory}/indices/Log.out

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def star_alignment(rna_sequence_files: list, star_index_directory: str, output_directory: str, species_name: str):
	"""
	Template: Align RNA data to indexed genome assembly using :script:`STAR`
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	working_dir = output_directory
	inputs = {'rna': rna_sequence_files,
		   	  'indices': [f'{output_directory}/indices/chrLength.txt',
						  f'{output_directory}/indices/chrName.txt',
						  f'{output_directory}/indices/chrNameLength.txt',
						  f'{output_directory}/indices/chrStart.txt',
						  f'{output_directory}/indices/Genome',
						  f'{output_directory}/indices/genomeParameters.txt',
						  f'{output_directory}/indices/SA',
						  f'{output_directory}/indices/SAindex']}
	outputs = {'bam': f'{output_directory}/rna_alignment/{species_abbreviation(species_name)}_Aligned.sortedByCoord.out.bam',
			   'sj': f'{output_directory}/rna_alignment/{species_abbreviation(species_name)}_SJ.out.tab',
			   'logs': [f'{output_directory}/rna_alignment/{species_abbreviation(species_name)}_Log.final.out',
			   			f'{output_directory}/rna_alignment/{species_abbreviation(species_name)}_Log.out',
						f'{output_directory}/rna_alignment/{species_abbreviation(species_name)}_Log.progress.out']}
	options = {
		'cores': 30,
		'memory': '40g',
		'walltime': '04:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate annotation
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/rna_alignment ] || mkdir -p {output_directory}/rna_alignment
	[ -d {output_directory}/tmp ] && rm -rf {output_directory}/tmp
	
	STAR \
		--runThreadN {options['cores']} \
		--runMode alignReads \
		--genomeDir {star_index_directory} \
		--readFilesIn {" ".join(rna_sequence_files)} \
		--readFilesCommand zcat \
		--outFileNamePrefix {output_directory}/rna_alignment/{species_abbreviation(species_name)}_ \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMstrandField intronMotif \
		--outSAMattrRGline ID:{species_abbreviation(species_name)}_RNA SM:{species_abbreviation(species_name)}_RNA \
		--outTmpDir {output_directory}/tmp \
		--outFilterScoreMinOverLread 0 \
		--outFilterMatchNminOverLread 0 \
		--outFilterMatchNmin 0
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, options=options, spec=spec)

def make_protein_db(reference_genome_file: str, gtf_annotation_file: str, output_directory: str):
	"""
	Template: Make protein database file in :format:`FASTA` format from reference genome and :format:`GTF` annotation file.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'reference': reference_genome_file,
		   	  'gtf':gtf_annotation_file}
	outputs = {'proteins': f'{output_directory}/proteinDB/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.protein.fasta'}
	protect = outputs['proteins']
	options = {
		'cores': 1,
		'memory': '20g',
		'walltime': '06:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/proteinDB ] || mkdir -p {output_directory}/proteinDB
	
	cd {output_directory}/proteinDB

	ln -s {reference_genome_file} {os.path.basename(reference_genome_file)}

	agat_sp_extract_sequences.pl \
		--gff {gtf_annotation_file} \
		--fasta {os.path.basename(reference_genome_file)} \
		--type cds \
		--protein \
		--output {output_directory}/proteinDB/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.protein.initial.prog.fasta
	
	awk \
		'{{gsub(/\s$/, "");
		print $0}}' \
		{output_directory}/proteinDB/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.protein.initial.prog.fasta \
		> {output_directory}/proteinDB/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.protein.prog.fasta
	
	rm {os.path.basename(reference_genome_file)}
	rm {os.path.basename(reference_genome_file)}.index
	rm {output_directory}/proteinDB/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.protein.initial.prog.fasta
	mv {output_directory}/proteinDB/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.protein.prog.fasta {outputs['proteins']}
	rm *.agat.log
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

# awk \
# 		'{{if ($0 ~ /^>/)
# 			{{gsub(/ /, ";")}}
# 		gsub(/\s$/, "");
# 		print $0}}' \
# 		{output_directory}/proteinDB/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.protein.initial.prog.fasta \
# 		> {output_directory}/proteinDB/{os.path.splitext(os.path.basename(reference_genome_file))[0]}.protein.prog.fasta

def braker1(genome_assembly_file: str, rna_alignment_bam: str, output_directory: str, species_name: str, genemark: str = '/home/jepe/software/GeneMark-ETP/bin/', prothint: str ='/home/jepe/software/ProtHint-2.6.0/bin'):
	"""
	Template: Runs BRAKER3 using RNA-sequence data to predict genes function and annotate genome assembly.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': genome_assembly_file,
		   	  'rna': rna_alignment_bam}
	outputs = {'gtf': f'{output_directory}/braker1/braker.gtf',
			   'bed': f'{output_directory}/braker1/genes.bed',
			   'coding': f'{output_directory}/braker1/braker.codingseq',
			   'aa': f'{output_directory}/braker1/braker.aa',
			   'evidence': f'{output_directory}/braker1/hintsfile.gff',
			   'cite': f'{output_directory}/braker1/what-to-cite.txt',
			   'other': [f'{output_directory}/braker1/genome_header.map',
						 f'{output_directory}/braker1/braker.log']}
	protect = [outputs['gtf'], outputs['coding'], outputs['aa'], outputs['evidence']]
	options = {
		'cores': 30,
		'memory': '300g',
		'walltime': '10:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate annotation
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	if [ -d {output_directory}/braker1 ]; then
		rm -rf {output_directory}/braker1
		mkdir -p {output_directory}/braker1
	else
		mkdir -p {output_directory}/braker1
	fi
	
	braker.pl \
		--genome {genome_assembly_file} \
		--bam {rna_alignment_bam} \
		--species {species_name.replace(' ', '_')}_{date.today().day}-{date.today().month}-{date.today().year} \
		--threads {options['cores']} \
		--workingdir {output_directory}/braker1 \
		--GENEMARK_PATH {genemark} \
		--PROTHINT_PATH {prothint}
	
	awk \
		'BEGIN{{FS=OFS="\\t"}} 
		{{if ($3 == "gene") 
			{{print $1, $4-1, $5}}}}' \
		{outputs['gtf']} \
		> {output_directory}/braker1/genes.prog.bed
	
	mv {output_directory}/braker1/genes.prog.bed {outputs['bed']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def braker2(genome_assembly_file: str, protein_database_file: str, output_directory: str, species_name: str, genemark: str = '/home/jepe/software/GeneMark-ETP/bin/', prothint: str ='/home/jepe/software/ProtHint-2.6.0/bin'):
	"""
	Template: Runs BRAKER3 using a protein database to predict genes function and annotate genome assembly.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': genome_assembly_file,
			  'protein': protein_database_file}
	outputs = {'gtf': f'{output_directory}/braker2/braker.gtf',
			   'bed': f'{output_directory}/braker2/genes.bed',
			   'coding': f'{output_directory}/braker2/braker.codingseq',
			   'aa': f'{output_directory}/braker2/braker.aa',
			   'evidence': f'{output_directory}/braker2/hintsfile.gff',
			   'cite': f'{output_directory}/braker2/what-to-cite.txt',
			   'other': [f'{output_directory}/braker2/genome_header.map',
						 f'{output_directory}/braker2/braker.log']}
	protect = [outputs['gtf'], outputs['coding'], outputs['aa'], outputs['evidence']]
	options = {
		'cores': 30,
		'memory': '300g',
		'walltime': '10:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate annotation
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	if [ -d {output_directory}/braker2 ]; then
		rm -rf {output_directory}/braker2
		mkdir -p {output_directory}/braker2
	else
		mkdir -p {output_directory}/braker2
	fi
	
	[ -e  ]

	braker.pl \
		--genome {genome_assembly_file} \
		--prot_seq {protein_database_file} \
		--species {species_name.replace(' ', '_')}_{date.today().day}-{date.today().month}-{date.today().year} \
		--threads {options['cores']} \
		--workingdir {output_directory}/braker2 \
		--GENEMARK_PATH {genemark} \
		--PROTHINT_PATH {prothint}
	
	awk \
		'BEGIN{{FS=OFS="\\t"}} 
		{{if ($3 == "gene") 
			{{print $1, $4-1, $5}}}}' \
		{outputs['gtf']} \
		> {output_directory}/braker2/genes.prog.bed
	
	mv {output_directory}/braker2/genes.prog.bed {outputs['bed']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def braker3(genome_assembly_file: str, rna_alignment_bam: str, protein_database_file: str, output_directory: str, species_name: str, genemark: str = '/home/jepe/software/GeneMark-ETP/bin/', prothint: str ='/home/jepe/software/ProtHint-2.6.0/bin'):
	"""
	Template: Runs BRAKER3 using both RNA-sequence data and a protein database to predict genes function and annotate genome assembly.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': genome_assembly_file,
		   	  'rna': rna_alignment_bam,
			  'protein': protein_database_file}
	outputs = {'gtf': f'{output_directory}/braker3/braker.gtf',
			   'bed': f'{output_directory}/braker3/genes.bed',
			   'coding': f'{output_directory}/braker3/braker.codingseq',
			   'aa': f'{output_directory}/braker3/braker.aa',
			   'evidence': f'{output_directory}/braker3/hintsfile.gff',
			   'cite': f'{output_directory}/braker3/what-to-cite.txt',
			   'other': [f'{output_directory}/braker3/genome_header.map',
						 f'{output_directory}/braker3/braker.log']}
	protect = [outputs['gtf'], outputs['coding'], outputs['aa'], outputs['evidence']]
	options = {
		'cores': 30,
		'memory': '300g',
		'walltime': '10:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate annotation
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	if [ -d {output_directory}/braker3 ]; then
		rm -rf {output_directory}/braker3
		mkdir -p {output_directory}/braker3
	else
		mkdir -p {output_directory}/braker3
	fi
	
	braker.pl \
		--genome {genome_assembly_file} \
		--bam {rna_alignment_bam} \
		--prot_seq {protein_database_file} \
		--species {species_name.replace(' ', '_')}_{date.today().day}-{date.today().month}-{date.today().year} \
		--threads {options['cores']} \
		--workingdir {output_directory}/braker3 \
		--GENEMARK_PATH {genemark} \
		--PROTHINT_PATH {prothint}
	
	awk \
		'BEGIN{{FS=OFS="\\t"}} 
		{{if ($3 == "gene") 
			{{print $1, $4-1, $5}}}}' \
		{outputs['gtf']} \
		> {output_directory}/braker3/genes.prog.bed
	
	mv {output_directory}/braker3/genes.prog.bed {outputs['bed']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def busco_protein(genome_assembly_file: str, genome_annotation_gtf: str, busco_dataset: str, output_directory: str, busco_download_path: str = '/faststorage/project/EcoGenetics/databases/BUSCO'):
	"""
	Template: Runs BUSCO analysis on protein sequences from an annotated gene set.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': genome_assembly_file,
		   	  'annotation': genome_annotation_gtf}
	outputs = {'protein': f'{output_directory}/busco/{os.path.basename(genome_assembly_file)}.protein.fa',
			   'stattxt': f'{output_directory}/busco/busco_{os.path.basename(genome_assembly_file)}.protein.fa/short_summary.specific.{busco_dataset}.busco_{os.path.basename(genome_assembly_file)}.protein.fa.txt',
			   'statjson': f'{output_directory}/busco/busco_{os.path.basename(genome_assembly_file)}.protein.fa/short_summary.specific.{busco_dataset}.busco_{os.path.basename(genome_assembly_file)}.protein.fa.json'}
	protect = [outputs['stattxt'], outputs['statjson']]
	options = {
		'cores': 30,
		'memory': '50g',
		'walltime': '10:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	mkdir -p {output_directory}/busco

	agat_sp_extract_sequences.pl \
		--fasta {genome_assembly_file} \
		--gff {genome_annotation_gtf} \
		--type gene \
		--protein \
		--output {output_directory}/busco/{os.path.basename(genome_assembly_file)}.protein.fa \

	busco \
		--cpu {options['cores']} \
		--force \
		--in {output_directory}/busco/{os.path.basename(genome_assembly_file)}.protein.fa \
		--mode proteins \
		--out busco_{os.path.basename(genome_assembly_file)}.protein.fa \
		--out_path {output_directory}/busco \
		--download_path {busco_download_path} \
		--lineage {busco_download_path}/lineages/{busco_dataset} \
		--tar
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)