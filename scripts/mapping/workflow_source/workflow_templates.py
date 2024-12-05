#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

########################## Functions ##########################

def species_abbreviation(species_name: str) -> str:
	"""Creates species abbreviation from species name.

	:param str species_name:
		Species name written as *genus* *species*"""
	genus, species = species_name.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

def get_sample_data(path: str, data_type: int) -> dict:
	"""Create dictionary of lists each containing all resequencing for for a sample in one read direction.
	
	:param str path:
		Path to sample directory containing resequencing data.
	:param int data_type:
		Sequencing type, single-end or paired-end (1 | 2)."""
	extformat = ('.fq', '.fa', '.fasta', '.fastq', '.fn')
	compress = ('', '.gz', '.gzip')
	if data_type == 1:
		reseq = {'sample_name': os.path.basename(path),
				 'read_files': sorted([os.path.join(path, file) for file in os.listdir(path) if file.endswith(tuple([i+j for i in extformat for j in compress]))])}
	elif data_type == 2:
		r1 = ('_1', '_R1')
		r2 = ('_2', '_R2')
		reseq = {'sample_name': os.path.basename(path),
				 'read1_files': sorted([os.path.join(path, file) for file in os.listdir(path) if file.endswith(tuple([i+j+k for i in r1 for j in extformat for k in compress]))]),
				 'read2_files': sorted([os.path.join(path, file) for file in os.listdir(path) if file.endswith(tuple([i+j+k for i in r2 for j in extformat for k in compress]))])}
	return reseq

########################## Mapping ##########################

def index_reference_genome(reference_genome_file: str, output_directory: str):
	"""
	Template: Index reference genome with :script:`bwa index` and :script:`samtools faidx`.
	
	Template I/O::
	
		inputs = {'reference': reference_genome_file}
		outputs = {'symlink': *, 'bwa': [*.amb, *.ann, *.pac, *.bwt, *.sa], 'fai': *.fai}
	
	:param str reference_genome_file:
		Path to reference genome.
	:param str output_directory:
		Path to output directory.
	"""
	inputs = {'reference': reference_genome_file}
	outputs = {'symlink': f'{output_directory}/reference/{os.path.basename(reference_genome_file)}',
			   'bwa': [f'{output_directory}/reference/{os.path.basename(reference_genome_file)}.amb',
					   f'{output_directory}/reference/{os.path.basename(reference_genome_file)}.ann',
					   f'{output_directory}/reference/{os.path.basename(reference_genome_file)}.pac',
					   f'{output_directory}/reference/{os.path.basename(reference_genome_file)}.bwt',
					   f'{output_directory}/reference/{os.path.basename(reference_genome_file)}.sa'],
			   'fai': f'{output_directory}/reference/{os.path.basename(reference_genome_file)}.fai'}
	protect = [outputs['symlink'], outputs['bwa'][0], outputs['bwa'][1], outputs['bwa'][2], outputs['bwa'][3], outputs['bwa'][4], outputs['fai']]
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '04:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/reference ] || mkdir -p {output_directory}/reference
	[ -e {output_directory}/reference/{os.path.basename(reference_genome_file)} ] && rm -f {output_directory}/reference/{os.path.basename(reference_genome_file)}
	ln -s {reference_genome_file} {output_directory}/reference/{os.path.basename(reference_genome_file)}
	
	bwa index \
		-p {output_directory}/reference/{os.path.basename(reference_genome_file)} \
		{output_directory}/reference/{os.path.basename(reference_genome_file)}
	
	samtools faidx \
		-o {output_directory}/reference/{os.path.basename(reference_genome_file)}.prog.fai \
		{output_directory}/reference/{os.path.basename(reference_genome_file)}
	
	mv {output_directory}/reference/{os.path.basename(reference_genome_file)}.prog.fai {outputs['fai']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def adapterremoval_pairedend(sample_name: str, read1_files: list, read2_files: list, output_directory: str, min_quality: int = 25, min_length: int = 20, adapter1: str = 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA', adapter2: str = 'AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG'):
	"""
	Template: Remove remnant adapter sequences from read data using :script:`AdapterRemoval`.
	
	Template I/O::
	
		inputs = {'read1': read1_files, 'read2': read2_files}
		outputs = {'pair1': *.pair1.truncated, 'pair2': *.pair2.truncated.,
			   'collapsed': [*.collapsed, *.collapsed.truncated],
			   'misc': [*.settings, *.singleton.truncated, *.discarded]}
	
	:param str sample_name:
		Name of sample.
	:param list read1_files:
		List of mate 1 read files.
	:param lsit read2_files:
		List of mate 2 read files.
	:param str output_directory:
		Path to output directory.
	:param int min_quality:
		Minimum quality score to include when trimming 5'/3' termini.
	:param int min_length:
		Minimum length of reads to keep after trimming.
	:param str adapter1:
		Adapter sequence expected to be found in mate 1 reads.
	:param str adapter2:
		Adapter sequence expected to be found in mate 2 reads.
	"""
	inputs = {'read1': read1_files,
		   	  'read2': read2_files}
	outputs = {'pair1': f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.pair1.truncated',
			   'pair2':f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.pair2.truncated',
			   'collapsed': [f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.collapsed',
			   				 f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.collapsed.truncated'],
			   'misc': [f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.settings',
			   			f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.singleton.truncated',
			   			f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.discarded']}
	options = {
		'cores': 16,
		'memory': '60g',
		'walltime': '06:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/adapterremoval/{sample_name} ] || mkdir -p {output_directory}/adapterremoval/{sample_name}
	
	AdapterRemoval \
		--threads {options['cores']} \
		--file1 {' '.join(read1_files)} \
		--file2 {' '.join(read2_files)} \
		--adapter1 {adapter1} \
		--adapter2 {adapter2} \
		--minquality {min_quality} \
		--minlength {min_length} \
		--basename {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog \
		--trimns \
		--trimqualities \
		--collapse
	
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.settings {outputs['misc'][0]}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.pair1.truncated {outputs['pair1']}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.pair2.truncated {outputs['pair2']}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.singleton.truncated {outputs['misc'][1]}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.collapsed {outputs['collapsed'][0]}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.collapsed.truncated {outputs['collapsed'][1]}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.discarded {outputs['misc'][2]}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def adapterremoval_singleend(sample_name: str, read_files: list, output_directory: str, min_qulaity: int = 25, min_length: int = 20, adapter1: str = 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA'):
	"""
	Template: Remove remnant adapter sequences from read data using :script:`AdapterRemoval`.
	
	Template I/O::
	
		inputs = {'read': read_files}
		outputs = {'truncated': *.pair1.truncated, 'misc': [*.settings, *.discarded]}
	
	:param str sample_name:
		Name of sample.
	:param list read_files:
		List of single-end read files.
	:param str output_directory:
		Path to output directory.
	:param int min_quality:
		Minimum quality score to include when trimming 5'/3' termini.
	:param int min_length:
		Minimum length of reads to keep after trimming.
	:param str adapter1:
		Adapter sequence expected to be found in reads.
	"""
	inputs = {'read': read_files}
	outputs = {'truncated': f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.pair1.truncated',
			   'misc': [f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.settings',
			   			f'{output_directory}/adapterremoval/{sample_name}/{sample_name}.discarded']}
	options = {
		'cores': 16,
		'memory': '60g',
		'walltime': '06:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/adapterremoval/{sample_name} ] || mkdir -p {output_directory}/adapterremoval/{sample_name}
	
	AdapterRemoval \
		--threads {options['cores']} \
		--file1 {' '.join(read_files)} \
		--adapter1 {adapter1} \
		--minquality {min_qulaity} \
		--minlength {min_length} \
		--basename {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog \
		--trimns \
		--trimqualities \
		--collapse
	
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.settings {outputs['misc'][0]}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.truncated {outputs['truncated']}
	mv {output_directory}/adapterremoval/{sample_name}/{sample_name}.prog.discarded {outputs['misc'][1]}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def alignment_pairedend(read1_file: str, read2_file: str, reference_genome_file: str, sample_name: str, output_directory: str):
	"""
	Template: Align and sort paired-end data to reference genome using :script:`bwa mem`.
	
	Template I/O::
	
		inputs = {'read1': read1_file, 'read2': read2_file, 'reference': reference_genome_file}
		outputs = {'alignment': *.paired_alignment.bam}
	
	:param str read1_file:
		Path to file containing mate 1 reads.
	:param str read2_file:
		Path to file containing mate 2 reads.
	:param str reference_genome_file:
		Path to indexed reference genome file.
	:param str sample_name:
		Name of sample.
	:param str output_directory:
		Path to output directory.
	"""
	inputs = {'read1': read1_file,
		   	  'read2': read2_file,
			  'reference': reference_genome_file}
	outputs = {'alignment': f'{output_directory}/alignment/{sample_name}/{sample_name}.paired_alignment.bam'}
	options = {
		'cores': 18,
		'memory': '60g',
		'walltime': '48:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/alignment/{sample_name}/tmp ] || mkdir -p {output_directory}/alignment/{sample_name}/tmp
	
	bwa mem \
		-t {options['cores']} \
		-R "@RG\\tID:{sample_name}\\tSM:{sample_name}" \
		{reference_genome_file} \
		{read1_file} \
		{read2_file} \
	| samtools sort \
		--threads {options['cores'] - 1} \
		-n \
		--output-fmt BAM \
		-T {output_directory}/alignment/{sample_name}/tmp \
		-o {output_directory}/alignment/{sample_name}/{sample_name}.paired_alignment.prog.bam \
		-
	
	mv {output_directory}/alignment/{sample_name}/{sample_name}.paired_alignment.prog.bam {outputs['alignment']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def alignment_collapsed(collapsed_read_files: list, reference_genome_file: str, sample_name: str, output_directory: str):
	"""
	Template: Align and sort collapsed paired-end data to reference genome using :script:`bwa mem`.
	
	Template I/O::
	
		inputs = {'collapse': collapsed_read_files, 'reference': reference_genome_file}
		outputs = {'alignment': *collapsed_alignment.bam}
	
	:param list collapsed_read_files:
		List of collapsed paired-end read files.
	:param str reference_genome_file:
		Path to indexed reference genome file.
	:param str sample_name:
		Name of sample.
	:param str output_directory:
		Path to output directory.
	"""
	inputs = {'collapse': collapsed_read_files,
		   	  'reference': reference_genome_file}
	outputs = {'alignment': f'{output_directory}/alignment/{sample_name}/{sample_name}.collapsed_alignment.bam'}
	options = {
		'cores': 18,
		'memory': '60g',
		'walltime': '48:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/alignment/{sample_name}/tmp ] || mkdir -p {output_directory}/alignment/{sample_name}/tmp
	
	bwa mem \
		-t {options['cores']} \
		-R "@RG\\tID:{sample_name}\\tSM:{sample_name}" \
		{reference_genome_file} \
		<(cat {' '.join(collapsed_read_files)}) \
	| samtools sort \
		--threads {options['cores'] - 1} \
		-n \
		--output-fmt BAM \
		-T {output_directory}/alignment/{sample_name}/tmp \
		-o {output_directory}/alignment/{sample_name}/{sample_name}.collapsed_alignment.prog.bam \
		-
	
	mv {output_directory}/alignment/{sample_name}/{sample_name}.collapsed_alignment.prog.bam {outputs['alignment']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def merge_alignments(alignment_files: list, sample_name: str, output_directory: str):
	"""
	Template: Merge alignment files using :script:`samtools merge`.
	
	Template I/O::
	
		inputs = {'alignments': alignment_files}
		outputs = {'merged': *.merged_alignment.bam}
	
	:param list alignment_files:
		List of paths to alignment files.
	:param str sample_name:
		Name of sample.
	:param str output_directory:
		Path to output directory
	"""
	inputs = {'alignments': alignment_files}
	outputs = {'merged': f'{output_directory}/alignment/{sample_name}/{sample_name}.merged_alignment.bam'}
	options = {
		'cores': 18,
		'memory': '30g',
		'walltime': '12:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/alignment/{sample_name} ] || mkdir -p {output_directory}/alignment/{sample_name}
	
	samtools merge \
		--threads {options['cores'] - 1} \
		-c \
		-p \
		-n \
		-o {output_directory}/alignment/{sample_name}/{sample_name}.merged_alignment.prog.bam \
		{' '.join(alignment_files)}
	
	mv {output_directory}/alignment/{sample_name}/{sample_name}.merged_alignment.prog.bam {outputs['merged']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def alignment_singleend(read_file: str, reference_genome_file: str, sample_name: str, output_directory: str, seed_length: int = 16500, missing_fraction: float = 0.1, max_gaps: int = 2):
	"""
	Template: Align and sort single-end data to reference genome using :script:`bwa aln`.
	
	Template I/O::
	
		inputs = {'read': read_file, 'reference': reference_genome_file}
		outputs = {'alignment': *.single_alignment.bam}
	
	:param str read_file:
		Path to single-end read file.
	:param str reference_genome_file:
		Path to indexed reference genome file.
	:param str sample_name:
		Name of sample.
	:param str output_directory:
		Path to output directory.
	:param int seed_length:
		Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled.
	:param float missing_fraction:
		The fraction of missing alignments given 2% uniform base error rate.
	:param int max_gaps:
		Maximum number of gap opens.
	"""
	inputs = {'read': read_file,
			  'reference': reference_genome_file}
	outputs = {'alignment': f'{output_directory}/alignment/{sample_name}/{sample_name}.single_alignment.bam'}
	options = {
		'cores': 18,
		'memory': '60g',
		'walltime': '48:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/alignment/{sample_name}/tmp ] || mkdir -p {output_directory}/alignment/{sample_name}/tmp
	
	bwa aln \
		-t {options['cores']} \
		-l {seed_length} \
		-n {missing_fraction} \
		-o {max_gaps} \
		{reference_genome_file} \
		{read_file} \
	| bwa samse \
		-r "@RG\\tID:{sample_name}\\tSM:{sample_name}" \
		{reference_genome_file} \
		- \
		{read_file} \
	| samtools sort \
		--threads {options['cores'] - 1} \
		-n \
		--output-fmt BAM \
		-T {output_directory}/alignment/{sample_name}/tmp \
		-o {output_directory}/alignment/{sample_name}/{sample_name}.single_alignment.prog.bam \
		- 
	
	mv {output_directory}/alignment/{sample_name}/{sample_name}.single_alignment.prog.bam {outputs['alignment']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def mark_duplicates_samtools(alignment_file: str, sample_name: str, output_directory: str):
	"""
	Template: Mark duplicate alignments using :script:`samtools markdup`.
	
	Template I/O::
	
		inputs = {'alignment': alignment_file}
		outputs = {'markdup': *.markdup.bam, 'bai': *.markdup.bam.bai, 'stats': *.markdup.bam.markdupstats}
	
	:param str alignment_file:
		Path to input alignment file.
	:param str sample_name:
		Name of sample.
	:param str output_directory:
		Path to output directory.
	"""
	inputs = {'alignment': alignment_file}
	outputs = {'markdup': f'{output_directory}/alignment/{sample_name}/{sample_name}.markdup.bam',
			   'bai': f'{output_directory}/alignment/{sample_name}/{sample_name}.markdup.bam.bai',
			   'stats': f'{output_directory}/alignment/{sample_name}/{sample_name}.markdup.bam.markdupstats'}
	options = {
		'cores': 18,
		'memory': '60g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/alignment/{sample_name}/tmp ] || mkdir -p {output_directory}/alignment/{sample_name}/tmp

	samtools fixmate \
		--threads {options['cores'] - 1} \
		-m \
		--output-fmt BAM \
		{alignment_file} \
		- \
	| samtools sort \
		--threads {options['cores'] - 1} \
		--output-fmt BAM \
		-T {output_directory}/alignment/{sample_name}/tmp \
		- \
	| samtools markdup \
		--threads {options['cores'] - 1} \
		--output-fmt BAM \
		-T {output_directory}/alignment/{sample_name}/tmp \
		-s \
		-f {output_directory}/alignment/{sample_name}/{sample_name}.markdup.bam.markdupstats \
		- \
		{output_directory}/alignment/{sample_name}/{sample_name}.markdup.bam
	
	samtools index \
		--threads {options['cores'] - 1} \
		-b \
		{output_directory}/alignment/{sample_name}/{sample_name}.markdup.bam \
		{output_directory}/alignment/{sample_name}/{sample_name}.markdup.prog.bam.bai

	mv {output_directory}/alignment/{sample_name}/{sample_name}.markdup.prog.bam.bai {outputs['bai']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def extract_unmapped_reads(alignment_file: str, sample_name: str, output_directory: str):
	"""
	Template: Extract unmapped reads from alignment based on bit-flag value 4.
	
	Template I/O::
	
		inputs = {'alignment': alignment_file}
		outputs = {'unmapped': *.unmapped.bam, 'bai': *.unmapped.bam.bai}
	
	:param str alignment_file:
		Path to input alignment file.
	:param str sample_name:
		Name of sample.
	:param str output_directory:
		Path to output directory.
	"""
	inputs = {'alignment': alignment_file}
	outputs = {'unmapped': f'{output_directory}/{sample_name}/{sample_name}.unmapped.bam',
			   'bai': f'{output_directory}/{sample_name}/{sample_name}.unmapped.bam.bai'}
	protect = [outputs['unmapped'], outputs['bai']]
	options = {
		'cores': 18,
		'memory': '30g',
		'walltime': '12:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/{sample_name} ] || mkdir -p {output_directory}/{sample_name}
	
	samtools view \
		--threads {options['cores'] - 1} \
		--require-flags 4 \
		--bam \
		--output {output_directory}/{sample_name}/{sample_name}.unmapped.bam \
		{alignment_file}

	samtools index \
		--threads {options['cores'] - 1} \
		--bai \
		--output {output_directory}/{sample_name}/{sample_name}.unmapped.prog.bam.bai \
		{output_directory}/{sample_name}/{sample_name}.unmapped.bam
		
	mv {output_directory}/{sample_name}/{sample_name}.unmapped.prog.bam.bai {outputs['bai']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def samtools_stats(alignment_file: str, output_directory: str):
	"""
	Template: Create various basic statistics on alignment using :script:`samtools`.
	
	Template I/O::
	
		inputs = {'alignment': alignment_file}
		outputs = {'stats': [*.idxstats, *.flagstat, *.coverage, *.stats]}
	
	:param str alignment_file:
		Path to input alignment file.
	:param str output_directory:
		Path to output directory.
	"""
	inputs = {'alignment': alignment_file}
	outputs = {'stats': [f'{output_directory}/{os.path.basename(alignment_file)}.idxstats',
						 f'{output_directory}/{os.path.basename(alignment_file)}.flagstat',
						 f'{output_directory}/{os.path.basename(alignment_file)}.coverage',
						 f'{output_directory}/{os.path.basename(alignment_file)}.stats']}
	protect = outputs['stats']
	options = {
		'cores': 18,
		'memory': '40g',
		'walltime': '12:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory} ] || mkdir -p {output_directory}
	
	samtools idxstats \
		--threads {options['cores'] - 1} \
		{alignment_file} \
		> {output_directory}/{os.path.basename(alignment_file)}.prog.idxstats
	
	samtools flagstat \
		--threads {options['cores'] - 1} \
		{alignment_file} \
		> {output_directory}/{os.path.basename(alignment_file)}.prog.flagstat

	samtools coverage \
		-o {output_directory}/{os.path.basename(alignment_file)}.prog.coverage \
		{alignment_file}
	
	samtools stats \
		--threads {options['cores'] - 1} \
		--coverage 1,1000,1 \
		{alignment_file} \
		> {output_directory}/{os.path.basename(alignment_file)}.prog.stats

	mv {output_directory}/{os.path.basename(alignment_file)}.prog.idxstats {outputs['stats'][0]}
	mv {output_directory}/{os.path.basename(alignment_file)}.prog.flagstat {outputs['stats'][1]}
	mv {output_directory}/{os.path.basename(alignment_file)}.prog.coverage {outputs['stats'][2]}
	mv {output_directory}/{os.path.basename(alignment_file)}.prog.stats {outputs['stats'][3]}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def samtools_filter(alignment_file: str, sample_name: str, output_directory: str, flags_excluded: int | None = 3844, flags_required: int | None = None, min_mq: int = 20):
	"""
	Template: Filter alignment file using :script:`samtools`.
	
	Template I/O::
	
		inputs = {'alignment': alignment_file}
		outputs = {'filtered': *.filtered.bam, 'bai': *.filtered.bam.bai}
	
	:param str alignment_file:
		PAth to input alignment file.
	:param str sample_name:
		Name of sample.
	:param str output_directory:
		Path to output directory.
	:param int | None flags_excluded:
		Do not output alignments with any of the speciefied bits set in the FLAG field.
	:param int | None flags_required:
		Only output alignments with all of the specified bits set in the FLAG field.
	:param int min_mq:
		Minimum MAPQ of alignments.
	"""
	flags = []
	if flags_excluded:
		flags.append(f'--exclude-flags {flags_excluded}')
	if flags_required:
		flags.append(f'--require-flags {flags_required}')
	inputs = {'alignment': alignment_file}
	outputs = {'filtered': f'{output_directory}/{sample_name}/{sample_name}.filtered.bam',
			   'bai': f'{output_directory}/{sample_name}/{sample_name}.filtered.bam.bai'}
	protect = [outputs['filtered'], outputs['bai']]
	options = {
		'cores': 18,
		'memory': '60g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/{sample_name} ] || mkdir -p {output_directory}/{sample_name}
	
	samtools view \
		--threads {options['cores'] - 1} \
		--bam \
		--min-MQ {min_mq} \
		--output {output_directory}/{sample_name}/{sample_name}.filtered.bam \
		{' '.join(flags)} \
		{alignment_file}

	samtools index \
		--threads {options['cores'] - 1} \
		--bai \
		--output {output_directory}/{sample_name}/{sample_name}.filtered.prog.bam.bai \
		{output_directory}/{sample_name}/{sample_name}.filtered.bam
	
	mv {output_directory}/{sample_name}/{sample_name}.filtered.prog.bam.bai {outputs['bai']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def qc_qualimap(alignment_file: str, output_directory: str):
	"""
	Template: Run :script:`qualimap` on alignment file.
	
	Template I/O::
	
		inputs = {'alignment': alignment_file}
		outputs = {'pdf': report.pdf, 'html': qualimapReport.html, 'raw':genome_results.txt}
	
	:param str alignment_file:
		PAth to input alignment file.
	:param str output_directory:
		PAth to output directory.
	"""
	inputs = {'alignment': alignment_file}
	outputs = {'pdf': f'{output_directory}/report.pdf',
			   'html': f'{output_directory}/qualimapReport.html',
			   'raw': f'{output_directory}/genome_results.txt'}
	protect = [outputs['pdf'], outputs['html'], outputs['raw']]
	options = {
		'cores': 32,
		'memory': '300g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory} ] || mkdir -p {output_directory}
	
	export _JAVA_OPTIONS="-Djava.awt.headless=true -Xmx{options['memory']}"
	
	qualimap bamqc \
		-nt {options['cores']} \
		-bam {alignment_file} \
		-outdir {output_directory} \
		-outformat PDF:HTML \
		-outfile report.prog.pdf \
		--java-mem-size={options['memory']}

	mv {output_directory}/report.prog.pdf {outputs['pdf']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def qualimap_multi(dataset: list, output_directory: str):
	"""
	Template: Compares QC report from multiple runs of :script:`qualimap`.
	
	Template I/O::
	
		inputs = {'raw': [*]}
		outputs = {'pdf': report.pdf, 'html': multisampleBAMQcReport.html}
	
	:param list dataset:
		List of paths to bamqc folders.
	:param str output_directory:
		Path to output directory.
	"""
	dataset_tabular = '\n'.join(['\t'.join(i) for i in dataset])
	inputs = {'raw': [f'{i[1]}/genome_results.txt' for i in dataset]}
	outputs = {'pdf': f'{output_directory}/report.pdf',
			   'html': f'{output_directory}/multisampleBamQcReport.html'}
	protect = [outputs['pdf'], outputs['html']]
	options = {
		'cores': 32,
		'memory': '300g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory} ] || mkdir -p {output_directory}

	export _JAVA_OPTIONS="-Djava.awt.headless=true -Xmx{options['memory']}"

	qualimap multi-bamqc \
		-d <(echo -e "{dataset_tabular}") \
		-outdir {output_directory} \
		-outfile report.prog.pdf \
		-outformat PDF:HTML
	
	mv {output_directory}/report.prog.pdf {outputs['pdf']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)