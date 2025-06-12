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

########################## RepeatMasking ##########################

def build_repeatmodeller_database(genome_assembly_file: str, working_directory: str, species_name: str):
	"""
	Template: Build ReapeatModeler database from :format:`FASTA` file using :script:`BuildDatabase` from RepeatModeler.
	
	Template I/O::
	
		inputs = {'assembly': genome_assembly_file}
		outputs = {'db_files': ['*.nhr', '*.nin', '*.njs', '*.nnd', '*.nog', '*.nsq', '*.translation']}
	
	:param str genome_assembly_file:
		Genome sequence :format:`FASTA` file to create database from.
	:param str output_directory:
		Output directory for resulting files.
	:param str species_name:
		Name of species being worked on.
	"""
	working_dir = working_directory
	inputs = {'assembly': genome_assembly_file}
	outputs = {'db_files': [f'{working_directory}/{species_name.replace(" ", "_")}.nhr',
						 	f'{working_directory}/{species_name.replace(" ", "_")}.nin',
						 	f'{working_directory}/{species_name.replace(" ", "_")}.njs',
						 	f'{working_directory}/{species_name.replace(" ", "_")}.nnd',
						 	f'{working_directory}/{species_name.replace(" ", "_")}.nni',
						 	f'{working_directory}/{species_name.replace(" ", "_")}.nog',
						 	f'{working_directory}/{species_name.replace(" ", "_")}.nsq',
						 	f'{working_directory}/{species_name.replace(" ", "_")}.translation']}
	options = {
		'cores': 2,
		'memory': '30g',
		'walltime': '01:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate repeatmasking
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {working_directory} ] || mkdir -p {working_directory}

	BuildDatabase \
		-name {species_name.replace(' ', '_')} \
		-engine rmblast \
		{genome_assembly_file}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, options=options, spec=spec)

def repeatmodeler(database: list, working_directory: str, species_name: str):
	"""
	Template: Models repetitive DNA using RepeatModeler. Furthermore adds species code to results and creates separate files of classified and unclassified elements.
	
	Template I/O::
	
		inputs = {'database': database}
		outputs = {'model_result': ['*-families.fa', '*-families.stk', '*-rmod.log'],
			   'all': '*-families.prefix.fa',
			   'unknown': '*-families.prefix.unknown.fa',
			   'known': '*-families.prefix.known.fa'}
	
	:param list database:
		List of files produced by **build_repeatmodeler_database**.
	:param str working_directory:
		Output directory for files. The script uses it as its working directory thus the name.
	:param str species_name:
		Name of species being worked on.
	"""
	working_dir = working_directory
	inputs = {'database': database}
	outputs = {'model_result': [f'{working_directory}/RM_DB_{species_name.replace(" ", "_")}/{species_name.replace(" ", "_")}-families.fa',
						 		f'{working_directory}/RM_DB_{species_name.replace(" ", "_")}/{species_name.replace(" ", "_")}-families.stk',
						 		f'{working_directory}/RM_DB_{species_name.replace(" ", "_")}/{species_name.replace(" ", "_")}-rmod.log'],
			   'all': f'{working_directory}/{species_name.replace(" ", "_")}-families.prefix.fa',
			   'unknown': f'{working_directory}/{species_name.replace(" ", "_")}-families.prefix.unknown.fa',
			   'known': f'{working_directory}/{species_name.replace(" ", "_")}-families.prefix.known.fa'}
	options = {
		'cores': 32,
		'memory': '192g',
		'walltime': '72:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate repeatmasking
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	RepeatModeler \
		-database ./RM_DB_{species_name.replace(' ', '_')}/{species_name.replace(' ', '_')} \
		-LTRStruct \
		-pa {int(options['cores']/4)}
	
	seqkit \
		fx2tab \
		RM_DB_{species_name.replace(' ', '_')}/{species_name.replace(' ', '_')}-families.fa \
	| awk \
		-v species_abbreviation={species_abbreviation(species_name)} \
		'{{print species_abbreviation"_"$0}}' \
	| seqkit \
		tab2fx \
		> {species_name.replace(' ', '_')}-families.prefix.prog.fa

	seqkit \
		fx2tab \
		{species_name.replace(' ', '_')}-families.prefix.prog.fa \
	| awk \
		'{{if ($0 ~ /Unknown/) {{print $0}}}}' \
	| seqkit \
		tab2fx
		> {species_name.replace(' ', '_')}-families.prefix.unknown.prog.fa

	seqkit \
		fx2tab \
		{species_name.replace(' ', '_')}-families.prefix.prog.fa \
	| awk \
		'{{if ($0 !~ /Unknown/) {{print $0}}}}' \
	| seqkit \
		tab2fx
		> {species_name.replace(' ', '_')}-families.prefix.known.prog.fa

	mv {species_name.replace(' ', '_')}-families.prefix.prog.fa {outputs['all']}
	mv {species_name.replace(' ', '_')}-families.prefix.unknown.prog.fa {outputs['unknown']}
	mv {species_name.replace(' ', '_')}-families.prefix.known.prog.fa {outputs['known']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, options=options, spec=spec)

def repeatmasker(genome_assembly_file: str, library_file: str, working_directory: str):
	"""
	Template: Masks repetitive DNA using RepeatMasker and produces :format:`GFF3` annotation file.
	
	Template I/O::
	
		inputs = {'assembly': genome_assembly_file,
			  'lib': library_file}
		outputs = {'gff': '*repeats.'run_name'.gff',
			   'repmaskout': ['*.cat.gz', '*.masked', '*.ori.out', '*.out', '*.tbl']}
	
	:param str genome_assembly_file:
		Genome assembly file in which to mask repeats.
	:param str library_file:
		Library file containing repeat sequences.
	:param str output_directory:
		Output directory for resulting files.
	:param str run_name:
		Name to identify the run. Adviced to indicate library file used.
	"""
	working_dir = working_directory
	inputs = {'assembly': genome_assembly_file,
		   	  'lib': library_file}
	outputs = {'gff': f'{working_directory}/{os.path.basename(genome_assembly_file)}.repeats.{os.path.basename(working_directory)}.gff',
			   'repmaskout': [f'{working_directory}/{os.path.basename(genome_assembly_file)}.cat.gz',
						 	  f'{working_directory}/{os.path.basename(genome_assembly_file)}.masked',
						 	  f'{working_directory}/{os.path.basename(genome_assembly_file)}.out',
						 	  f'{working_directory}/{os.path.basename(genome_assembly_file)}.tbl']}
	options = {
		'cores': 32,
		'memory': '192g',
		'walltime': '08:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate repeatmasking
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {working_directory} ] || mkdir -p {working_directory}

	RepeatMasker \
		-e rmblast \
		-pa {int(options['cores']/4)} \
		-dir {working_directory} \
		-xsmall \
		-lib {library_file} \
		{genome_assembly_file}
	
	awk \
		-F " " \
		'BEGIN{{OFS = "\\t"; print "##gff-version 3";}}
		{{if (NR > 3)
			{{if ($9 == "C")
				{{strand = "-";}}
			else
				{{strand = "+";}}
			if ($12 ~ /\(/)
				{{start = $14;}}
			else
				{{start = $12;}}
			print $5, "RepeatMasker", "repeat_region", $6, $7, ".", strand, ".", "ID="$15";Name="$10";Class="$11";Family="$11";Target="$10" "start" "$13;
			}}
		}}' \
		{working_directory}/{os.path.basename(genome_assembly_file)}.out \
		> {working_directory}/{os.path.basename(genome_assembly_file)}.repeats.{os.path.basename(working_directory)}.prog.gff

	mv {working_directory}/{os.path.basename(genome_assembly_file)}.repeats.{os.path.basename(working_directory)}.prog.gff {outputs['gff']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, options=options, spec=spec)

def combine_repeatmasker_runs(repeatmasker_run1: list, repeatmasker_run2: list, library_file: str, working_directory: str, species_name: str):
	"""
	Template: Combines the output of two different RepeatMasker runs to one complete set of files.
	
	Template I/O::
	
		inputs = {'run1': ['*.cat.gz', '*.masked', '*.ori.out', '*.fasta.out', '*.tbl', '*.gff'],
			  'run2': ['*.cat.gz', '*.masked', '*.ori.out', '*.fasta.out', '*.tbl', '*.gff'],
			  'lib': library_file}
		outputs = {'cat': '*.fasta.fullmask.cat.gz',
			   'out': '*.fasta.fullmask.out',
			   'gff': '*.repeats.gff',
			   'other': ['*.fasta.fullmask.ori.out',
				     '*.fasta.fullmask.tbl']}
	
	:param str repeatmasker_run:
		List of files from first RepeatMasker run.
	:param str repeatmasker_run2:
		List of files from second RepeatMasker run.
	:param str library_file:
		Library file containing repeat sequences.
	:param str output_directory:
		Output directory for resulting files.
	:param str species_name:
		Name of species being worked on.
	"""
	working_dir = working_directory
	inputs = {'run1': repeatmasker_run1,
			  'run2': repeatmasker_run2,
			  'lib': library_file}
	outputs = {'cat': f'{working_directory}/{species_abbreviation(species_name)}.fasta.fullmask.cat.gz',
			   'out': f'{working_directory}/{species_abbreviation(species_name)}.fasta.fullmask.out',
			   'gff': f'{working_directory}/{species_abbreviation(species_name)}.repeats.gff',
			   'other': [f'{working_directory}/{species_abbreviation(species_name)}.fasta.fullmask.ori.out',
						 f'{working_directory}/{species_abbreviation(species_name)}.fasta.fullmask.tbl']}
	options = {
		'cores': 2,
		'memory': '32g',
		'walltime': '02:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate repeatmasking
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {working_directory} ] || mkdir -p {working_directory}

	cat \
		{inputs['run1'][0]} \
		{inputs['run2'][0]} \
		> {species_abbreviation(species_name)}.fasta.fullmask.prog.cat.gz

	cat \
		{inputs['run1'][2]} \
		<(tail -n +4 {inputs['run2'][2]}) \
		> {species_abbreviation(species_name)}.fasta.fullmask.prog.out
	
	mv {species_abbreviation(species_name)}.fasta.fullmask.prog.cat.gz {outputs['cat']}
	mv {species_abbreviation(species_name)}.fasta.fullmask.prog.out {outputs['out']}
		
	ProcessRepeats \
		-lib {library_file} \
		{species_abbreviation(species_name)}.fasta.fullmask.cat.gz
	
	awk \
		-F " " \
		'BEGIN{{OFS = "\\t"; print "##gff-version 3";}}
		{{if (NR > 3)
			{{if ($9 == "C")
				{{strand = "-";}}
			else
				{{strand = "+";}}
			if ($12 ~ /\(/)
				{{start = $14;}}
			else
				{{start = $12;}}
			print $5, "RepeatMasker", "repeat_region", $6, $7, ".", strand, ".", "ID="$15";Name="$10";Class="$11";Family="$11";Target="$10" "start" "$13;
			}}
		}}' \
		{species_abbreviation(species_name)}.fasta.fullmask.out \
		> {species_abbreviation(species_name)}.repeats.prog.gff

	mv {species_abbreviation(species_name)}.repeats.prog.gff {outputs['gff']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, options=options, spec=spec)

def mask_assembly(genome_assembly_file: str, annotation_file: str, working_directory: str, species_name: str):
	"""
	Template: Creates softmasked version of genome assembly file using :script:`bedtools`.
	
	Template I/O::
	
		inputs = {'assembly': genome_assembly_file,
			  'annotation': annotation_file}
		outputs = {'masked': '*.softmasked.fna',
			   'bed': '*repeats.bed'}
	
	:param str genome_assembly_file:
		Genome assembly file to be soft-masked.
	:param str annotation_file:
		Annotation file containing repeat regions to soft-mask.
	:param str output_directory:
		Output directory for resulting files.
	:param str species_name:
		Name of species being worked on.
	"""
	working_dir = working_directory
	inputs = {'assembly': genome_assembly_file,
		   	  'annotation': annotation_file}
	outputs = {'masked': f'{working_directory}/{species_abbreviation(species_name)}.softmasked.fna',
			   'bed': f'{working_directory}/{species_abbreviation(species_name)}.repeats.bed'}
	options = {
		'cores': 2,
		'memory': '32g',
		'walltime': '02:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate repeatmasking
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	bedtools maskfasta \
		-soft \
		-fi {genome_assembly_file} \
		-bed {annotation_file} \
		-fo {species_abbreviation(species_name)}.softmasked.prog.fna

	awk \
		'{{gsub(/\s$/, ""); print $0}}' \
		{species_abbreviation(species_name)}.softmasked.prog.fna \
		> {species_abbreviation(species_name)}.softmasked.nowhitespace.prog.fna

	awk \
		-F "\\t" \
		'BEGIN{{OFS = "\\t"}}
		{{if ($0 ~ /^[^#]/)
			{{print $1, ($4 - 1), $5}}
		}}' \
		{annotation_file} \
		> {species_abbreviation(species_name)}.repeats.prog.bed
	
	mv {species_abbreviation(species_name)}.softmasked.nowhitespace.prog.fna {outputs['masked']}
	rm {species_abbreviation(species_name)}.softmasked.prog.fna
	mv {species_abbreviation(species_name)}.repeats.prog.bed {outputs['bed']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(working_dir=working_dir, inputs=inputs, outputs=outputs, options=options, spec=spec)