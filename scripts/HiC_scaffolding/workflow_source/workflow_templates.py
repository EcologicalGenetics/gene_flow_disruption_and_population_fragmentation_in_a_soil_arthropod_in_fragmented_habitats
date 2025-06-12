#!/bin/env python3
from gwf import AnonymousTarget
import os, sys, glob

########################## Functions ##########################

def species_abbreviation(species_name: str) -> str:
	"""Creates species abbreviation from species name.

	:param str species_name:
		Species name written as *genus* *species*"""
	genus, species = species_name.replace(' ', '_').split('_')
	genus = genus[0].upper() + genus[1:3]
	species = species[0].upper() + species[1:3]
	return genus + species

#--------------------------------------------------------------------------
########################## Hi-C scaffolding YaHS ##########################
#--------------------------------------------------------------------------

def index_reference(reference_genome_file: str, output_directory: str):
	"""
	Template: Index reference genome using both :script:`bwa index` and :script:`samtools faidx`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'reference': reference_genome_file}
	outputs = {'bwa': [f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.amb',
					   f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.ann',
					   f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.pac',
					   f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.bwt',
					   f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.sa'],
			   'fai': f'{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.fai'}
	protect = [outputs['bwa'][0], outputs['bwa'][1], outputs['bwa'][2], outputs['bwa'][3], outputs['bwa'][4], outputs['fai']]
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '02:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/HiC_scaffolding/YaHS/reference ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS/reference
	ln -s {reference_genome_file} {output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}

	bwa index \
		-p {output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)} \
		{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}
	
	samtools faidx \
		-o {output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.prog.fai \
		{output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}
	
	mv {output_directory}/HiC_scaffolding/YaHS/reference/{os.path.basename(reference_genome_file)}.prog.fai {outputs['fai']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def hic_alignment_to_draft_assembly(hic_seqeuence_files: list, draft_genome_file: str, reference_indices: list, output_directory: str, species_name: str):
	"""
	Template: Align Hi-C data to draft genome using :script:`bwa mem`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'hic_reads': hic_seqeuence_files,
		   	  'genome': draft_genome_file,
			  'indices': reference_indices}
	outputs = {'bam': f'{output_directory}/HiC_scaffolding/YaHS/alignment/{species_abbreviation(species_name)}.HiC_to_draft.bam'}
	options = {
		'cores': 36,
		'memory': '100g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/HiC_scaffolding/YaHS/alignment/tmp ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS/alignment/tmp
	
	bwa mem \
        -t {options['cores']} \
        -R '@RG\\tID:{species_abbreviation(species_name)}.HiC\\tSM:HiC_to_draft' \
        -S \
        -P \
        -5 \
        -T 0 \
        {draft_genome_file} \
        {hic_seqeuence_files[0]} \
        {hic_seqeuence_files[1]} \
    | samtools sort \
        -@ {int(options['cores']) - 1} \
		-n \
        -O BAM \
        -T {output_directory}/HiC_scaffolding/YaHS/alignment/tmp \
        -o {output_directory}/HiC_scaffolding/YaHS/alignment/{species_abbreviation(species_name)}.HiC_to_draft.prog.bam
	
	mv {output_directory}/HiC_scaffolding/YaHS/alignment/{species_abbreviation(species_name)}.HiC_to_draft.prog.bam {outputs['bam']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mark_duplicates_picard(alignment_bam_file: str, output_directory: str):
	"""
	Template: Mark duplicates in :format:`BAM` alignment file using PICARD :script:`MarkDuplicates`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'bam': alignment_bam_file}
	outputs = {'markdup': f'{os.path.splitext(alignment_bam_file)[0]}.markdup.bam',
			   'metric': f'{os.path.splitext(alignment_bam_file)[0]}.markdup.txt'}
	protect = outputs['metric']
	options = {
		'cores': 1,
		'memory': '100g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/HiC_scaffolding/YaHS/alignment/tmp ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS/alignment/tmp
	
	export _JAVA_OPTIONS="-Xmx{options['memory']}"

	picard MarkDuplicates \
		--INPUT {alignment_bam_file} \
		--OUTPUT {os.path.splitext(alignment_bam_file)[0]}.markdup.prog.bam \
		--METRICS_FILE {os.path.splitext(alignment_bam_file)[0]}.markdup.prog.txt \
		--TMP_DIR {output_directory}/HiC_scaffolding/YaHS/alignment/tmp
	
	mv {os.path.splitext(alignment_bam_file)[0]}.markdup.prog.bam {outputs['markdup']}
	mv {os.path.splitext(alignment_bam_file)[0]}.markdup.prog.txt {outputs['metric']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def hic_scaffolding_yahs(draft_genome_file: str, hic_to_draft_bam_file: str, output_directory: str, species_name: str):
	"""
	Template: Genome scaffolding with Hi-C data using :script:`YaHS`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': draft_genome_file,
		   	  'bam': hic_to_draft_bam_file}
	outputs = {'final_agp': f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}_scaffolds_final.agp',
               'final_fa': f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}_scaffolds_final.fa',
               'bin': f'{output_directory}/HiC/YaHS/{species_abbreviation(species_name)}.bin'}
	protect = [outputs['final_agp'], outputs['final_fa'], outputs['bin']]
	options = {
		'cores': 1,
		'memory': '100g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/HiC_scaffolding/YaHS ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS
	
	yahs \
        -r 500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000 \
        -o {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)} \
        {draft_genome_file} \
        {hic_to_draft_bam_file}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def yahs_conversion_manual_curation(hic_bin_file: str, scaffolds_final_agp_file: str, draft_assembly_fai_index_file: str, output_directory: str, species_name: str):
    """
    Template: Convert Hi-C alignment file to format required by juicer_tools.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'bin': hic_bin_file,
              'scaffolds': scaffolds_final_agp_file,
              'index': draft_assembly_fai_index_file}
    outputs = {'jbat': [f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.txt',
                        f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.liftover.agp',
                        f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.assembly',
                        f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.assembly.agp',
                        f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.log']}
    protect = outputs['jbat']
    options = {
        'cores': 1,
        'memory': '30g',
        'walltime': '12:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory}/HiC_scaffolding/YaHS/tmp ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS/tmp
    [ -d {output_directory}/HiC_scaffolding/YaHS ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS

    juicer pre \
        -a \
        -o {output_directory}/HiC_sacffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog \
        {hic_bin_file} \
        {scaffolds_final_agp_file} \
        {draft_assembly_fai_index_file} \
        > {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.log 2>&1
    
    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.txt {outputs['jbat'][0]}
    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.liftover.agp {outputs['jbat'][1]}
    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.assembly {outputs['jbat'][2]}
    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.assembly.agp {outputs['jbat'][3]}
    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.log {outputs['jbat'][4]}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def chrom_sizes(scaffolds_final_fa_file: str, chrom_size_script: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/get_length.py'):
    """
    Template: Create file containing two columns. In column 1 the name of each sequence, in column 2 the length of each sequence.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'scaffolds': scaffolds_final_fa_file}
    outputs = {'sizes': f'{scaffolds_final_fa_file}.chrom_sizes'}
    protect = outputs['sizes']
    options = {
        'cores': 1,
        'memory': '16g',
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
    
    python {chrom_size_script} \
        {scaffolds_final_fa_file} \
        {scaffolds_final_fa_file}.prog.chrom_sizes
    
    mv {scaffolds_final_fa_file}.prog.chrom_sizes {outputs['sizes']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def contact_matrix_manual_curation(JBAT_text_file: str, JBAT_log_file: str, output_directory: str, species_name: str, juicer_tools: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/juicer_tools.2.20.00.jar'):
    """
    Template: Generate Hi-C contact matrix using :script:`juicer pre`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'jbat': JBAT_text_file,
              'sizes': JBAT_log_file}
    outputs = {'hic': f'{output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.hic'}
    protect = outputs['hic']
    options = {
        'cores': 32,
        'memory': '100g',
        'walltime': '48:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory}/HiC_scaffolding/YaHS/tmp ] || mkdir -p {output_directory}/HiC_scaffolding/YaHS/tmp

    java -jar -Xmx{options['memory']} {juicer_tools} pre \
        -t {output_directory}/HiC_scaffolding/YaHS/tmp \
        -j {options['cores']} \
        {JBAT_text_file} \
        {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.hic \
        <(cat {JBAT_log_file} | grep PRE_C_SIZE | awk '{{print $2" "$3}}')

    mv {output_directory}/HiC_scaffolding/YaHS/{species_abbreviation(species_name)}.JBAT.prog.hic {outputs['hic']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

#----------------------------------------------------------------------------
########################## Hi-C scaffolding juicer ##########################
#----------------------------------------------------------------------------

def setup_for_juicer(hic_read1_files: list, hic_read2_files: list, draft_genome_assembly_file: str, output_directory: str, split_size: int = 8000000):
	"""
	Template: Set up directory structure required by the :script:`juicer` pipeline.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'hic_read1': hic_read1_files,
		   	  'hic_read2': hic_read2_files,
		   	  'genome': draft_genome_assembly_file}
	outputs = {'reference': f'{output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}',
			   'fastq1': [f'{output_directory}/juicer/fastq/{i + 1}_R1.{"fastq.gz" if j.endswith(".gz") else ".fastq"}' for i, j in enumerate(hic_read1_files)],
			   'fastq2': [f'{output_directory}/juicer/fastq/{i + 1}_R2.{"fastq.gz" if j.endswith(".gz") else ".fastq"}' for i, j in enumerate(hic_read2_files)],
			   'bwa': [f'{output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}.amb',
					   f'{output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}.ann',
					   f'{output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}.pac',
					   f'{output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}.bwt',
					   f'{output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}.sa'],
			   'fai': f'{output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}.fai',
			   'sizes': f'{output_directory}/juicer/chrom.sizes'}
	protect = [outputs['reference'], outputs['bwa'][0], outputs['bwa'][1], outputs['bwa'][2], outputs['bwa'][3], outputs['bwa'][4], outputs['fai'], outputs['sizes']]
	options = {
		'cores': 1,
		'memory': '20g',
		'walltime': '36:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/juicer ] || mkdir -p {output_directory}/juicer
	[ -d {output_directory}/juicer/aligned ] && rm -rf {output_directory}/juicer/aligned
	[ -d {output_directory}/juicer/debug ] && rm -rf {output_directory}/juicer/debug
	[ -d {output_directory}/juicer/HIC_tmp ] && rm -rf {output_directory}/juicer/HIC_tmp
	[ -d {output_directory}/juicer/fastq ] && rm -rf {output_directory}/juicer/fastq
	mkdir -p {output_directory}/juicer/fastq
	[ -d {output_directory}/juicer/references ] && rm -rf {output_directory}/juicer/references
	mkdir -p {output_directory}/juicer/references
	[ -d {output_directory}/juicer/restriction_sites ] || mkdir -p {output_directory}/juicer/restriction_sites
	[ -d {output_directory}/juicer/splits ] && rm -rf {output_directory}/juicer/splits 
	mkdir -p {output_directory}/juicer/splits
	
	[ -e {outputs['reference']} ] || ln -s {draft_genome_assembly_file} {outputs['reference']}

	read1=({" ".join(hic_read1_files)})
	fastq1=({" ".join(outputs['fastq1'])})
	read2=({" ".join(hic_read2_files)})
	fastq2=({" ".join(outputs['fastq2'])})

	for (( i = 0; i < ${{#read1[@]}}; i++ )); do
		[ -e "${{fastq1[$i]}}" ] || ln -s "${{read1[$i]}}" "${{fastq1[$i]}}"
	done

	for (( i = 0; i < ${{#read2[@]}}; i++ )); do
		[ -e "${{fastq2[$i]}}" ] || ln -s "${{read2[$i]}}" "${{fastq2[$i]}}"
	done
	
	bwa index \
		-p {output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)} \
		{output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}
	
	samtools faidx \
		-o {output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}.prog.fai \
		{output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}
	
	mv {output_directory}/juicer/references/{os.path.basename(draft_genome_assembly_file)}.prog.fai {outputs['fai']}

	split \
		-a 3 \
		-l 8000000 \
		-d \
		--additional-suffix=_R1.fastq \
		--filter='gzip > $FILE.gz' \
		<(zcat {" ".join(hic_read1_files)}) \
		{output_directory}/juicer/splits/
	
	split \
		-a 3 \
		-l 8000000 \
		-d \
		--additional-suffix=_R2.fastq \
		--filter='gzip > $FILE.gz' \
		<(zcat {" ".join(hic_read2_files)}) \
		{output_directory}/juicer/splits/
	
	awk \
		'BEGIN{{OFS = "\\t"}}
		{{print $1, $2}}' \
		{outputs['fai']} \
		> {output_directory}/juicer/chrom.prog.sizes

	mv {output_directory}/juicer/chrom.prog.sizes {outputs['sizes']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def juicer_pipeline(draft_genome_assembly_file: str, chrom_sizes_file: str, output_directory: str, species_name: str, juicer: str = f'{os.path.dirname(os.path.realpath(__file__))}/software/juicer-1.6/scripts/juicer.sh'):
	"""
	Template: Starts the :script:`juicer` pipeline.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': draft_genome_assembly_file,
		   	  'sizes': chrom_sizes_file}
	outputs = {'nodups': f'{output_directory}/juicer/aligned/merged_nodups.txt'}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '01:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	bash {juicer} \
		-d {output_directory}/juicer \
		-D {os.path.dirname(os.path.dirname(juicer))} \
		-p {chrom_sizes_file} \
		-A EcoGenetics \
		-s none \
		-z {draft_genome_assembly_file} \
		-q short \
		-Q 12:00:00 \
		-l normal \
		-L 24:00:00 \
		-t 40 \
		> {species_abbreviation(species_name)}_juicer.log
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def assembly_3ddna(draft_assembly_file: str, merged_nodups_file: str, output_directory: str, input_size: int = 15000, edit_rounds: int = 3):
	"""
	Template: 3D de novo assembly
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': draft_assembly_file,
		   	  'nodups': merged_nodups_file}
	outputs = {'final_fasta': f'{output_directory}/juicer/3ddna_in{input_size}_r{edit_rounds}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.final.fasta',
			   'final_asm': f'{output_directory}/juicer/3ddna_in{input_size}_r{edit_rounds}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.final.asm',
			   'final_assembly': f'{output_directory}/juicer/3ddna_in{input_size}_r{edit_rounds}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.final.assembly',
			   'final_cprops': f'{output_directory}/juicer/3ddna_in{input_size}_r{edit_rounds}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.final.cprops',
			   'final_hic': f'{output_directory}/juicer/3ddna_in{input_size}_r{edit_rounds}/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.final.hic'}
	options = {
		'cores': 40,
		'memory': '80g',
		'walltime': '48:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {output_directory}/juicer/3ddna_in{input_size}_r{edit_rounds}/tmp ] || mkdir -p {output_directory}/juicer/3ddna_in{input_size}_r{edit_rounds}/tmp
	
	cd {output_directory}/juicer/3ddna_in{input_size}_r{edit_rounds}

	export _JAVA_OPTIONS="-Djava.io.tmpdir={output_directory}/juicer/3ddna_in{input_size}_r{edit_rounds}/tmp -Xmx{options['memory']}"

	3d-dna \
		--rounds {edit_rounds} \
		--input {input_size} \
		{draft_assembly_file} \
		{merged_nodups_file}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def finalize_3ddna(reviewed_assembly_file: str, draft_assembly_file: str, merged_nodups_file: str, number_of_chromosomes: int, final_hic_file: str, post_review: str = f'{os.path.dirname(os.path.dirname(os.path.realpath(sys.executable)))}/share/3d-dna/run-asm-pipeline-post-review.sh'):
	"""
	Template: Finalizes draft assembly using reviewed :format:`assembly` file from :script:`JuiceBox`.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'asm': reviewed_assembly_file,
		   	  'assembly': draft_assembly_file,
			  'nodups': merged_nodups_file,
			  'hic': final_hic_file}
	outputs = {'fasta': f'{os.path.dirname(final_hic_file)}/finalize/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}_HiC.fasta',
			   'asm': f'{os.path.dirname(final_hic_file)}/finalize/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}_HiC.assembly',
			   'final_fasta': f'{os.path.dirname(final_hic_file)}/finalize/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.nodebris.fasta'}
	options = {
		'cores': 40,
		'memory': '80g',
		'walltime': '12:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	[ -d {os.path.dirname(final_hic_file)}/finalize ] || mkdir -p {os.path.dirname(final_hic_file)}/finalize

	cd {os.path.dirname(final_hic_file)}/finalize

	bash {post_review} \
		--sort-output \
		-s finalize \
		-r {reviewed_assembly_file} \
		{draft_assembly_file} \
		{merged_nodups_file}

	seqkit range \
		-j {options['cores']} \
		-r 1:{number_of_chromosomes} \
		-o {os.path.dirname(final_hic_file)}/finalize/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.nodebris.prog.fasta \
		{os.path.dirname(final_hic_file)}/finalize/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}_HiC.fasta
	
	mv {os.path.dirname(final_hic_file)}/finalize/{os.path.splitext(os.path.basename(draft_assembly_file))[0]}.nodebris.prog.fasta {outputs['final_fasta']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def busco_genome(genome_assembly_file: str, busco_dataset: str, busco_download_path: str = '/faststorage/project/EcoGenetics/databases/BUSCO'):
	"""
	Template: Runs BUSCO analysis on genome assembly.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genome': genome_assembly_file}
	outputs = {'stattxt': f'{os.path.dirname(genome_assembly_file)}/busco_{os.path.basename(genome_assembly_file)}/short_summary.specific.{busco_dataset}.busco_{os.path.basename(genome_assembly_file)}.txt',
			   'statjson': f'{os.path.dirname(genome_assembly_file)}/busco_{os.path.basename(genome_assembly_file)}/short_summary.specific.{busco_dataset}.busco_{os.path.basename(genome_assembly_file)}.json'}
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
	
	cd {os.path.dirname(genome_assembly_file)}

	busco \
		--cpu {options['cores']} \
		--force \
		--in {genome_assembly_file} \
		--mode genome \
		--out busco_{os.path.basename(genome_assembly_file)} \
		--out_path {os.path.dirname(genome_assembly_file)} \
		--download_path {busco_download_path} \
		--lineage {busco_download_path}/lineages/{busco_dataset} \
		--tar \
		--offline

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def bundle_sequences(assembly_fasta_file: str, n_chromosomes_to_keep: int, n_insertion_size: int = 1000):
	"""
	Template: template_description
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'assembly': assembly_fasta_file}
	outputs = {'bundled': f'{os.path.splitext(assembly_fasta_file)[0]}.bundled{os.path.splitext(assembly_fasta_file)[1]}',
			   'bed': f'{os.path.splitext(assembly_fasta_file)[0]}.bundled.bed'}
	options = {
		'cores': 1,
		'memory': '10g',
		'walltime': '02:00:00'
	}
	spec = f"""
	# Sources environment
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate assembly
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	awk \
		-v chromosomes_to_keep={n_chromosomes_to_keep} \
		-v spacing={n_insertion_size} \
		'BEGIN{{
			RS = ">"
			ORS = ""
			FS = "\\n"
			OFS = "\\n"
			ndebris = 0
			ngaps = 0
			start = 0
			end = 0
			if (! spacing)
			{{
				spacing = 1000
			}}
		}}
		{{
			if (NR > 1 && NR <= chromosomes_to_keep + 1)
			{{
				print ">" $1 "\\n"
				for (i = 2; i <= NF; i++)
				{{
					print $i
				}}
				print "\\n"
			}}
			if (NR > chromosomes_to_keep + 1)
			{{
				seqid = $1
				if (ndebris == 0)
				{{
					print ">bundled_sequences\\n"
				}}
				ndebris = 1
				if (NR > chromosomes_to_keep + 2)
				{{
					for (i = 1; i <= spacing; i++)
					{{
						print "N"
					}}
					ngaps += 1
					end += spacing
					print "bundled_sequences\\t" start "\\t" end "\\tN_gap_" ngaps "\\n" > "{os.path.splitext(assembly_fasta_file)[0]}.bundled.prog.bed"
					start = end
				}}
				for (i = 2; i <= NF; i++)
				{{
					end += length($i)
					print $i
				}}
				print "bundled_sequences\\t" start "\\t" end "\\t" seqid "\\n" > "{os.path.splitext(assembly_fasta_file)[0]}.bundled.prog.bed"
				start = end 
			}}
		}}
		END{{
			print "\\n"
		}}' \
		{assembly_fasta_file} \
	| fold \
		> {os.path.splitext(assembly_fasta_file)[0]}.bundled.prog{os.path.splitext(assembly_fasta_file)[1]}
	
	mv {os.path.splitext(assembly_fasta_file)[0]}.bundled.prog{os.path.splitext(assembly_fasta_file)[1]} {outputs['bundled']}
	mv {os.path.splitext(assembly_fasta_file)[0]}.bundled.prog.bed {outputs['bed']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)