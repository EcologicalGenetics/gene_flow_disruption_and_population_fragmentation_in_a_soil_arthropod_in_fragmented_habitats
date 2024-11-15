#!/bin/env python3
from gwf import AnonymousTarget
import os, glob
import subprocess as sp


########################## Functions ##########################

def species_abbreviation(species_name: str) -> str:
	"""Function: Creates species abbreviation from species name.

	:param str species_name:
		Species name written as *genus* *species*"""
	genus, species = species_name.replace(' ', '_').split('_') # first changes space to underscore, then splits at underscore. could have split at space?
	genus = genus[0].upper() + genus[1:3] # returns first character ([0]) as upper case? apparently python 1:3 does not include the last number. it means "up until"
	species = species[0].upper() + species[1:3]
	return genus + species

def awk_parse_pops_from_vcf(inputfile, outputfile):
    cmd = """awk '/^#CHROM/ {for(i=10; i<=NF; i++) print $i, "\t", i-9; exit}' """ + inputfile + " > " + outputfile
    os.system(cmd)

def awk_create_pairs(inputfile, outputfile):
    cmd = """awk -v OFS='\t' 'NR==FNR {id[NR] = $1; counter[NR] = $2; next} {for (i=1; i<FNR; i++) if (counter[i] <= $2) print id[i], $1, counter[i], $2}' """ + inputfile + " " + inputfile +" > " + outputfile	 	
    os.system(cmd)


def prepare_data(vcf_file: str, output_prefix: str, working_dir: str, pop_pairs_file: str):
	"""
	Template: Create fai indexed genome file for bedtools.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf_file': vcf_file}
	outputs = {'outfile_dat': f'{working_dir}/adata_prep/{output_prefix}_out2.tsv'}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '24:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate migration_fsc
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	# make wd
	echo Making wd in {working_dir}/adata_prep/
	mkdir -p {working_dir}
	mkdir -p {working_dir}/adata_prep/

	
	echo Prepping additional data for 2dSFS

	awk 'BEGIN{{ OFS="\t" }}
	{{if ($0 ~ /^#/)
			{{next}}
	else
			{{for(i=10;i<=NF;i++)
					if (substr($i,length($i)-2,3)=="0.0") {{sub($i, "0", $i)}} else {{sub($i,substr($i,length($i)-3,4), $i)}}
			}}; print
	}}' {inputs['vcf_file']} > {outputs['outfile_dat'].replace("out2", "out")}
	#out
	
	sed -i 's/://g' {outputs['outfile_dat'].replace("out2", "out")}

	awk 'BEGIN{{ OFS="\t" }}
			{{ for(i=10;i<=NF;i++)
					sub($i, $i*100, $i)
			; print $0 }}
	' {outputs['outfile_dat'].replace("out2", "out")} > {outputs['outfile_dat'].replace("out2", "out1")}
	# out1

	sed -i '/^#/d' {outputs['outfile_dat'].replace("out2", "out1")}
	#out1

	awk 'BEGIN {{ FS="\t"; OFS="\t" }} {{$1=$2=$3=$4=$5=$6=$7=$8=$9="";gsub(",+",",",$0)}}1' {outputs['outfile_dat'].replace("out2", "out1")} > {outputs['outfile_dat']} 
	sed -i 's/\t\\+/\t/g;s/^\t//' {outputs['outfile_dat']}
	sed -i '1d' {outputs['outfile_dat']} 
	echo "Data prep finished"

	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def create_input_dict_2dSFS(pair_file: str, exclude_list: list, include_list: list) -> list:
	"""Function: Creates input dictionary for running a map function of 2dSFS for each pop pair
	"""
	dict_list = []
	file = open(pair_file)
	for line in file:
		#print(line)
		infos = line.strip("\n").split("\t")
		#print(infos)
		pop1 = infos[0]
		#print(pop1)
		pop2 = infos[1]
		#print(pop2)
		if any(pop1 in s for s in exclude_list):
			#print(pop1 + " in exclude_list")
			continue
		if any(pop2 in s for s in exclude_list):
			#print(pop2 + " in exclude_list")
			continue
		conditions = [any(pop1 in s for s in include_list), any(pop2 in s for s in include_list)]
		if all(conditions):
			#print("yes")
			c1 = infos[2]
			c2 = infos[3]
			name = f'2dSFS_{pop1}_vs_{pop2}_{c1}_vs_{c2}'
			output_name = f'{name}/{name}_sfs2d_jointMAFpop1_0.obs'
			new_dict = {"pop1": pop1, "pop2": pop2, "c1": c1, "c2": c2, "name": name, "outname": output_name}
					#print(new_dict)
			dict_list.append(new_dict)
		#else:
		#	print("no")
	
	return dict_list


def create_run_name_fsc_pair(idx: str, target: AnonymousTarget) -> str:
	#pair_name=target.inputs['name_pops']
	return f'pair_2DSFS_map_target_{os.path.basename(target.outputs["outfile_path"]).replace("-","_")}'


def pair_2DSFS_map_target(pop1: str, pop2: str, c1: int, c2: int, name: str, out2_file: str, working_directory: str, outname: str):
	"""
	Template: Create netural bed file from genes and repeats bed for species.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {}
	##inputs = {'sfs_data': out2_file}
	outputs = {'outfile_path': f'{working_directory}/{outname}'} # OBS define outout
	options = {
		'cores': 1,
		'memory': '1g',
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate migration_fsc

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	mkdir -p {working_directory}/{name}
	##cp inputs['sfs_data'] /scratch/$SLURM_JOBID/out2
	## echo {out2_file}
	##../../../workflow_source/scripts/2dsfs.sh /scratch/$SLURM_JOBID/out2 {pop1} {pop2} {c1} {c2} {name} {working_directory} {outputs['outfile_path']}
	
	#mv {outputs['outfile_path'].replace(".obs", ".obs.tmp")} {outputs['outfile_path']}
	
	if [ -f {outputs['outfile_path']}]
	then
		echo Checked if outfile exist, and it does. File: {outputs['outfile_path']}
	else
		echo Outout file does not exist
	fi


	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def pair_2DSFS_map_target_noSingletons(pop1: str, pop2: str, c1: int, c2: int, name: str, out2_file: str, working_directory: str, outname: str, working_directory_alldat: str):
	"""
	Template: Create netural bed file from genes and repeats bed for species.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'sfs_to_copy': f'{working_directory_alldat}/{outname}'}
	##inputs = {'sfs_data': out2_file}
	outputs = {'sfs_newrun_outfile_path': f'{working_directory}/{outname.replace("_sfs2d_jointMA", "_NoSingl_sfs2d_jointMA")}'} # OBS define outout
	options = {
		'cores': 1,
		'memory': '1g',
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate migration_fsc

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	mkdir -p {working_directory}/{name}
	##cp inputs['sfs_data'] /scratch/$SLURM_JOBID/out2
	## echo {out2_file}
	##../../../workflow_source/scripts/2dsfs.sh /scratch/$SLURM_JOBID/out2 {pop1} {pop2} {c1} {c2} {name} {working_directory} {outputs['sfs_newrun_outfile_path']}
	
	#mv {outputs['sfs_newrun_outfile_path'].replace(".obs", ".obs.tmp")} {outputs['sfs_newrun_outfile_path']}
	cp {inputs['sfs_to_copy']} {outputs['sfs_newrun_outfile_path']}
	
	if [ -f {outputs['sfs_newrun_outfile_path']}]
	then
		echo Checked if outfile exist, and it does. File: {outputs['sfs_newrun_outfile_path']}
	else
		echo Outout file does not exist
	fi

	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def create_run_name_fsc(idx: str, target: AnonymousTarget) -> str:
	#pair_name=target.inputs['name_pops']
	return f'run_fsc_migr_{os.path.basename(target.outputs["parameter_value_likelihood_file"]).replace("-","_")}'


def setup_run_FSC_map_target(SFS_file: str, migration_divide: int, name_pops: str):
	"""
	Template: Setup folders and files for running FSC and start the run.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'SFS_file': SFS_file}
	outputs = {'parameter_value_likelihood_file': f'{os.path.dirname(inputs['SFS_file'])}/{migration_divide}/{os.path.basename(inputs['SFS_file']).replace("_sfs2d_","_" + str(migration_divide) + "migDiv_sfs2d_").replace("_jointMAFpop1_0.obs", "")}/{os.path.basename(inputs['SFS_file']).replace("_sfs2d_","_" + str(migration_divide) + "migDiv_sfs2d_").replace(".obs",".bestlikelyhoods").replace("_jointMAFpop1_0", "")}'} 
	options = {
		'cores': 1,
		'memory': '1g',
		#'walltime': '11:00:00'
		'walltime': '1-12:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate migration_fsc

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	# Make files and folders ready:
	
	# get folder name
	diris={os.path.dirname(inputs['SFS_file'])}
	echo Input directory of SFS is: $diris
	
	# make subfolder for specific migration divide
	mkdir -p $diris/{migration_divide}
	mig_div_var={migration_divide}

	# Copy obs sfs and modify name to include migration divide
	obs_file=`ls $diris/*sfs2d_jointMAFpop1_0.obs`
	obs_file_base=`basename $obs_file`
	echo SFS is moved and renamed for new run like this: $diris/{migration_divide}/${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}}

	if [ ! -f $diris/{migration_divide}/${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}} ]
	then
		echo .obs file did not exist in correct folder, copying and renaming.
		cp $obs_file $diris/{migration_divide}/${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}}
	fi
	ls $diris/{migration_divide}/
	#cp $diris/*sfs2d_jointMAFpop1_0.obs $diris/{migration_divide}/
	
	
	# Copy est and tpl files
	est_file=../../../workflow_source/scripts/pop1_pop2_sfs2d.est
	tpl_file=../../../workflow_source/scripts/pop1_pop2_sfs2d.tpl
	est_base=`basename $est_file`
	tpl_base=`basename $tpl_file`
	
	# modify them and change names
	if [ ! -f $diris/{migration_divide}/${{est_base/"_sfs2d"/"_"$mig_div_var"migDiv_sfs2d"}} ]
	then
		cp $est_file $diris/{migration_divide}/${{est_base/"_sfs2d"/"_"$mig_div_var"migDiv_sfs2d"}}
	fi
	if [ ! -f $diris/{migration_divide}/${{tpl_base/"_sfs2d"/"_"$mig_div_var"migDiv_sfs2d"}} ]
	then
		cp $tpl_file $diris/{migration_divide}/${{tpl_base/"_sfs2d"/"_"$mig_div_var"migDiv_sfs2d"}}
	fi

	
	# change migration tpl
	sed -i -e 's/600/{migration_divide}/g' $diris/{migration_divide}/${{tpl_base/"_sfs2d"/"_"$mig_div_var"migDiv_sfs2d"}}

	# now fix pop names:
	rename "pop1_pop2" {name_pops} $diris/{migration_divide}/pop1_pop2*


	# run FastSimCoal
	cd $diris/{migration_divide}
	#/faststorage/project/EcoGenetics/people/anneaa/programs/fastsimcoal/fsc28_linux64/fsc28 -t *.tpl -e *.est -n 100000 -m -M -L 40 -c 2
	/faststorage/project/EcoGenetics/people/anneaa/programs/fastsimcoal/fsc28_linux64/fsc28 -t *.tpl -e *.est -n 100000 -m -M -L 40 -y 3

	echo FastSimCoal run has finished
	mv {outputs['parameter_value_likelihood_file'].replace(".bestlikelyhoods", ".bestlhoods")} {outputs['parameter_value_likelihood_file']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def create_run_name_fsc_noSingl(idx: str, target: AnonymousTarget) -> str:
	#pair_name=target.inputs['name_pops']
	return f'run_fsc_migr_NoSingl_{os.path.basename(target.outputs["parameter_value_likelihood_file"]).replace("-","_")}'


def setup_run_FSC_map_target_nosingletons(SFS_file: str, migration_divide: int, name_pops: str):
	"""
	Template: Setup folders and files for running FSC and start the run.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'SFS_file': SFS_file}
	outputs = {'parameter_value_likelihood_file': f'{os.path.dirname(inputs['SFS_file'])}/{migration_divide}/{os.path.basename(inputs['SFS_file']).replace("_sfs2d_","_" + str(migration_divide) + "migDiv_sfs2d_").replace("_jointMAFpop1_0.obs", "")}/{os.path.basename(inputs['SFS_file']).replace("_sfs2d_","_" + str(migration_divide) + "migDiv_sfs2d_").replace(".obs",".bestlikelyhoods").replace("_jointMAFpop1_0", "")}'} 
	options = {
		'cores': 1,
		'memory': '1g',
		#'walltime': '11:00:00'
		'walltime': '1-12:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate migration_fsc

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	# Make files and folders ready:
	
	# get folder name
	diris={os.path.dirname(inputs['SFS_file'])}
	echo Input directory of SFS is: $diris
	
	# make subfolder for specific migration divide
	mkdir -p $diris/{migration_divide}
	mig_div_var={migration_divide}

	# Copy obs sfs and modify name to include migration divide
	obs_file=`ls $diris/*sfs2d_jointMAFpop1_0.obs`
	obs_file_base=`basename $obs_file`
	echo SFS is moved and renamed for new run like this: $diris/{migration_divide}/${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}}

	if [ ! -f $diris/{migration_divide}/${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}} ]
	then
		echo .obs file did not exist in correct folder, copying and renaming.
		cp $obs_file $diris/{migration_divide}/${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}}
	fi
	ls $diris/{migration_divide}/
	#cp $diris/*sfs2d_jointMAFpop1_0.obs $diris/{migration_divide}/
	
	
	# Copy est and tpl files
	est_file=../../../workflow_source/scripts/pop1_pop2_sfs2d.est
	tpl_file=../../../workflow_source/scripts/pop1_pop2_sfs2d.tpl
	est_base=`basename $est_file`
	tpl_base=`basename $tpl_file`
	
	# modify them and change names
	if [ ! -f $diris/{migration_divide}/${{est_base/"_sfs2d"/"_NoSingl_"$mig_div_var"migDiv_sfs2d"}} ]
	then
		cp $est_file $diris/{migration_divide}/${{est_base/"_sfs2d"/"_NoSingl_"$mig_div_var"migDiv_sfs2d"}}
	fi
	if [ ! -f $diris/{migration_divide}/${{tpl_base/"_sfs2d"/"_NoSingl_"$mig_div_var"migDiv_sfs2d"}} ]
	then
		cp $tpl_file $diris/{migration_divide}/${{tpl_base/"_sfs2d"/"_NoSingl_"$mig_div_var"migDiv_sfs2d"}}
	fi

	
	# change migration tpl
	sed -i -e 's/600/{migration_divide}/g' $diris/{migration_divide}/${{tpl_base/"_sfs2d"/"_NoSingl_"$mig_div_var"migDiv_sfs2d"}}

	# now fix pop names:
	rename "pop1_pop2" {name_pops} $diris/{migration_divide}/pop1_pop2*


	# run FastSimCoal
	cd $diris/{migration_divide}
	#/faststorage/project/EcoGenetics/people/anneaa/programs/fastsimcoal/fsc28_linux64/fsc28 -t *.tpl -e *.est -n 100000 -m -M -L 40 -c 2
	/faststorage/project/EcoGenetics/people/anneaa/programs/fastsimcoal/fsc28_linux64/fsc28 -t *.tpl -e *.est -n 100000 -m -M -L 40 -y 3 --nosingleton

	echo FastSimCoal run has finished
	mv {outputs['parameter_value_likelihood_file'].replace(".bestlikelyhoods", ".bestlhoods")} {outputs['parameter_value_likelihood_file']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

