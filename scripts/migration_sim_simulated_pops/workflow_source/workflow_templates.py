#!/bin/env python3
from gwf import AnonymousTarget
import os, glob, yaml
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


def python_make_yml(SFS_file: str, describer: str):
	# Define input and output file paths
	inputs = SFS_file
	# {'SFS_file': 'path/to/input_sfs_file'}
	output_file = f'{os.path.dirname(inputs)}/{os.path.basename(inputs).replace(".obs", ".yaml")}'
	# {'new_sfs_list_file': 'path/to/output_file.yaml'}

	# Step 1: Read the first field of the first line to get `obs_count`
	with open(inputs, 'r') as sfs_file:
		first_line = sfs_file.readline()
		obs_count = int(first_line.split()[0])

	# Step 2: Prepare output file in YAML format
	with open(output_file, 'w') as out_file:
		out_file.write("# List of the 10 simulated SFS observations, as replicates of the simulated SFS.\n")
		out_file.write("replicate_sfs_list:\n")

	# Step 3: Check if there’s only one observation or more
	if obs_count == 1:
		print(f"Only one observation in input SFS file: {inputs}")
	else:
		print("Input SFS file has more than one observation. Making new names for split files.")

		# Step 4: Loop to create directories and file names for each observation
		sfs_dir = os.path.dirname(inputs)
		basename_file = os.path.basename(inputs)
		describer_is = describer  # Replace this with the actual describer pattern

		with open(output_file, 'a') as out_file:
			for obs_n in range(1, obs_count + 1):
				#print(obs_n)
				
				# Generate new path
				rep_name = f"rep_{obs_n}_{describer_is}"
				new_path = f"{sfs_dir}/rep_{obs_n}/{basename_file.replace(describer_is, rep_name)}"
				
				# Write to YAML file
				out_file.write(f"  - {new_path}\n")
				
				# Create directory
				os.makedirs(os.path.join(sfs_dir, f"rep_{obs_n}"), exist_ok=True)




def split_file(input_filename: str, path_out: str, file_name_base: str):
    with open(input_filename, 'r') as file:
        lines = file.readlines()

    # Retrieve the header from line 2
    header = lines[1]
    file_index = 1  # Initialize file counter

    # Initialize an empty list to collect lines for the current split file
    current_data = ["1 observations\n", header]  # Start with "1 observations" and header

    for line in lines[2:]:
        # Check if line starts with "d1_0", indicating a new section
        if line.startswith("d1_0"):
            # Only write the file if `current_data` contains more than "1 observations" and the header
            if len(current_data) > 2:
                # Create directory for each file, if not exists
                os.makedirs(f'{path_out}/rep_{file_index}', exist_ok=True)
                with open(f'{path_out}/rep_{file_index}/{file_name_base.replace("rep_X", "rep_" + str(file_index))}', 'w') as outfile:
                    outfile.writelines(current_data)
                file_index += 1  # Increment file index after writing the file

            # Reset current_data with the header for the new section
            current_data = ["1 observations\n", header, line]
        else:
            # Add line to the current data list
            current_data.append(line)

    # Write the remaining data to the last file if it has content
    if len(current_data) > 2:
        os.makedirs(f'{path_out}/rep_{file_index}', exist_ok=True)
        with open(f'{path_out}/rep_{file_index}/{file_name_base.replace("rep_X", "rep_" + str(file_index))}', 'w') as outfile:
            outfile.writelines(current_data)





def create_run_fsc(idx: str, target: AnonymousTarget) -> str:
	#pair_name=target.inputs['name_pops']
	return f'run_fsc_migr_{os.path.basename(target.outputs["parameter_value_likelihood_file"]).replace("-","_")}'


def setup_run_FSC_target(SFS_file: str, migration_divide: int, name_pops: str):
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
		'walltime': '11:59:00'
		#'walltime': '1-12:00:00'
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
	obs_file=`ls $diris/*sfs2d_*_jointMAFpop1_0.obs`
	obs_file_base=`basename $obs_file`
	echo Original observed file is: $obs_file
	echo SFS is moved and renamed for new run like this: $diris/{migration_divide}/${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}}

	if [ ! -f $diris/{migration_divide}/${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}} ]
	then
		echo .obs file did not exist in correct folder, copying and renaming.
		cp $obs_file $diris/{migration_divide}/${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}}
	else
		echo .obs file already exist in correct folder. $diris/{migration_divide}/${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}}
	fi
	# get name from sfs file
	new_sfs_base=${{obs_file_base/"_sfs2d_"/"_"$mig_div_var"migDiv_sfs2d_"}}
	sfs_inherit_name="${{new_sfs_base%%_joint*}}"
	echo $sfs_inherit_name
	
	
	# Copy est and tpl files
	est_file=../../../workflow_source/scripts/pop1_pop2_sfs2d.est
	tpl_file=../../../workflow_source/scripts/pop1_pop2_sfs2d.tpl
	est_base=`basename $est_file`
	tpl_base=`basename $tpl_file`
	
	# modify them and change names
	#if [ ! -f $diris/{migration_divide}/${{est_base/"_sfs2d"/"_"$mig_div_var"migDiv_sfs2d"}} ]
	if [ ! -f $diris/{migration_divide}/$sfs_inherit_name".est" ]
	then
		#cp $est_file $diris/{migration_divide}/${{est_base/"_sfs2d"/"_"$mig_div_var"migDiv_sfs2d"}}
		cp $est_file $diris/{migration_divide}/$sfs_inherit_name.est
	fi
	#if [ ! -f $diris/{migration_divide}/${{tpl_base/"_sfs2d"/"_"$mig_div_var"migDiv_sfs2d"}} ]
	if [ ! -f $diris/{migration_divide}/$sfs_inherit_name".tpl" ]
	then
		#cp $tpl_file $diris/{migration_divide}/${{tpl_base/"_sfs2d"/"_"$mig_div_var"migDiv_sfs2d"}}
		cp $tpl_file $diris/{migration_divide}/$sfs_inherit_name.tpl
	fi

	echo "\n Listing directory: $diris/{migration_divide}/"
	ls $diris/{migration_divide}/
	echo "\n"

	# change migration tpl
	#sed -i -e 's/600/{migration_divide}/g' $diris/{migration_divide}/${{tpl_base/"_sfs2d"/"_"$mig_div_var"migDiv_sfs2d"}}
	sed -i -e 's/600/{migration_divide}/g' $diris/{migration_divide}/$sfs_inherit_name.tpl
	echo migration divide changed to {migration_divide}

	# now fix pop names:	
	files_name_correct=`ls $diris/{migration_divide}/{name_pops}*|head -n1`
	if [ -f $files_name_correct ]
	then
		echo Files already have desired population names: {name_pops}
	else
		echo Renaming the files with new population names: "pop1_pop2" changed to {name_pops}
		rename "pop1_pop2" {name_pops} $diris/{migration_divide}/pop1_pop2*
	fi


	# run FastSimCoal
	cd $diris/{migration_divide}
	pwd
	echo Ready to run FastSimCoal
	/home/anneaa/EcoGenetics/people/anneaa/programs/fastsimcoal/fsc28_linux64/fsc28 -t *.tpl -e *.est -n 100000 -m -M -L 40 -y 3
	#-y 3 is default
	
	echo FastSimCoal run has finished
	mv {outputs['parameter_value_likelihood_file'].replace(".bestlikelyhoods", ".bestlhoods")} {outputs['parameter_value_likelihood_file']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



