#!/bin/env python3
import sys, os
sys.path.insert(0, os.path.realpath('../../../workflow_source/'))
from workflow_source import *

gwf = migration_simulation_simsfs()


# creates 2dSFS and estimates migration rates based on those using FastSimCoal


# conda create -n migration_fsc python R gwf pyyaml r-ggplot2 r-viridis

# OBS: NEED TO CHECK IF PAIRS ARE CORRECT AFTER NEW POPS ADDED? Do they still get found correctly in the VCF that assumedly only has the "old pops"?
# ALSO: Check which populations fail, and remove them (via the not run list), since they most likely are new anyway
# cd /home/anneaa/EcoGenetics/general_workflows/population_genetics/migration_sim_simulated_pops/configurations/collembola/Entomobrya_nicoleti_simulated
# conda activate migration_fsc
# gwf run
# gwf status
# gwf status -f summary
# gwf logs pi_remodelling_file
# gwf logs --stderr pair_2DSFS_map_target_519
#⨯ pi_calculation_all_positions failed
#⨯ pi_add_context cancelled

#⨯ shouldrun         0
#- submitted     23252
#↻ running          54
#✓ completed     15843
#⨯ failed            0
#⨯ cancelled         0

#36 min 2 cores 1 gb
#1t min 2 cores 1 gb
# similar with 1 core. huge time difference between jobs. some with 2 cores are still running 8h, while most finish quickly with both 1 or 2 coeres.
# could use -y parameter, indication that after X rounds of not inproving likelyhood, go back to the previous parameters improving likelyhood and retry.

#FastSimCoal: new analysis
#Gene flow:
#change:
#generation time:#
#	intervals of 100 from 0 to split time.

#look at suggested references for similar SFS for different demographic history.

#compare simulated SFS with the observed one. 
#distribution of site frequency spectres?

#more clear in methodology about poolsec
#- coverage of each cromosome
