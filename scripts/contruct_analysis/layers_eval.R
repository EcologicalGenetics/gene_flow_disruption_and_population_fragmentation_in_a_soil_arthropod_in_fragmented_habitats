#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Run spatial models with varying K to inspect layer contributions

# Example usage cmd: Rscript layers_eval.R mydata.Rdata

# Prep environgment
pacman::p_load("conStruct")
load(args[1]) # Load Rdata to load

#---- Running construct with different K to compare layer contributions
# https://cran.r-project.org/web/packages/conStruct/vignettes/model-comparison.html


#setwd as construct saves pictures/files there
setwd("/faststorage/project/ostrich_thermal/people/madsfs/ecogenetics/entnic/construct/output/layer400k")


my.run_K1 <- conStruct(spatial = TRUE, 
                    K = 1, 
                    freqs = Ent_Freq, 
                    geoDist = Ent_geo, 
                    coords = Ent_GPS,
                    prefix = "K1_LayEval",
                    save.files = TRUE,
                    make.figs = TRUE,
                    n.chains = 3,
                    n.iter = 400000) 

my.run_K2 <- conStruct(spatial = TRUE, 
                    K = 2, 
                    freqs = Ent_Freq, 
                    geoDist = Ent_geo, 
                    coords = Ent_GPS,
                    prefix = "K2_LayEval",
                    save.files = TRUE,
                    make.figs = TRUE,
                    n.chains = 3,
                    n.iter = 400000)

my.run_K3 <- conStruct(spatial = TRUE, 
                    K = 3, 
                    freqs = Ent_Freq, 
                    geoDist = Ent_geo, 
                    coords = Ent_GPS,
                    prefix = "K3_LayEval",
                    save.files = TRUE,
                    make.figs = TRUE,
                    n.chains = 3,
                    n.iter = 400000)

my.run_K4 <- conStruct(spatial = TRUE, 
                    K = 4, 
                    freqs = Ent_Freq, 
                    geoDist = Ent_geo, 
                    coords = Ent_GPS,
                    prefix = "K4_LayEval",
                    save.files = TRUE,
                    make.figs = TRUE,
                    n.chains = 3,
                    n.iter = 400000)

my.run_K5 <- conStruct(spatial = TRUE, 
                    K = 5, 
                    freqs = Ent_Freq, 
                    geoDist = Ent_geo, 
                    coords = Ent_GPS,
                    prefix = "K5_LayEval",
                    save.files = TRUE,
                    make.figs = TRUE,
                    n.chains = 3,
                    n.iter = 400000)

