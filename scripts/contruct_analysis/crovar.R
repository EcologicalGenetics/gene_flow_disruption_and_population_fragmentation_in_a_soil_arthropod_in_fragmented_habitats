#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Example usage cmd: Rscript crovar.R mydata.Rdata Ent_Nic_K3

# Prep environgment
pacman::p_load("conStruct")
load(args[1]) # Load Rdata to load
myprefix = args[2] # Set prefix


#---- Running construct with cross validation
# https://cran.r-project.org/web/packages/conStruct/vignettes/model-comparison.html


#setwd as construct saves pictures/files there
setwd("/faststorage/project/ostrich_thermal/people/madsfs/ecogenetics/entnic/construct/output/croval400k")

my.xvals <- x.validation(train.prop = 0.9,
                         n.reps = 16,
                         K = 1:5,
                         freqs = Ent_Freq,
                         data.partitions = NULL,
                         geoDist = Ent_geo,
                         coords = Ent_GPS,
                         prefix = myprefix,
                         n.iter = 400000,
                         make.figs = TRUE,
                         save.files = FALSE,
                         parallel = TRUE,
                         n.nodes = 16)


#Format results from the output list
sp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$sp)}),init=NULL)
nsp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$nsp)}),init=NULL)


# first, get the 95% confidence intervals for the spatial and nonspatial
#   models over values of K (mean +/- 1.96 the standard error)

sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
