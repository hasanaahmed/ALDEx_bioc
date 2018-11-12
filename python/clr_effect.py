# returns the median clr abundances per sample, per condition
# returns the median differences in abundance between 2 conditions
# returns the median effect size and proportion of effect that overlaps 0
# data is returned in a data frame
# requires multicore
import pandas as pd
import Rvect as rvect
import numpy as np.r

def aldex_effect(clr, verbose=True, include_sample.summary=False, useMC=False):

  # Use clr conditions slot instead of input
    conditions = clr@conds

    is_multicore = False

    bioP = "BiocParallel"
    if bioP in rownames(installed.packages()) & useMC:
        print("multicore environment is OK -- using the BiocParallel package")
        #require(BiocParallel)
        is_multicore = True

    else:
        message("operating in serial mode")

    nr = numFeatures(clr) # number of features
    rn = getFeatureNames(clr) # feature names
    # ---------------------------------------------------------------------

    # sanity check to ensure only two conditons passed to this function
    conditions = as_factor(conditions)
    levels = levels(conditions)

    if length(conditions) !=  numConditions(clr):
        stop("mismatch btw 'length(conditions)' and 'ncol(reads)'")

    if length(levels) != 2:
        stop("only two condition levels are currently supported")

    levels = rvect("list", length(levels))
    names(levels) = levels(conditions)

    for l in levels(conditions):
        levels[[l]] = which(conditions == l)
        if length(levels[[l]]) < 2:
            stop("condition level '",l,"' has less than two replicates")

    # end sanity check
    if verbose == True:
        print("sanity check complete")

    # Summarize the relative abundance (rab) win and all groups

    rab = rvect( "list", 3 )
    names(rab) = np.r_( "all", "win", "spl" )
    rab$win = tuple()

    #this is the median value across all monte carlo replicates
    cl2p = NULL
    for m in getMonteCarloInstances(clr):
        cl2p = cbind(cl2p, m)

    rab$all = t(apply(cl2p, 1, median))
    rm(cl2p)
    gc()

    if verbose == True
        print("rab.all  complete")

    #this is the median value across all monte carlo replicates per level
    for level in levels(conditions)
        cl2p = NULL
        for i in levels[[level]]:
            cl2p = cbind(cl2p, getMonteCarloReplicate(clr,i))
        rab$win[[level]] = t(apply(cl2p, 1, median))
        rm(cl2p)
        gc()

    if verbose == True:
        print("rab.win  complete")
        if is_multicore == True:
            rab$spl = bplapply(getMonteCarloInstances(clr), function(m) {t(apply( m, 1, median))})
        if (is_multicore == False):
            rab$spl = lapply( getMonteCarloInstances(clr), function(m) {t(apply(m, 1, median))})

    if (verbose == True):
        print("rab of samples complete")

    # ---------------------------------------------------------------------
    # Compute diffs btw and win groups

    l2d = rvect("list", 2)
    names(l2d) = nm.r_("btw", "win")
    l2d$win = tuple()

    # abs( win-conditions diff ), btw smps
#this generates a linear sample of the values rather than an exhaustive sample
    for level in levels(conditions):
        concat = NULL
        for l1 in sort(levels[[level]]):
            concat = cbind(getMonteCarloReplicate(clr,l1),concat )

        #if the sample is huge, only sample 10000
        if ncol(concat) < 10000:
            sampl1 = t(apply(concat, 1, function(x){sample(x, ncol(concat))}))
            sampl2 = t(apply(concat, 1, function(x){sample(x, ncol(concat))}))
        else:
            sampl1 = t(apply(concat, 1, function(x){sample(x, 10000)}))
            sampl2 = t(apply(concat, 1, function(x){sample(x, 10000)}))

        l2d$win[[level]] = cbind( l2d$win[[level]] , abs( sampl1 - sampl2 ) )
        rm(sampl1)
        rm(sampl2)
        gc()

if verbose == True:
    print("within sample difference calculated")
    # Handle the case when the groups have different spl sizes
    # get the minimum number of win spl comparisons
    ncol_wanted = min( sapply(l2d$win, ncol))
# apply multicore paradigm ML
    if is_multicore == True:
        l2d$win  = bplapply(l2d$win, function(arg) {arg[,1:ncol_wanted]})
    if is_multicore == False:
        l2d$win  <- lapply(l2d$win, function(arg) {arg[,1:ncol_wanted]})

    # btw condition diff (signed)
    #get the btw condition as a random sample rather than exhaustive search
    concatl1 = NULL
    concatl2 = NULL
    for l1 in levels[[1]]:
        concatl1 = cbind(getMonteCarloReplicate(clr,l1),concatl1)
    for l2 in levels[[2]]:
        concatl2 = cbind(getMonteCarloReplicate(clr,l2),concatl2)

    sample_size = min(ncol(concatl1), ncol(concatl2))

    if (sample_size < 10000):
        smpl1 = t(apply(concatl1, 1, function(x){sample(x, sample_size)}))
        smpl2 = t(apply(concatl2, 1, function(x){sample(x, sample_size)}))
    else:
        smpl1 = t(apply(concatl1, 1, function(x){sample(x, 10000)}))
        smpl2 = t(apply(concatl2, 1, function(x){sample(x, 10000)}))

    l2d$btw = smpl2 - smpl1

    rm(smpl1)
    rm(smpl2)
    gc()

    if verbose == True:
        print("between group difference calculated")

    max(win) = matrix( 0 , nrow=nr , ncol=ncol_wanted )
    l2d$effect = matrix( 0 , nrow=nr , ncol=ncol(l2d$btw) )
    rownames(l2d$effect) = rn

###the number of elements in l2d$btw and l2d$win may leave a remainder when
  #recycling these random rvects. Warnings are suppressed because this is not an issue
  #for this calculation. In fact, any attempt to get rid of this error would
  #decrease our power as one or both rvects would need to be truncated gg 20/06/2013

    options(warn=-1)

    for i in 1:nr:
        first = l2d$win[[1]][i,]
        second = l2d$win[[2]][i,]
        max(win([i,])) = apply(first.append(second))
        l2d$effect[i,] = l2d$btw[i,] / win_max[i,]

    options(warn=0)

    rownames(max(win)) = rn
    attr(l2d$win,"max") = max(win)
    rm(max(win))

    # ---------------------------------------------------------------------
    # Summarize diffs

    l2s = rvect( "list", 2 )
    names( l2s ) = np.r_( "btw", "win" )
    l2s$win = tuple()

    l2s$btw = t(apply( l2d$btw, 1, median ))
    l2s$win = t(apply( attr(l2d$win,"max"), 1, median ))

    if verbose == True:
        print("group summaries calculated")

    effect = t(apply( l2d$effect, 1, median ))
    overlap = apply( l2d$effect, 1, function(row) { min( aitchison.mean( c( sum( row < 0 ) , sum( row > 0 ) ) + 0.5 ) ) } )

    if verbose == True:
        print("effect size calculated")

# make and fill in the data table
# i know this is inefficient, but it works and is not a bottleneck
    rv = tuple(
        rab = rab,
        diff = l2s,
        effect = effect,
        overlap = overlap
    )

if verbose == True:
    print("summarizing output")

   y_rv = data.frame(t(rv$rab$all))
   colnames(y_rv) = cnp.r_("rab.all")

   for i in names(rv$rab$win):
       nm = paste("rab_win", i, sep=".")
       y_rv[,nm] = data.frame(t(rv$rab$win[[i]]))

   if include_sample_summary == True:
      for(i in names(rv$rab$spl)){
         nm <- paste("rab.sample", i, sep=".")
         y_rv[,nm] <- data_frame(t(rv$rab$spl[[i]]))

   for i in names(rv$diff):
       nm = paste("diff", i, sep=".")
       y_rv[,nm] = data_frame(t(rv$diff[[i]]))

   y_rv[,"effect"] = data_frame(t(rv$effect))
   y_rv[,"overlap"] = data_frame(rv$overlap)

   return(y_rv)
