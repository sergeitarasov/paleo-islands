# Run from the project root: Rscript smaps/run-smaps-posterior-parallel.R pilot
# After the pilot passes, replace pilot with full.
# Six local workers; each replicate uses matching model and rates in both clades.
args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args)) args[1] else "pilot"
stopifnot(mode %in% c("pilot","full"))
source("smaps/setup-smaps.R")
# Freeze source identities so resumed runs cannot mix package/input revisions.
package.files <- list.files(file.path(find.package("divextra"),"R"),full.names=TRUE)
package.md5 <- tools::md5sum(package.files)
names(package.md5) <- basename(package.files)
# Completed maps from before the interpolation-check fix remain valid.
# Their cache had already passed the unchanged accuracy thresholds.
previous.package.md5 <- NULL
if (file.exists(file.path(run.dir,"validated-previous-package.rds")))
  previous.package.md5 <- readRDS(file.path(run.dir,"validated-previous-package.rds"))
input.files <- c("asr/results/sisyphini.gr6.rds","asr/results/sisyphini.gr7.rds",
 "asr/results/nesovin.gr6.rds","asr/results/nesovin.gr7.rds","asr/results/model-weights.rds",
 list.files("mcmc/inputs",full.names=TRUE),
 list.files("mcmc/results/roihu-5000-20261023",pattern="^gr[67]-chain-[1-4]\\.rds$",full.names=TRUE))
input.md5 <- tools::md5sum(input.files)
if (mode == "full") {
  pilot <- readRDS("smaps/runs/posterior-aic-saved-asr-pilot/pilot-success.rds")
  stopifnot(identical(pilot$package.md5,package.md5), identical(pilot$input.md5,input.md5))
}
saveRDS(list(package.md5=package.md5,input.md5=input.md5),file.path(run.dir,"run-provenance.rds"))

plan <- data.frame(model=c(rep(6L,N.GR6),rep(7L,N.GR7)),
                   map=c(seq_len(N.GR6),seq_len(N.GR7)))
# Interleave models so all workers cover both; no preferential replacement of failures.
plan <- plan[order(plan$map,plan$model),]
if (mode == "pilot") {
  # Original posterior rows/seeds, including all three previously failed draws.
  plan <- data.frame(model=c(6L,7L,6L,7L,7L,6L),map=c(1L,1L,20L,18L,46L,2L))
  # A pilot rerun must regenerate maps, not merely reuse the prior pilot.
  marker <- file.path(run.dir,"pilot-success.rds")
  if (file.exists(marker)) unlink(marker)
}
write.csv(plan,file.path(run.dir,"map-plan.csv"),row.names=FALSE)
workers <- vector("list",6)
for (worker in 1:6) {
  workers[[worker]] <- parallel::mcparallel(tryCatch({
    for (j in seq(worker,nrow(plan),by=6)) {
      model <- plan$model[j]
      i <- plan$map[j]
      if (model == 6L) {
        row <- rows.nesosisyphus.gr6[i]
        rates <- posterior.gr6[row,]
        full.pars <- lik.gr6(rates,pars.only=TRUE)
      } else {
        row <- rows.nesosisyphus.gr7[i]
        rates <- posterior.gr7[row,]
        full.pars <- lik.gr7(rates,pars.only=TRUE)
      }
      for (clade in c("nesosisyphus","nesovin")) {
        filename <- file.path(maps.dir,sprintf("%s-gr%d-%04d.rds",clade,model,i))
        if (mode == "full" && file.exists(filename)) {
          previous <- readRDS(filename)
          stopifnot(identical(previous$posterior.draw$rates,rates),
                    previous$posterior.draw$posterior.row == row,
                    (identical(previous$posterior.draw$package.md5, package.md5) ||
                     identical(previous$posterior.draw$package.md5, previous.package.md5)),
                    identical(previous$posterior.draw$input.md5, input.md5))
          next
        }
        if (clade == "nesosisyphus") {
          tree <- sisyphini
          tip.weights <- if (model == 6L) sisyphini.gr6$tip.states else sisyphini.gr7$tip.states
          tips <- tips.nesosisyphus
          root.prior <- if (model == 6L) root.nesosisyphus.gr6 else root.nesosisyphus.gr7
          seed <- (model-5)*1000000+i
        } else {
          tree <- nesovin
          tip.weights <- if (model == 6L) nesovin.gr6$tip.states else nesovin.gr7$tip.states
          tips <- tips.nesovin
          root.prior <- if (model == 6L) root.nesovin.gr6 else root.nesovin.gr7
          seed <- (model-3)*1000000+i
        }
        message("Starting ",clade," gr",model," map ",i," posterior row ",row)
        timing <- system.time(sim <- make.simmap.classe.td(tree=tree,pars=full.pars,
          tip.weights=tip.weights,k=16,state.labels=state.labels,
          backbone.method="tree-conditioned",root.prior=root.prior,
          side.lineage.mode="extinct_clades",extinction.sampler="direct",side.lineage.tips=tips,
          expand.side.lineages=TRUE,
          max.taxa=20000,max.events=200000,nsim=1,seed=seed))
        stopifnot(all(sim$tips$fate[sim$tips$generated] == "extinct"))
        sim <- keep.lineage.simmap(sim,tips=tips)
        stopifnot(setequal(sim$tips$label[!sim$tips$generated],tips),
                  max(abs(vapply(sim$tree$maps,sum,0)-sim$tree$edge.length)) < 1e-7)
        # Also validates root-map and mapped-segment durations.
        lineage.intervals.simmap(sim)
        sim$posterior.draw <- list(model=model,map=i,posterior.row=row,
          chain=(row-1)%/%4000+1,sweep=(row-1)%%4000+1001,rates=rates,
          full.pars=full.pars,root.prior=root.prior,seed=seed,
          seconds=unname(timing["elapsed"]),package.md5=package.md5,input.md5=input.md5,divextra.version=as.character(packageVersion("divextra")))
        saveRDS(sim,paste0(filename,".tmp"))
        stopifnot(file.rename(paste0(filename,".tmp"),filename))
        message("Completed ",clade," gr",model," map ",i," in ",timing["elapsed"]," seconds")
      }
    }
    TRUE
  }, error=function(e) {
    writeLines(conditionMessage(e),file.path(run.dir,paste0("worker-",worker,"-error.txt")))
    message("Worker ",worker,": ",conditionMessage(e))
    stop(e)
  }),mc.set.seed=FALSE)
}
status <- parallel::mccollect(workers)
saveRDS(status,file.path(run.dir,"worker-status.rds"))
for (worker in 1:6) stopifnot(isTRUE(status[[worker]]))

# Pilot: report all 12 maps and record the exact inputs/package that passed ----
if (mode == "pilot") {
  report <- data.frame()
  for (j in seq_len(nrow(plan))) {
    for (clade in c("nesosisyphus","nesovin")) {
      sim <- readRDS(file.path(maps.dir,sprintf("%s-gr%d-%04d.rds",clade,plan$model[j],plan$map[j])))
      report <- rbind(report,data.frame(clade=clade,model=plan$model[j],map=plan$map[j],
        posterior.row=sim$posterior.draw$posterior.row,seconds=sim$posterior.draw$seconds,
        observed=sum(!sim$tips$generated),extinct=sum(sim$tips$fate=="extinct")))
    }
  }
  stopifnot(nrow(report)==12)
  summary.dir <- "smaps/summaries/posterior-aic-saved-asr-pilot"
  dir.create(summary.dir,recursive=TRUE,showWarnings=FALSE)
  write.csv(report,file.path(summary.dir,"pilot-results.csv"),row.names=FALSE)
  saveRDS(list(package.md5=package.md5,input.md5=input.md5,report=report),
          file.path(run.dir,"pilot-success.rds"))
  print(report)
  message("PILOT PASSED. Start the full analysis with: bash smaps/run-smaps-local.sh full")
  quit(save="no")
}

# Collect models separately, then combine using their AIC allocations ----
for (clade in c("nesosisyphus","nesovin")) {
  combined <- list()
  for (model in c(6L,7L)) {
    n <- if (model == 6L) N.GR6 else N.GR7
    maps <- vector("list",n)
    for (i in seq_len(n)) {
      maps[[i]] <- readRDS(file.path(maps.dir,sprintf("%s-gr%d-%04d.rds",clade,model,i)))
    }
    saveRDS(maps,file.path(maps.dir,sprintf("smaps-%s-gr%d.rds",clade,model)))
    combined <- c(combined,maps)
  }
  saveRDS(combined,file.path(maps.dir,paste0("smaps-",clade,"-full.rds")))
  # Observed-only copies for subsequent arrival-time analyses.
  observed <- vector("list",length(combined))
  for (i in seq_along(combined)) observed[[i]] <- prune.unobserved.simmap(combined[[i]])
  saveRDS(observed,file.path(maps.dir,paste0(clade,"-observed.rds")))
}
source("smaps/plot-sunken-richness-posterior.R")
