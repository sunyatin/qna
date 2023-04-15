require(ggplot2); require(ggdist); require(reshape2); require(viridis); require(paletteer);
require(scico); require(scales); require(ggbeeswarm); require(ggbump); require(dplyr); require(ggmap)
options(dplyr.summarise.inform = FALSE)

model_comparison <- function(Model.Paths, # list
                                    Model.Runs, # list
                                    Model.Names, # list
                                    N, # number of sims to accept | can be NA
                                    distance.metrics.files, # vector of files with stat distance metrics & weights
                                    distance.metrics.comparison, # single spar file of stats with which to compare models
                                    #
                                    generation.time,
                                    Sprime.min.match.rate,
                                    psmc..xlim,
                                    #
                                    outpng = NA,
                                    do.scaled.distance = TRUE, # scale the distances? (not useful if already using distances of scaled values)
                                    functions.of.summary = c("min", "mean", "median"),
                                    function.of.dispersion = "mad",
                                    verbose = FALSE,
                                    import.AFS.2D = TRUE,
                                    abc.tolerance = 0.1) 
{
  # NB:
  # "split.by.statistics", "all.raw.stats" & "relative.errors" elements in the output list will be EMPTY if there is more than one `distance.metrics.files`
  
  # Output
  #> "results"                    final results of the analysis: scaled distances per model average across stats
  #> "accepted.runs"              list of the accepted runs for each model
  #> "split.by.statistics"        scaled distances per model NOT averaged across the stats  /!\ only if single distance.metrics.file provided!
  #                                 NOTE! this is a list, whose size is equal to the number of elements in `functions.of.summary`
  #> "all.raw.stats"              values of the stats across each model and each run  /!\ only if single distance.metrics.file provided!
  #> "relative.percentage.errors" relative % errors of the stats across each model and each run  /!\ only if single distance.metrics.file provided!
  
  cat("Warning. Distances are NOT weighted for model comparison ss.\n")
  if (length(Model.Paths) != length(Model.Runs)) stop("Model.Paths & Model.Runs must have equal length")
  if (length(Model.Runs) != length(Model.Names)) stop("Model.Names & Model.Runs must have equal length")
  Res <- c()
  RunsAccepted <- vector("list", length(distance.metrics.files))
  names(RunsAccepted) <- distance.metrics.files
  
  #--------------------------------- Import ---------------------------------#
  
  Model.Data <- vector("list", length(Model.Paths))
  names(Model.Data) <- unlist(Model.Names)
  for (i in seq_along(Model.Paths)) {
    if (verbose) cat("\n\n================================================================> ")
    if (verbose) cat(Model.Names[[i]],"\n")
    Model.Paths[[i]] <- gsub("\\/$", "", Model.Paths[[i]], perl = TRUE)
    if (verbose) {
      Model.Data[[Model.Names[[i]]]] <- abc.import(Model.Paths[[i]],
                                                   files = Model.Runs[[i]],
                                                   do.AFS.2D = import.AFS.2D, 
                                                   read.params = F, 
                                                   generation.time = generation.time, 
                                                   Sprime.min.match.rate = Sprime.min.match.rate)
    } else {
      invisible(capture.output(Model.Data[[Model.Names[[i]]]] <- abc.import(Model.Paths[[i]],
                                                                            files = Model.Runs[[i]],
                                                                            do.AFS.2D = import.AFS.2D, 
                                                                            read.params = F, 
                                                                            generation.time = generation.time, 
                                                                            Sprime.min.match.rate = Sprime.min.match.rate)))
    }
  }
  
  cat("Number of sims in RT model: ",length(Model.Data[["z. This study"]]$NAME),"\n")
  
  #--------------------------------- Analyze ---------------------------------#
  
  # iterate over all the distance metric files
  for (spar_idx in seq_along(distance.metrics.files)) {
    if (verbose) cat("\n\n================================================================> ")
    if (verbose) cat(distance.metrics.files[spar_idx],"\n")
    cat(" ",distance.metrics.files[spar_idx],"\n")

    DISTANCES <- c()
    MODELS <- c()
    STATS <- c()
    RE <- c() # Relative Error for all stats
    RunsAccepted[[spar_idx]] <- vector("list", length(Model.Paths))
    names(RunsAccepted[[spar_idx]]) <- Model.Names
    for (i in seq_along(Model.Paths)) {
      Model <- Model.Data[[Model.Names[[i]]]]
      # SELECT RUNS
      invisible(capture.output(ABC <- abc.analyze(obs.data, 
                                                  Model, 
                                                  parfile = distance.metrics.files[spar_idx], 
                                                  N = N, 
                                                  dispersion.fun = as.character(function.of.dispersion), 
                                                  psmc..xlim = psmc..xlim,
                                                  verbose = FALSE)))
      RunsAccepted[[spar_idx]][[i]] <- ABC$ACCEPTED$name
      # CALCULATE ALL DISTANCES FOR MODEL COMPARISON
      if (N==1) {
        invisible(capture.output(ABC_COMPARE <- abc.analyze(obs.data, 
                                                            Model, 
                                                            parfile = distance.metrics.comparison, 
                                                            N = N, 
                                                            dispersion.fun = as.character(function.of.dispersion), 
                                                            psmc..xlim = psmc..xlim,
                                                            verbose = FALSE)))
        Acc <- ABC_COMPARE$DISTANCES$raw
        Acc <- Acc[ABC$ACCEPTED$index,,drop=F]
        if (length(spar_idx)==1) {
          STATS <- rbind(STATS, data.frame(model = Model.Names[[i]], 
                                           run = ABC$ACCEPTED$name,
                                           ABC_COMPARE$sumstats.merged.sim[ABC$ACCEPTED$index,,drop=F]))
          re <- sweep(ABC_COMPARE$sumstats.merged.sim[ABC$ACCEPTED$index,,drop=F], 2, ABC_COMPARE$sumstats.merged.obs, FUN = "-")
          re <- sweep(re, 2, ABC_COMPARE$sumstats.merged.obs, FUN = "/")
          RE <- rbind(RE, data.frame(model = Model.Names[[i]], run = ABC$ACCEPTED$name, re))
        }
      } else {
        invisible(capture.output(ABC_COMPARE <- abc.analyze(obs.data, 
                                                            extract(Model, ABC$ACCEPTED$index), 
                                                            parfile = distance.metrics.comparison, 
                                                            N = N, 
                                                            dispersion.fun = as.character(function.of.dispersion), 
                                                            psmc..xlim = psmc..xlim,
                                                            verbose = FALSE)))
        Acc <- ABC_COMPARE$DISTANCES$raw
        if (length(spar_idx)==1) {
          STATS <- rbind(STATS, data.frame(model = Model.Names[[i]], 
                                           run = ABC$ACCEPTED$name,
                                           ABC_COMPARE$sumstats.merged.sim))
          re <- sweep(ABC_COMPARE$sumstats.merged.sim, 2, ABC_COMPARE$sumstats.merged.obs, FUN = "-")
          re <- sweep(re, 2, ABC_COMPARE$sumstats.merged.obs, FUN = "/")
          RE <- rbind(RE, data.frame(model = Model.Names[[i]], run = ABC$ACCEPTED$name, re))
        }
      }
      DISTANCES <- rbind(DISTANCES, Acc)
      MODELS <- c(MODELS, rep(Model.Names[[i]], nrow(Acc)))
    }
    
    # 
    UNSCALED.DISTANCES <- data.frame(model = MODELS, DISTANCES)
    
    #--------------------------------- Scale? ---------------------------------#

    if (do.scaled.distance) {
      Means <- colMeans(DISTANCES, na.rm = TRUE)
      Sds <- apply(DISTANCES, 2, rlang::as_function(function.of.dispersion), na.rm = TRUE)
      if (any(Sds==0)) stop("Problem! Some dispersion values are null. If you used the MAD, try SD instead.")
      NORM.D <- sweep(DISTANCES, 2, Means, FUN = "-")
      NORM.D <- sweep(NORM.D, 2, Sds, FUN = "/")
      NORM.OBS <- 1/Sds * (rep(0, length(Means)) - Means)
      DISTANCES <- abs(sweep(NORM.D, 2, NORM.OBS, FUN = "-"))
    }
    
    #--------------------------------- Punctualize ---------------------------------#
    
    # checks
    if (verbose) {
      cat("\n\nNumber of simulations per model:\n")
      print(table(MODELS))
      cat("\n\n")
    }
    
    # ABC
    if (verbose) cat("\nABC: nRows: ",nrow(DISTANCES),"\n")
    PP <- suppressWarnings(abc::postpr(target =   rep(0, ncol(DISTANCES)),
                                       index =    MODELS,
                                       sumstat =  DISTANCES,
                                       tol =      abc.tolerance,
                                       method =   "mnlogistic"))
    PP <- summary(PP, print = F)
    if (names(PP)[1]=="rejection") {
      PP <- PP$mnlogistic$Prob
    } else {
      PP <- PP$Prob
    }
    
    # Punctualize
    BY.MODEL <- lapply(functions.of.summary, function(fun) aggregate(DISTANCES, by = list(model = MODELS), FUN = rlang::as_function(fun), na.rm = T))
    names(BY.MODEL) <- functions.of.summary
    
    #--------------------------------- Plot ---------------------------------#
    
    if (!is.na(outpng)) {
      
      DF <- reshape2::melt(data.frame(model = MODELS, DISTANCES), id.vars = "model")
      BP <- reshape2::melt(BY.MODEL, id.vars = "model")
      if (length(unique(BP$L1))!=1) stop("For plotting, you must specify a **single** summary function.")

      DF$model <- factor(as.character(DF$model), levels = sort(unlist(Model.Names)))
      BP$model <- factor(as.character(BP$model), levels = sort(unlist(Model.Names)))
      DF$variable <- factor(gsub(":", ".", as.character(DF$variable)), levels = gsub(":", ".", colnames(DISTANCES)))
      BP$variable <- factor(gsub(":", ".", as.character(BP$variable)), levels = gsub(":", ".", colnames(DISTANCES)))

      values <- c("#c43b60", "#aeaca0", "#006CA2", "#31c0d0", "#ACCDDC", "#005778", 
                  "black", "#bed8b3", "#335893", "#31dab0", "#908f88", "#3a3393")
      scg <- scale_color_manual(values = rev(values))
      scgf <- scale_fill_manual(values = rev(values))
      min.y <- min(DF$value, na.rm=T)*3
      
      g1=ggplot(BP, aes(x = as.numeric(variable), y = value, color = model)) +
        theme_ggdist() + ylab("Distance to observed") + xlab("") +
        geom_bump(size = 2, alpha = .8) + 
        geom_point(size = 2, col = "white") +
        geom_point(size = .3) +
        scg +
        scale_y_log10() +
        theme(axis.text.x = element_blank(),
              legend.position = "top",
              axis.line.x = element_blank()) +
        labs(color = "Models")
      
      g2=ggplot(DF) +
        theme_ggdist() + ylab("Distance to observed") + xlab("") +
        geom_leg(aes(x=variable, y=min.y, xend=variable, yend=Inf), col = "#EAEAEA", size = 25, lineend = "round") +
        geom_boxplot(aes(x=variable, y=value, fill=model, col=model), fatten = NULL, alpha = .8, position = position_dodge(.5), width = .3) +
        scg + scgf +
        scale_y_log10() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank())
      
      G <- cowplot::plot_grid(g1, g2, nrow = 2, rel_heights = c(.5, .5))
      invisible(capture.output(suppressWarnings(ggsave(G, filename = outpng, width = 15, height = 9))))
      
    }
    
    
    #--------------------------------- Statistics ---------------------------------#
    
    # Summarize each model across statistics
    SplitByStats <- list()
    res <- c()
    ## Iterating over all the functions.of.summary
    for (i in seq_along(BY.MODEL)) {
      X <- BY.MODEL[[i]]
      rownames(X) <- X[,1]
      if (length(spar_idx)==1) SplitByStats <- c(SplitByStats, list(X[,-1,drop=F]))
      pp <- c(PP[X[,1]])
      # ranks
      X <- data.frame(spar = distance.metrics.files[spar_idx],
                        fun.selection = names(BY.MODEL)[i],
                        model = X[,1], 
                      # summary across statistics:
                        avg.rank = rowMeans(apply(X[,-1,drop=F], 2, rank) / nrow(X)),
                        mean = apply(X[,-1,drop=F], 1, mean), 
                        median = apply(X[,-1,drop=F], 1, median),
                        abc.pp = pp)
      #X$median.rank <- rank(X$median) / nrow(X)
      X$median.ratio2best <- X$median / min(X$median)
      #X$mean.rank <- rank(X$mean) / nrow(X)
      X$mean.ratio2best <- X$mean / min(X$mean)
      rownames(X) <- NULL
      res <- rbind(res, X)
    }
    
    Res <- rbind(Res, res)
    
  }

  if (length(spar_idx)==1) {
    RE <- RE
    RE[,-c(1,2)] <- apply(RE[,-c(1,2)], 2, function(x) return(round(100*x, 1)))
  }
  return(list("results"=Res, 
              "accepted.runs"=RunsAccepted, 
              "split.by.statistics"=SplitByStats, 
              "all.raw.stats"=STATS, 
              "relative.percentage.errors"=RE,
              "unscaled.distances"=UNSCALED.DISTANCES))
  
}



#___