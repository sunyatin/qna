
# Version Control:
# 0     100521   creation
# 1     170521   exclude run if no *.stats file + changed psmc comparison functions + abc.adjust()
# 2     220621   now possibility to pass multiple input files as `target`in abc.import() +
#                changed importation of PSMC (since sumstats.py v12.0, the names are given as *.LABEL.psmc2.stats instead of *.INDEX.psmc2.stats)
#                changed importation of cov stats (since sumstats.py v12.0, names as *.AscertainmentMode.cov.stats instead of *.cov.stats)
# 3     281121   new option: read.median.tract.length=FALSE
# 4     040522   added import for f4-ratios
# 5     270622   added option to abc.import to exclude .cov.stats results if IC95.low < 0 (not by default)
# 6     010722   now imports empirical IQR95 & I removed the length.SE, match.SE, alpha.SE whose calculation made little sense
# 7     191222   switched args order `MAPE <- function(obs, sim)` to `MAPE <- function(sim, obs)` for homogeneity

#################################################################################
#################################################################################
#############################            Core           #########################
#################################################################################
#################################################################################


nsim <- function(data) cat(length(data$NAME))

extract_best_psmc <- function(psmc_file,
                             mut.rate,
                             generation.time,
                             binsize = 100, # in base pairs
                             select.round = -1,
                             verbose = FALSE,
                             warning = FALSE)
{

  # If select.round == -1, the function will select the last iteration.
  # If select.round == -NA, the function will select the run with best goodness-of-fit (GOF).
  # GOF is measured by looking at the RI field (= Relative Information) which is a K-L distance.
  # For more details, cf Li (2011) SM: Remark 9-10.
  psmc = scan(psmc_file, what = character(), sep = "\n", quiet = TRUE)
  rd = which(psmc=="//")
  DF = RI = c()
  for (i in 1:(length(rd)-1)) {
    theta = psmc[rd[i+1]-1]
    if (!grepl("^PA", theta)) stop("Error1.")
    theta = strsplit(theta, "\t")[[1]][2]
    theta = as.numeric(strsplit(theta, " ")[[1]][2])
    df = psmc[(rd[i]+9):(rd[i+1]-2)]
    if (!all(sapply(df, function(x) grepl("^RS", x)))) stop("Error2.")
    ri = psmc[rd[i]+5]
    if (!grepl("^RI", ri)) stop("Error 3.")
    RI = c(RI, as.numeric(strsplit(ri, "\t")[[1]][2]))
    df = t(sapply(df, function(x) strsplit(x, "\t")[[1]][3:4]))
    DF = rbind(DF, data.frame(round=i, x=as.numeric(df[,1]), y=as.numeric(df[,2]), theta=theta))
  }
  DF$No = DF$theta / (4*mut.rate*binsize)
  DF$time_yBP = 2 * DF$x * DF$No * generation.time
  DF$Ne = DF$y * DF$No
  if (is.na(select.round)) {
    select.round = which.min(RI)
    cat(select.round,"\n")
  }
  if (!is.na(select.round) && select.round==-1) select.round = max(as.numeric(DF$round))
  if (verbose) {
    if (!is.na(select.round)) print(paste0("Selecting run: ",select.round)) else cat("Using GOF to select round: ")
  }
  if (!select.round%in%DF$round) return(list("GOFs"=RI, "PSMC"=data.frame(time_yBP=NA, Ne=NA)))
  DF = subset(DF, round==select.round)
  if (select.round!=which.min(RI)) {
    if (warning) {
      warning(paste0(psmc_file,": selected run ",select.round," but the one with lowest GOF is ",which.min(RI)))
    }
  }
  return( list("GOFs"=RI, "PSMC"=DF[,c("time_yBP", "Ne")]) )
}


abc.get_psmc = function(prefix,
                    index,
                    generation.time,
                    mut.rate = 1.2e-8,
                    run = -1)
{
  file0 <- paste0(prefix,".",index,".psmc2.stats")
  file1 <- paste0(prefix,".",index,".psmc.stats")
  if (file.exists(file0)) {
    dt <- read.table(file0, sep="\t", header=F, comment.char="#")
    dt <- data.frame("time_yBP" = dt[,1] * generation.time, "Ne" = dt[,2])
  } else if (file.exists(file1)) {
    dt <- extract_best_psmc(file1,
                            generation.time,
                            mut.rate = mut.rate,
                            select.round = run)$PSMC
  } else {
    dt <- data.frame("time_yBP"=NA, "Ne"=NA)
  }
  return(dt)
}


abc.import <- function(target,
                   generation.time,
                   read.params = TRUE,
                   Sprime.min.match.rate = 0.0,
                   params.to.exclude = c("seq_length", "n_seq", "mut_rate", "recomb_rate", "generation_time"),
                   ascertainment.mode = 0,
                   do.AFS.2D = FALSE,
                   read.median.tract.length = FALSE,
                   # optional, for old *.psmc.stats format:
                   psmc..mut.rate = 1.2e-8,
                   psmc..run = -1,
                   # to extract specific runs
                   files = NA,
                   #
                   ld.exclude.negative.ic95 = FALSE)
{
  # `target` can be either (i) the location of a directory or (ii) the specific run to analyze

  if (dir.exists(target[1])) {
    cat("Reading files in a directory.\n")
    prefixes <- list.files(target, pattern = ".par", full.names = TRUE)
    prefixes <- prefixes[!grepl("failed", prefixes)]
    prefixes <- gsub(".par$", "", prefixes)
    prefixes <- prefixes[file.exists(paste0(prefixes,".stats"))]
  } else {
    cat("Reading specific file(s).\n")
    # x <- sapply(target, function(x) { if (!file.exists(paste0(x,".stats"))) stop(paste0(x,": *.par or *.stats file do(es) not exist.")) })
    prefixes = target
  }

  if (!is.na(files[1])) {
    files <- files[!duplicated(files)]
    cat("/!\ Subsetting to",length(files),"specific run(s).\n")
    xx <- sapply(prefixes, function(x) rev(strsplit(x, "/", fixed=T)[[1]])[1])
    xx <- sapply(xx, function(x) gsub(".par", "", x))
    files <- as.character(files)
    notin <- files[!files%in%xx]
    if (length(notin)>0) {
      cat("Warning. These files are not in the directory:\n")
      print(notin)
      cat("\n\n")
    }
    prefixes <- prefixes[xx%in%files]
  }

  # initialize some lists and vectors
  sumstats <- c("D", "DCFS", "AFS.CEU", "AFS.YRI", "AFS.2D", "Fst", "PI", "SPRIME", "CRF", "LD", "PSMC.CEU", "PSMC.YRI", "f4")
  X <- as.list(rep(NA, length(sumstats)))
  names(X) <- sumstats
  X <- list("NAME"=NA, "PAR"=NA, "SS"=X)

  TEMPLATE <- as.list(rep(NA, length(prefixes)))
  names(TEMPLATE) <- prefixes

  # IDs
  X$NAME <- prefixes

  # parameters
  if (read.params) {
    cat("> Parameters.\n")
    S <- TEMPLATE
    np <- NA
    for (prefix in prefixes) {
      file <- paste0(prefix,".par")
      if (file.exists(file)) {
        par <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
        par <- subset(par, ! V1 %in% params.to.exclude)
        np <- unique(c(np, nrow(par)))
        S[[prefix]] <- c(par[,2])
        names(S[[prefix]]) <- par[,1]
      }
    }
    if (length(na.omit(np))!=1) stop("Cannot merge parameters: there is an unequal number of parameters between runs.")
    X$PAR <- do.call(rbind, S)
    rownames(X$PAR) <- NULL
    if (all(sapply(paste0(prefixes,".par"), file.exists)==FALSE)) X$PAR <- NA
  }
  if (all(is.na(c(X$PAR)))) X$PAR <- NA

  # stats
  cat("> Statistics.\n")
  stats <- c("AFS.CEU", "AFS.YRI", "AFS.2D", "DCFS", "D", "Fst", "PI", "f4")
  S <- as.list(rep(NA, length(stats))); names(S) <- stats
  S <- lapply(prefixes, function(x) S); names(S) <- prefixes
  for (prefix in prefixes) {
    file <- paste0(prefix,".stats")
    if (file.exists(file)) {
      dt <- scan(file, what = "chr", sep = "\n", quiet = TRUE)
      for (x in dt) {
        x <- strsplit(x, " ")[[1]]
        if (x[1]=="AFS_CEU") S[[prefix]]$AFS.CEU <- as.numeric(x[3:length(x)])
        if (x[1]=="AFS_YRI") S[[prefix]]$AFS.YRI <- as.numeric(x[3:length(x)])
        if (do.AFS.2D==TRUE && x[1]=="2D_AFS") S[[prefix]]$AFS.2D <- as.numeric(x[3:length(x)])
        if (x[1]=="DCFS") S[[prefix]]$DCFS <- as.numeric(x[3:length(x)])
        if (grepl("^D\\_", x[1])) { S[[prefix]]$D <- as.numeric(x[c(3,4,5)]); names(S[[prefix]]$D) <- c("D", "SE", "Z") }
        if (x[1]=="Fst") S[[prefix]]$Fst <- as.numeric(x[3])
        if (x[1]=="pi_CEU_per_kb") PI.CEU <- as.numeric(x[3])
        if (x[1]=="pi_YRI_per_kb") PI.YRI <- as.numeric(x[3])
        if (grepl("^F4ratio\\_", x[1])) { S[[prefix]]$f4 <- as.numeric(x[c(3,4,5)]); names(S[[prefix]]$f4) <- c("f4", "SE", "Z") }
      }
      PI <- c(PI.CEU, PI.YRI); names(PI) <- c("CEU", "YRI"); S[[prefix]]$PI <- PI
    }
  }
  for (stat in stats) {
    SS <- lapply(S, "[[", stat)
    X$SS[[stat]] <- do.call(rbind, SS)
    rownames(X$SS[[stat]]) <- NULL
  }

  # LD
  cat("> LD.\n")
  S <- TEMPLATE
  for (prefix in prefixes) {
    file <- ifelse(file.exists(paste0(prefix,".",ascertainment.mode,".cov.stats")),
                paste0(prefix,".",ascertainment.mode,".cov.stats"),
                paste0(prefix,".cov.stats"))
    if (file.exists(file)) {
      dt <- read.table(file, header = FALSE, nrows = 10)
      if (nrow(dt)==0) { S[[prefix]] <- c("mean"=NA, "SE"=NA); next() }
      rownames(dt) <- dt[,1]
      add <- c("mean"=dt["Time_gBSample.jk.MEAN",2], "SE"=dt["Time_gBSample.jk.SE",2])
      ic95 <- dt["Time_gBSample.jk.IC95.low",2]
      add <- if (!is.na(add[1]) && any(add<0)) add * NA else add * generation.time
      if (ld.exclude.negative.ic95) {
        if (!is.na(add[1])) {
          if (!is.na(ic95)) {
            if (ic95<0) add <- add * NA
          }
        }
      }
      S[[prefix]] <- add
    }
  }
  X$SS$LD <- do.call(rbind, S)
  rownames(X$SS$LD) <- NULL
  if (all(sapply(paste0(prefixes,".0.cov.stats"), file.exists)==FALSE) &&
      all(sapply(paste0(prefixes,".1.cov.stats"), file.exists)==FALSE) &&
      all(sapply(paste0(prefixes,".2.cov.stats"), file.exists)==FALSE) &&
      all(sapply(paste0(prefixes,".cov.stats"), file.exists)==FALSE)) X$SS$LD <- NA

  # SPRIME
  cat("> Sprime.\n")
  S <- TEMPLATE
  for (prefix in prefixes) {
    file0 <- paste0(prefix,".sprime.stats")
    file1 <- paste0(prefix,".sprime2.stats")
    if (file.exists(file1)) {
      dt <- read.table(file1, header = TRUE)
      if (nrow(dt)==0) { S[[prefix]] <- c("length.mean" = NA, #"length.SE" = NA,
                                          "match.mean" = NA, #"match.SE" = NA,
                                          "length.2.5p" = NA, "length.97.5p" = NA,
                                          "match.2.5p" = NA, "match.97.5p" = NA); next() }
      if (!Sprime.min.match.rate %in% dt[,1]) stop("Sprime.min.match.rate not in *.sprime.stats")
      dt <- subset(dt, MinFragmentMatchRate == Sprime.min.match.rate)
      Col <- ifelse(read.median.tract.length, "LengthBpMedian", "LengthBpMean")
      S[[prefix]] <- c("length.mean" = dt[,Col],
                       #"length.SE"   = 1/1.96*(dt$LengthBpMean-dt$LengthBpQ2.5),
                       "match.mean"  = dt$MatchRateMean,
                       #"match.SE"    = 1/1.96*(dt$MatchRateMean-dt$MatchRateQ2.5),
                       "length.2.5p"  = dt$LengthBpQ2.5,
                       "length.97.5p" = dt$LengthBpQ97.5,
                       "match.2.5p"   = dt$MatchRateQ2.5,
                       "match.97.5p"  = dt$MatchRateQ97.5)
    } else if (file.exists(file0)) {
      dt <- read.table(file0, header = TRUE)
      if (nrow(dt)==0) { S[[prefix]] <- c("length.mean" = NA, #"length.SE" = NA,
                                          "match.mean" = NA, #"match.SE" = NA,
                                          "length.2.5p" = NA, "length.97.5p" = NA,
                                          "match.2.5p" = NA, "match.97.5p" = NA); next() }
      Col <- ifelse(read.median.tract.length, "Median", "Mean")
      S[[prefix]] <- c("length.median" = dt[1,Col],
                       #"length.SE"    = 1/1.96*(dt$Mean[1]-dt$Q2.5p[1]),
                       "match.mean"   = dt$Mean[2],
                       #"match.SE"     = 1/1.96*(dt$Mean[2]-dt$Q2.5p[2]),
                       "length.2.5p"  = dt$Q2.5p[1],
                       "length.97.5p" = dt$Q97.5p[1],
                       "match.2.5p"   = dt$Q2.5p[2],
                       "match.97.5p"  = dt$Q97.5p[2])
    }
  }
  X$SS$SPRIME <- do.call(rbind, S)
  rownames(X$SS$SPRIME) <- NULL
  if (all(sapply(paste0(prefixes,".sprime.stats"), file.exists)==FALSE) &&
      all(sapply(paste0(prefixes,".sprime2.stats"), file.exists)==FALSE)) X$SS$SPRIME <- NA

  # CRF
  cat("> CRF.\n")
  S <- TEMPLATE
  for (prefix in prefixes) {
    file <- paste0(prefix,".crf.stats")
    if (file.exists(file)) {
      dt <- read.table(file, header = TRUE)
      if (nrow(dt)==0) { S[[prefix]] <- c( "length.mean" = NA, #"length.SE" = NA,
                                           "alpha.mean" = NA, #"alpha.SE" = NA,
                                           "length.2.5p" = NA, "length.97.5p" = NA,
                                           "alpha.2.5p" = NA, "alpha.97.5p" = NA); next() }
      Col <- ifelse(read.median.tract.length, "Median", "Mean")
      S[[prefix]] <- c( "length.mean" = dt[1,Col],
                        #"length.SE"   = 1/1.96*(dt$Mean[1]-dt$Q2.5p[1]),
                        "alpha.mean"  = dt$Mean[2],
                        #"alpha.SE"    = 1/1.96*(dt$Mean[2]-dt$Q2.5p[2])
                        "length.2.5p"  = dt$Q2.5p[1],
                        "length.97.5p" = dt$Q97.5p[1],
                        "alpha.2.5p"   = dt$Q2.5p[2],
                        "alpha.97.5p"  = dt$Q97.5p[2]
                        )
    }
  }
  X$SS$CRF <- do.call(rbind, S)
  rownames(X$SS$CRF) <- NULL
  if (all(sapply(paste0(prefixes,".crf.stats"), file.exists)==FALSE)) X$SS$CRF <- NA

  # PSMC
  cat("> PSMC.\n")
  S1 <- S2 <- TEMPLATE
  for (prefix in prefixes) {
    ix1 <- ifelse(file.exists(paste0(prefix,".CEU.psmc2.stats")), "CEU", "10")
    ix2 <- ifelse(file.exists(paste0(prefix,".YRI.psmc2.stats")), "YRI", "0")
    S1[[prefix]] <- as.data.frame(t(abc.get_psmc(prefix, index = ix1, generation.time, psmc..mut.rate, psmc..run)))
    S2[[prefix]]<- as.data.frame(t(abc.get_psmc(prefix, index = ix2, generation.time, psmc..mut.rate, psmc..run)))
    names(S1[[prefix]]) <- paste0("P.",1:ncol(S1[[prefix]]))
    names(S2[[prefix]]) <- paste0("P.",1:ncol(S2[[prefix]]))
  }
  S1 <- as.matrix(data.table::rbindlist(S1, use.names = TRUE, fill = TRUE)); dimnames(S1) <- NULL
  S2 <- as.matrix(data.table::rbindlist(S2, use.names = TRUE, fill = TRUE)); dimnames(S2) <- NULL
  X$SS$PSMC.CEU <- if (all(c(is.na(S1)))) NA else list("x.yBP"=S1[seq(1, nrow(S1), by=2),,drop=F], "y"=S1[seq(2, nrow(S1), by=2),,drop=F])
  X$SS$PSMC.YRI <- if (all(c(is.na(S2)))) NA else list("x.yBP"=S2[seq(1, nrow(S2), by=2),,drop=F], "y"=S2[seq(2, nrow(S2), by=2),,drop=F])
  #if (all(sapply(paste0(prefixes,".CEU.psmc2.stats"), file.exists)==FALSE)) X$SS$PSMC.CEU <- NA
  #if (all(sapply(paste0(prefixes,".YRI.psmc2.stats"), file.exists)==FALSE)) X$SS$PSMC.YRI <- NA

  # final check
  ns <- length(X$NAME)
  cat("Number of simulations:",ns,"\n")
  ck <- names(X$SS)[!grepl("PSMC",names(X$SS))]
  if (!all(unlist(sapply(X$SS[ck], function(x) nrow(x)))==ns)) stop("Unequal number of simulations.")
  if (length(X$SS$PSMC.CEU)>1) if (nrow(X$SS$PSMC.CEU$x.yBP)!=ns || nrow(X$SS$PSMC.CEU$y)!=ns) stop("PSMC.CEU: unequal number of simulations.")
  if (length(X$SS$PSMC.YRI)>1) if (nrow(X$SS$PSMC.YRI$x.yBP)!=ns || nrow(X$SS$PSMC.YRI$y)!=ns) stop("PSMC.YRI: unequal number of simulations.")
  if (length(X$PAR)>1) if (nrow(X$PAR)!=ns) stop("Unequal number of simulations.")
  #X$SS[ck] <- lapply(X$SS[ck], function(x) if (length(x)==1 && is.na(x)) return(NA) else if (all(is.na(c(x[,1])))) return(NA) else return(x))
  X$SS[ck] <- lapply(X$SS[ck], function(x) if (length(x)==1 && is.na(x)) return(NA) else return(x))
  cat("Check: ok.\n")

  return(X)
}


is_empty <- function(obs.data,
                     sim.data,
                     stat,
                     control.obs = TRUE)
{
  if (control.obs) if (!stat %in% names(obs.data$SS)) return(TRUE)
  if (!stat %in% names(sim.data$SS)) return(TRUE)
  if (length(sim.data$SS[[stat]])==1 && is.na(sim.data$SS[[stat]])) return(TRUE)
  if (control.obs) if (length(obs.data$SS[[stat]])==1 && is.na(obs.data$SS[[stat]])) return(TRUE)
  return(FALSE)
}


abc.analyze <- function(obs.data,
                        sim.data,
                        parfile,
                        N = NA,
                        dispersion.fun = "sd",
                        psmc..xlim = c(NA, NA),
                        verbose = TRUE)
{
  if (verbose) cat("/!\\ PSMC statistics not reported in `all.sim`\n")
  if (verbose) cat("If you want to compute the distance on subvalues of a statistic that are on different scales,\nyou must provide these sublevels separately in the parfile, e.g. D:D and D:Z.\n\n")

  parfile <- read.table(parfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  rownames(parfile) <- parfile[,1]

  raw.d <- matrix(NA, nrow = length(sim.data$NAME), ncol = nrow(parfile), dimnames = list(NULL, parfile[,1]))
  all.sim <- all.obs <- c()
  for (i in 1:nrow(parfile)) {
    stat <- parfile[i,1]
    distfun <- parfile[i,2]
    wg <- parfile[i,3]
    sn <- strsplit(stat, ":", fixed = TRUE)[[1]]
    if (is_empty(obs.data, sim.data, sn[1])) stop("The statistic ",stat," is empty.")

    # case where PSMC
    if (grepl("PSMC", sn[1])) {
      if (verbose) cat(stat,"\n")
      r <- sapply(1:nrow(sim.data$SS[[sn[1]]]$x.yBP), function(j) {
        do.call(distfun, list(obs = cbind(obs.data$SS[[sn[1]]]$x.yBP[1,], obs.data$SS[[sn[1]]]$y[1,]),
                              sim = cbind(sim.data$SS[[sn[1]]]$x.yBP[j,], sim.data$SS[[sn[1]]]$y[j,]),
                              xlim = psmc..xlim))
      })
      raw.d[,stat] <- r
      next()
    }

    # case where we take all values of the statistic
    if (length(sn)==1) {
      if (verbose) cat(sn[1],":",ncol(obs.data$SS[[sn[1]]]),"col(s)\n")
      s <- sim.data$SS[[sn[1]]]
      colnames(s) <- if (is.null(colnames(s))) paste0(sn[1],":",1:ncol(s)) else paste0(sn[1],":",colnames(s))
      o <- obs.data$SS[[sn[1]]]
      colnames(o) <- if (is.null(colnames(o))) paste0(sn[1],":",1:ncol(o)) else paste0(sn[1],":",colnames(o))
      r <- do.call(distfun, list(s, o))
      all.sim <- cbind(all.sim, s)
      all.obs <- cbind(all.obs, o)

    # case where we take one specific value of the statistic
    } else if (length(sn)==2) {
      if (verbose) cat("____________________",sn[1],":",sn[2],"\n")
      if (grepl("^\\d", sn[2])) {
        c <- as.integer(sn[2])
        if (! c %in% 1:ncol(sim.data$SS[[sn[1]]])) stop(c," not in the space of ",sn[1])
      } else {
        c <- sn[2]
        if (! c %in% colnames(sim.data$SS[[sn[1]]])) stop(c," not in the colnames of ",sn[1])
      }
      s <- sim.data$SS[[sn[1]]][,c,drop=F]
      colnames(s) <- paste0(sn[1],":",sn[2])
      o <- obs.data$SS[[sn[1]]][,c,drop=F]
      colnames(o) <- paste0(sn[1],":",sn[2])
      r <- do.call(distfun, list(s, o))
      all.sim <- cbind(all.sim, s)
      all.obs <- cbind(all.obs, o)
    } else {
      stop("Error.")
    }
    raw.d[,stat] <- r
  }

  # remove duplicated columns in the all.sim/all.obs table
  dup <- duplicated(colnames(all.obs))
  if (sum(dup)>1) warning(paste0("Found ",sum(dup)," duplicated colnames"))
  all.sim <- all.sim[,!dup,drop=F]
  all.obs <- all.obs[,!dup,drop=F]

  # normalize distances
  means <- colMeans(raw.d, na.rm = TRUE)
  if (dispersion.fun=="sd") {
    fun <- sd
  } else if (dispersion.fun=="mad") {
    fun <- mad
  } else {
    stop("This dispersion.fun is not implemented.")
  }
  sds <- apply(raw.d, 2, fun, na.rm = TRUE)
  if (any(sds==0)) {
    print(colnames(raw.d)[sds==0])
    print(raw.d[,sds==0])
    stop("Problem! There are some dispersion values that are null. If you used the MAD, try SD instead.")
  }
  norm.d <- sweep(raw.d, 2, means, FUN="-")
  norm.d <- sweep(norm.d, 2, sds, FUN="/")
  norm.obs <- ( rep(0, length(means)) - means ) / sds
  norm.d <- abs(sweep(norm.d, 2, norm.obs, FUN="-"))

  # weighted normalized distances
  if (!all(colnames(norm.d) %in% rownames(parfile))) stop("Some colnames of norm.d not in the parameter-weight file!")
  wg <- parfile[colnames(norm.d),3]
  wgt.d <- sweep(norm.d, 2, wg, FUN="*")

  # total distances
  tot.d <- rowSums(wgt.d)

  # best runs
  accepted <- data.frame("index"=seq_along(sim.data$NAME), "name"=sim.data$NAME, "total.distance"=tot.d, stringsAsFactors = FALSE)
  accepted <- accepted[order(accepted$total.distance),]
  if (is.na(N)) N <- nrow(accepted)
  accepted <- accepted[1:N,]
  rownames(accepted) <- accepted$index

  names(parfile) <- c("sumstats", "distfun", "weight")
  r <- list("ACCEPTED" = accepted,
       "NAME" = sim.data$NAME,
       "DISTANCE" = tot.d,
       "DISTANCES" = list("raw"=raw.d, "normalized"=norm.d, "normalized.weighted"=wgt.d),
       "sumstats.merged.sim" = all.sim,
       "sumstats.merged.obs" = all.obs,
       "ARGUMENTS" = parfile)

  if (verbose) cat("\n")
  return(r)
}


getmode <- function(x, na.rm = TRUE)
{
  density_estimate <- density(x, na.rm = na.rm)
  return(density_estimate$x[which.max(density_estimate$y)])
}


abc.adjust <- function(param, distance, size = 1000L, seed = 42)
{
  if (F) {
    # add the null
    x <- c(0, distance)
    # normalization
    x <- -1 * (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    # scale to 0-1 to construct the weights
    y <- (x-min(x))/(max(x)-min(x))
    y <- y[-1]
  } else {
    win <- max(distance) + sd(c(0, distance))
    y <- stats::dnorm(distance, sd = win/2)
  }
  par <- apply(param, 2, function(x) {
    set.seed(seed)
    sample(x, size = size, replace = TRUE, prob = y)
  })
  return(par)
}

# extract simulations based on their indices
extract <- function(data, indices)
{
  if (length(data$PAR)==1 && is.na(data$PAR)) {
		data$PAR <- NA
	} else {
		data$PAR <- data$PAR[indices,,drop=F]
	}
  data$NAME <- data$NAME[indices]
  data$SS <- lapply(data$SS, function(x) {
    if (length(x)==1 && is.na(x)) return(NA)
    if (is.list(x)) return(lapply(x, function(xx) xx[indices,,drop=F]))
    return(x[indices,,drop=F])
  })
  return(data)
}

# performs goodness-of-fit tests
assess_fit <- function(obs.data, sim.data, file.of.stats, output.png, nb.replicates = 1000L, tol = 0.01)
{
  # `file.of.stats` must have format with main.stat[:sub.stat] on each line
  require(abc)
  warning("Currently, assess_fit is implemented with no weighting for the vectorial statistics (eg. AFS, DCFS)!")
  ST <- unlist(read.table(file.of.stats, header = FALSE, stringsAsFactors = FALSE))
  RES <- data.frame(statistic = ST, P = NA)
  for (i in seq_along(ST)) {
    st <- ST[i]
    cat(st," ")
    if (grepl(":", st, fixed=T)) {
      st <- strsplit(st, ":")[[1]]
      O <- as.matrix(obs.data$SS[[st[1]]][,st[2],drop=F])
      S <- as.matrix(sim.data$SS[[st[1]]][,st[2],drop=F])
    } else {
      O <- as.matrix(obs.data$SS[[st[1]]])
      S <- as.matrix(sim.data$SS[[st[1]]])
    }
    gf <- suppressWarnings(abc::gfit(target = O,
                                     sumstat = S,
                                     nb.replicate = nb.replicates,
                                     tol = tol,
                                     statistic = median))
    RES[i,2] <- summary(gf)$pvalue
  }
  RES$statistic <- factor(RES$statistic, level = RES$statistic[order(RES$P)])
  RES <- RES[order(RES$P),]
  gf <- ggplot(GF) +
    geom_point(aes(x=statistic, y=P, col=P>0.05), size = 4) +
    geom_segment(aes(x=statistic, xend=statistic, y=0., yend=P, col=P>0.05), size = 1.) +
    coord_flip() +
		ggdist::theme_ggdist() +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    scale_color_viridis(discrete = T, begin = .2, end = .8, direction = -1) +
    ylab("P-value") + xlab("") +
    ggtitle(paste0("Goodness-of-fit test\nNum replicates: ",nb.replicates," - Tolerance: ",tol)); gf
  ggsave(gf, filename = output.png)
  return(RES)
}

#################################################################################
#################################################################################
#############################         Distances         #########################
#################################################################################
#################################################################################


psmc.align <- function(obs, sim, xlim, psmc.approx.method = "constant")
{
  if (!is.na(xlim[1])) obs <- na.omit(obs[obs[,1]>=xlim[1] & obs[,1]<=xlim[2],])
    if (all(is.na(sim[,2]))) {
    return(list("obs.y"=obs[,2], "sim.y"=rep(NA, nrow(sim))))
  } else {
    y <- approx(x = sim[,1], y = sim[,2], xout = obs[,1], method = psmc.approx.method, rule = 2, f = 0, ties = "ordered")$y
    return(list("obs.y"=obs[,2], "sim.y"=y, "x"=obs[,1]))
  }
}

psmc.CorDistance <- function(obs, sim, xlim)
{
  x <- psmc.align(obs, sim, xlim)
  if (all(is.na(x$sim.y)) | var(x$sim.y)==0) return(NA)
  return(TSdist::CorDistance(x$obs.y, x$sim.y))
}

psmc.CrossCorDistance <- function(obs, sim, xlim)
{
  # Note. Not a proper distance metric because
  # it does not verify the triangular inequality.
  lag.max <- NULL
  x <- psmc.align(obs, sim, xlim)
  if (all(is.na(x$sim.y)) | var(x$sim.y)==0) return(NA)
  ccr <- ccf(x$obs.y, x$sim.y, plot = FALSE, type = c("correlation"), lag.max = lag.max)
  ccr <- max(abs(ccr$acf))
  return(1-ccr)
}

###

l2 <- function(sim, obs)
{
  if (nrow(obs)!=1) stop("The observed data should contain a single dataset.")
  y <- sweep(x=sim, 2, obs, FUN="-")**2
  y <- sqrt(rowSums(y))
  return(y)
}

MRAE <- function(sim, obs) # Mean Relative Absolute Error
{
  if (nrow(obs)!=1) stop("The observed data should contain a single dataset.")
  y <- sweep(x=sim, 2, obs, FUN="-")
  y <- sweep(y, 2, obs, FUN="/")
	y <- abs(y)
  y <- rowMeans(y)
  return(y)
}

RMSE <- function(sim, obs) # Mean Relative Absolute Error
{
  if (nrow(obs)!=1) stop("The observed data should contain a single dataset.")
  y <- sweep(x=sim, 2, obs, FUN="-")**2
  y <- sqrt(rowMeans(y))
  return(y)
}

#___
