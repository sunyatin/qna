
abc.par <- function(sim.data, 
                    analysis.results, 
                    priors, 
                    output, 
                    nbins = 20, 
                    size = 1000L, 
                    true.values = NA, 
                    do.adjusted = FALSE, 
                    ks.threshold = 0.1)
{
  
  graphics.off()
  
  R <- analysis.results
  accepted <- R$ACCEPTED
  d <- matrix(R$DISTANCE[accepted[,1]], ncol = 1); colnames(d) = "distance"
  
  priors <- priors[!duplicated(priors[,1]),,drop=F]
  rownames(priors) <- as.character(priors[,1])
  priors <- subset(priors, grepl("unif", V2))
  PAR <- sim.data$PAR[accepted[,1],,drop=FALSE]
  ix <- intersect(as.character(priors[,1]), colnames(PAR))
  PAR <- PAR[,ix,drop=F]
  priors <- priors[ix,,drop=F]
  
  # adjusted parameters
  if (do.adjusted) {
    invisible(capture.output(ADJ <- suppressWarnings(abc::abc(target = 0.0,
                                                              param = PAR,
                                                              sumstat = d,
                                                              tol = 1.0,
                                                              method = "neuralnet",
                                                              numnet = 100,
                                                              sizenet = 3,
                                                              hcorr = TRUE,
                                                              transf = "logit",
                                                              logit.bounds = priors[,3:4],
                                                              kernel = "epanechnikov"))))
    ADJ <- ADJ$adj.values
    
    # RT adjustment function
    #hist(c(d), main = "", xlab = "Distances", col = "lightgray")
    RT <- abc.adjust(PAR, c(d), size = size, seed = 42)
  }
  
  SUMMARY <- SUMMARY.2 <- SUMMARY.3 <- SD <- c()
  for (p in 1:nrow(priors)) {
    pname <- as.character(priors[p,1])
    if (grepl("log10unif", priors[p,2])) {
      PAR[,pname] <- log10(PAR[,pname])
      if (do.adjusted) ADJ[,pname] <- log10(ADJ[,pname])
      if (do.adjusted) RT[,pname] <- log10(RT[,pname])
      priors[p,3:4] <- log10(priors[p,3:4])
    }
    SUMMARY <- rbind(SUMMARY, data.frame(V1=pname, mode=getmode(PAR[,pname]), mean=mean(PAR[,pname])))
    if (do.adjusted) SUMMARY.2 <- rbind(SUMMARY.2, data.frame(V1=pname, mode=getmode(ADJ[,pname]), mean=mean(ADJ[,pname])))
    if (do.adjusted) SUMMARY.3 <- rbind(SUMMARY.3, data.frame(V1=pname, mode=getmode(RT[,pname]), mean=mean(RT[,pname])))
    SD <- rbind(SD, data.frame(V1=pname, sd=sd(PAR[,pname])))
  }
  PAR.DF <- melt(PAR)
  PAR.DF <- data.frame(V1=PAR.DF$Var2, V2=PAR.DF$value)
  if (do.adjusted) ADJ.DF <- melt(ADJ)
  if (do.adjusted) ADJ.DF <- data.frame(V1=ADJ.DF$Var2, V2=ADJ.DF$value)
  if (do.adjusted) RT.DF <- melt(RT)
  if (do.adjusted) RT.DF <- data.frame(V1=RT.DF$Var2, V2=RT.DF$value)
  
  # 20.12.21 => KS test of uniformity
  #print(priors)
  #KS <- c()
  #CHI <- c()
  AllTests <- c()
  for (i in 1:ncol(PAR)) {
    par <- colnames(PAR)[i]
    #if (priors[priors$V1==par,"V2"]=="unif") transf <- function(x) return(x) else transf <- function(x) return(log10(x))
    # things are already transformed!
    transf <- function(x) return(x)
    set.seed(42)
    Unif <- runif(n = 1e3L, transf(priors[priors$V1==par,"V3"]), transf(priors[priors$V1==par,"V4"]))
    ks <- ks.test(Unif, transf(PAR[,i]), exact = F)$p.value
    #KS <- rbind(KS, data.frame(V1 = par, ks = ks))
    if (length(PAR[,i]) >= 30) {
      chi <- spgs::chisq.unif.test(transf(PAR[,i]),
                                   interval = transf(c(priors[priors$V1==par,"V3"], priors[priors$V1==par,"V4"])),
                                   min.bin.size = 5)$p.value
    } else {
      chi <- NA
    }
    #CHI <- rbind(KS, data.frame(V1 = par, chi = chi))
    AllTests <- rbind(AllTests, data.frame(V1 = par, ks = ks, chi = chi))
  }
  AllTests$ks <- sapply(AllTests$ks, function(ks) as.character(round(ks, 3)))
  AllTests$chi <- sapply(AllTests$chi, function(chi) as.character(round(chi, 3)))
  rownames(AllTests) <- AllTests$V1
  AllTests$text <- sapply(1:nrow(AllTests), function(ii) paste0(AllTests[ii,2],"\n",AllTests[ii,3]))
  #PAR.DF$tests <- sapply(PAR.DF$V1, function(ii) paste0(AllTests[ii,2],"\n",AllTests[ii,3]))
  PAR.DF$NotUniform <- sapply(PAR.DF$V1, function(ii) {
    if (is.na(AllTests[ii,3])) return( AllTests[ii,2]<ks.threshold )
    return( AllTests[ii,2]<ks.threshold & AllTests[ii,3]<ks.threshold )
  }  )
  
  g1 <- ggplot(priors) +
    geom_rect(aes(xmin=V3, xmax=V4, ymin=0, ymax=1, col=V2=="unif"), alpha=.9, size = 1.2) +
    geom_histogram(data=PAR.DF, aes(x=V2, fill=NotUniform), bins=nbins, alpha=.7) +
    facet_wrap(~V1, scales="free", ncol=6) +
    geom_vline(data = SUMMARY, aes(xintercept = mode)) +
    geom_vline(data = SUMMARY, aes(xintercept = mean), linetype = "dotted", size = 1.) +
    theme_classic() + xlab("") + ylab("") +
    theme(legend.position = "none") +
    #scale_fill_manual(values = c("lightgray", "dodgerblue4")) +
    ggtitle("Unadjusted values") +
    geom_text(data = AllTests, x=-Inf, y=Inf, aes(label = text), hjust = 0, vjust = 1)
  
  if (do.adjusted) g2 <- ggplot(priors) +
    geom_rect(aes(xmin=V3, xmax=V4, ymin=0, ymax=1, fill=V2=="unif"), alpha=.9) +
    geom_histogram(data=ADJ.DF, aes(x=V2), fill="red", bins=nbins, alpha=.6) +
    facet_wrap(~V1, scales="free", ncol=6) +
    geom_vline(data = SUMMARY.2, aes(xintercept = mode)) +
    geom_vline(data = SUMMARY.2, aes(xintercept = mean), linetype = "dotted", size = 1.) +
    theme_classic() + xlab("") + ylab("") +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("lightgray", "dodgerblue4")) +
    ggtitle("Neural net adjusted values")
  
  if (do.adjusted) g3 <- ggplot(priors) +
    geom_rect(aes(xmin=V3, xmax=V4, ymin=0, ymax=size*0.1, fill=V2=="unif"), alpha=.9) +
    geom_histogram(data=RT.DF, aes(x=V2), fill="red", bins=nbins, alpha=.6) +
    facet_wrap(~V1, scales="free", ncol=6) +
    geom_vline(data = SUMMARY.3, aes(xintercept = mode)) +
    geom_vline(data = SUMMARY.3, aes(xintercept = mean), linetype = "dotted", size = 1.) +
    theme_classic() + xlab("") + ylab("") +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("lightgray", "dodgerblue4")) +
    ggtitle("Weighted values")
  
  if (!is.na(true.values[1])) {
    Par <- data.frame(V1=colnames(true.values), x=c(true.values), stringsAsFactors = FALSE)
    Par <- subset(Par, V1 %in% priors$V1)
    Priors <- priors; rownames(Priors) <- Priors[,1]
    Par[Priors[Par$V1,2]=="log10unif","x"] <- log10(Par[Priors[Par$V1,2]=="log10unif","x"])
    g1 <- g1 + geom_vline(data=Par, aes(xintercept = x), col = "black", size = 1.5)
    if (do.adjusted) g2 <- g2 + geom_vline(data=Par, aes(xintercept = x), col = "black", size = 1.5)
    if (do.adjusted) g3 <- g3 + geom_vline(data=Par, aes(xintercept = x), col = "black", size = 1.5)
  }
  
  invisible(pdf(output, width = 20, height = 12))
  invisible(print(g1))
  if (do.adjusted) invisible(print(g2))
  if (do.adjusted) invisible(print(g3))
  dev.off()
  
  for (p in 1:nrow(SUMMARY)) {
    pname <- as.character(SUMMARY[p,1])
    if (grepl("log10unif", priors[pname,2])) {
      SUMMARY[p,2:3] <- sapply(SUMMARY[p,2:3], function(x) 10**x)
      if (do.adjusted) SUMMARY.2[p,2:3] <- sapply(SUMMARY.2[p,2:3], function(x) 10**x)
      if (do.adjusted) SUMMARY.3[p,2:3] <- sapply(SUMMARY.3[p,2:3], function(x) 10**x)
      SD[p,2:3] <- sapply(SD[p,2], function(x) 10**x)
    }
  }
  
  if (do.adjusted) {
    x <- data.frame(par=SUMMARY$V1, unadj.mode=SUMMARY$mode, unadj.mean=SUMMARY$mean,
                    adj.mode=SUMMARY.2$mode, adj.mean=SUMMARY.2$mean,
                    rt.mode=SUMMARY.3$mode, rt.mean=SUMMARY.3$mean,
                    unadj.sd=SD$sd)
    write.table(x, paste0(output,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    #
    cat("\n\nNumber of simulations:\n", file = paste0(output,".txt"), append = TRUE)
    cat(nrow(R$ACCEPTED), file = paste0(output,".txt"), append = TRUE)
    #
    cat("\n\nAccepted (ordered, closest 1st):\n", file = paste0(output,".txt"), append = TRUE)
    cat(paste(R$ACCEPTED[,2], collapse=" "), file = paste0(output,".txt"), append = TRUE)
    #
    cat("\n\nAccepted (simplified & ordered, closest 1st):\n", file = paste0(output,".txt"), append = TRUE)
    cat(paste(sapply(R$ACCEPTED[,2], function(x) rev(strsplit(x, "/", fixed=T)[[1]])[1]), collapse=" "),
        file = paste0(output,".txt"), append = TRUE)
    #
    cat("\n\nDistances:\n", file = paste0(output,".txt"), append = TRUE)
    cat(paste(R$ACCEPTED[,3], collapse=" "), file = paste0(output,".txt"), append = TRUE)
    #
    cat("\n\nArguments:\n", file = paste0(output,".txt"), append = TRUE)
    write.table(R$ARGUMENTS, paste0(output,".txt"), quote = FALSE, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  
  #ggsave(g, filename = output, width=20, height=12)
}


get_weighted_mode <- function(x, wgt, na.rm = TRUE)
{
  if (na.rm) xw <- na.omit(cbind(x, wgt))
  xw[,2] <- xw[,2] / sum(xw[,2], na.rm = T)
  density_estimate <- density(xw[,1], weights = xw[,2], na.rm = F)
  return(density_estimate$x[which.max(density_estimate$y)])
}


abc.par.2 <- function(sim.data, 
                      analysis.results, 
                      output, 
                      priors, 
                      directory.priors = NA, 
                      size = 1000L, 
                      nnet = 200L)
{
  
  graphics.off()
  if (!is.na(directory.priors)) { if (!dir.exists(directory.priors)) { directory.priors <- NA } }
  
  R <- analysis.results
  accepted <- R$ACCEPTED
  d <- matrix(R$DISTANCE[accepted[,1]], ncol = 1); colnames(d) = "distance"
  
  # prior distributions
  EST <- priors
  rownames(EST) <- as.character(EST[,1])
  priors <- priors[!duplicated(priors[,1]),,drop=F]
  rownames(priors) <- as.character(priors[,1])
  priors <- subset(priors, grepl("unif", V2))
  PAR <- sim.data$PAR[accepted[,1],,drop=FALSE]
  ix <- intersect(as.character(priors[,1]), colnames(PAR))
  PAR <- PAR[,ix,drop=F]
  priors <- priors[ix,,drop=F]
  
  # effective prior parameters
  if (!is.na(directory.priors)) {
    cat("\n\nImporting prior parameters.\n")
    lp <- list.files(directory.priors, pattern = ".par")
    lp <- lp[sapply(lp, function(x) grepl("\\.par$", x, fixed=F))]
    lp <- lp[!grepl("failed", lp)]
    lp <- lp[!grepl("CRF", lp)]
    PREPARS <- lapply(paste0(directory.priors,"/",lp), function(x) read.table(x, header=F)[,2])
    PREPARS <- do.call(rbind, PREPARS)
    colnames(PREPARS) <- read.table(paste0(directory.priors,"/",lp[1]), header=F, stringsAsFactors=F)[,1]
    PREPARS <- PREPARS[,colnames(PAR)]
  }
  
  # transforming
  cat("Transforming parameters.\n")
  for (p in 1:nrow(priors)) {
    pname <- as.character(priors[p,1])
    if (grepl("log10unif", priors[p,2])) {
      priors[p,3:4] <- log10(priors[p,3:4])
      PAR[,pname] <- log10(PAR[,pname])
      if (!is.na(directory.priors)) PREPARS[,pname] <- log10(PREPARS[,pname])
    }
  }
  
  # adjusted parameters
  invisible(capture.output(ADJ <- suppressWarnings(abc::abc(target = 0.0,
                                                              param = PAR,
                                                              sumstat = d,
                                                              tol = 1.0,
                                                              method = "neuralnet",
                                                              numnet = nnet,
                                                              sizenet = 5,
                                                              hcorr = TRUE,
                                                              transf = "logit",
                                                              logit.bounds = priors[,3:4],
                                                              kernel = "epanechnikov"))))

  WGT <- ADJ$weights
  ADJ <- ADJ$adj.values
  AdjMode <- apply(ADJ, 2, function(x) get_weighted_mode(x, WGT, na.rm=T))
  AdjMode <- t(as.matrix(AdjMode))
  colnames(AdjMode) <- colnames(ADJ)
  
  # RT weights
  RT <- abc.adjust(PAR, c(d), size = size, seed = 42)
  
  cat("Plotting.\n")
  pdf(paste0(output,".pdf"))
  par(mfrow = c(3, 3),
      oma = c(2, 2, 1, 0), 
      mar = c(3, 1, 1, 1))
  for (j in 1:ncol(PAR)) {
    pname <- colnames(PAR)[j]
    ###
    if (!is.na(directory.priors)) {
      ds <- hist(PREPARS[,pname], plot = F)
      ds$counts <- ds$counts / max(ds$counts)
      plot(ds, 
           col = "grey",
           main = pname,
           axes = T,
           xlab = NULL,
           ylab = NULL,
           border = F,
           xlim = c(priors[pname,3], priors[pname,4]),
           ylim = c(0, 1))
    } else {
      plot(NA, 
           col = "grey",
           main = pname,
           axes = T,
           xlab = NULL,
           ylab = NULL,
           xlim = c(priors[pname,3], priors[pname,4]),
           ylim = c(0, 1))
    }
    ### rejection parameters
    ds <- density(PAR[,pname]); ds$y <- ds$y / max(ds$y)
    lines(ds, col = "black")
    axis(side = 2, labels = F)
    ### ABC adjusted parameters
    xw <- na.omit(cbind(ADJ[,pname], WGT))
    xw[,2] <- xw[,2] / sum(xw[,2])
    ds <- density(xw[,1], weights = xw[,2]); ds$y <- ds$y / max(ds$y)
    lines(ds, col = "red")
    ### RT weight-adjusted parameters
    ds <- density(RT[,pname]); ds$y <- ds$y / max(ds$y)
    lines(ds, col = "brown")
  }
  dev.off()
  
  # export the posterior parameter files
  cat("Generating posterior files.\n")
  remplace <- function(neo, EST, output, fun = getmode) {
    est <- EST[,1:2,drop=F]
    EST[,1:2] <- apply(EST[,1:2], 2, as.character)
    for (j in 1:ncol(neo)) {
      pname <- colnames(neo)[j]
      if (!EST[pname,2] %in% c("unif", "log10unif")) next()
      est[pname,2] <- as.character(fun(neo[,pname]))
      if (EST[pname,2] == "log10unif") est[pname,2] <- as.character( 10**as.numeric(as.character(est[pname,2])) )
    }
    write.table(est, output, row.names=F, col.names=F, sep="\t", quote=F)
  }
  remplace(PAR, EST, paste0(output,".rej.par"))
  remplace(AdjMode, EST, paste0(output,".adj.par"), fun = identity)
  remplace(RT, EST, paste0(output,".wgt.par"))
  
}

#___