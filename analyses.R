setwd("~/gitdir/qna")
require(paletteer); require(scico); require(scales); require(ggridges); require(plotrix); require(ggnewscale); require(ggsci)
rm(list=ls(all=TRUE))
source("scripts/abc.R")
source("scripts/plot.R")
source("scripts/abc.par.2.R")
source("scripts/model_comparison.R")

#_____________________________________________________________________________________________________________________#
# fun

do.abc.par <- function(LX, R, output, prior.file, dir.priors) 
{
  priors <- read.table(prior.file, header = FALSE, fill = TRUE, col.names = paste0("V",1:20), stringsAsFactors=F)
  priors[priors$V2=="discrete",2:4] <- data.frame("unif", 
                                                  t(apply(priors[priors$V2=="discrete",4:ncol(priors)], 1, range, na.rm=T)), 
                                                  stringsAsFactors=F)
  priors <- subset(priors, V2!="math")
  priors <- priors[,1:4]
  priors[,3:4] <- apply(priors[,3:4], 2, as.numeric)
  abc.par(LX, R, paste0(output,".par.pdf"), priors = priors, nbins = 15, size = 2000L)
  abc.par.2(LX, R, paste0(output,".par2"), priors = priors, directory.priors = dir.priors, size = 2000L)
}

compare <- function(DirRT, 
                    AccNamesRT, 
                    DirsPublished, 
                    N, 
                    spar, 
                    fun.summary,
                    comparison.spar, 
                    outprefix, 
                    tg, 
                    Sprime.min.match.rate, 
                    psmc..xlim,
                    verbose = FALSE,
                    do.pdf = FALSE ) 
{
  extensions <- c("png")
  if (do.pdf) extensions <- c(extensions, "pdf")
  for (ext in extensions) Res <- model_comparison(Model.Paths = c(list(DirRT), as.list(DirsPublished)),
                            Model.Runs = c(list(AccNamesRT), as.list(rep(NA, length(DirsPublished)))),
                            Model.Names = c(list("z. This study"), as.list(names(DirsPublished))),
                            N = N,
                            distance.metrics.file = spar,
                            distance.metrics.comparison = comparison.spar,
                            do.scaled.distance = TRUE,
                            outpng = paste0(outprefix,".compare.",ext),
                            #
                            generation.time = tg,
                            Sprime.min.match.rate = Sprime.min.match.rate,
                            psmc..xlim = psmc..xlim,
                            #
                            functions.of.summary = fun.summary, 
                            function.of.dispersion = "mad",
                            verbose = verbose,
                            import.AFS.2D = FALSE)
  save(Res, file = paste0(outprefix,".RData"))
  write.table(Res$results, file = paste0(outprefix,".txt"), row.names = F)
  return(Res)
}

plot_set <- function(indir, 
                     outprefix, 
                     spar, 
                     N,
                     do.pdf = FALSE )
{
  cat("tg: ",tg,"\n")
  extensions <- c("png")
  if (do.pdf) extensions <- c(extensions, "pdf")
  DATA <- abc.import(indir, 
                     do.AFS.2D = F, 
                     read.params = F, 
                     generation.time = tg, 
                     Sprime.min.match.rate = Sprime.min.match.rate)
  # Analyze
  R <- abc.analyze(obs.data, 
                   DATA, 
                   parfile = spar, 
                   N = N, 
                   dispersion.fun = "mad", 
                   psmc..xlim = psmc..xlim)
  #---
  R.accepted <- R$ACCEPTED
  #---
  for (ext in extensions) {
    abc.plot(obs.data, 
             DATA, 
             R.accepted, 
             output = paste0(outprefix,".tiledA.",ext), 
             title = spar, 
             line.size = line.size, 
             point.size = point.size, 
             scg = scg,
             sfg = sfg,
             plot.each = T, 
             tile = T, 
             show.LD.obs = F, 
             psmc..ylim = psmc..ylim,
             xlim.D = c(-0.02, .2),
             squeeze = squish,
             layout = "A")
    abc.plot(obs.data, 
             DATA, 
             R.accepted, 
             output = paste0(outprefix,".tiledB.",ext), 
             title = spar, 
             line.size = line.size, 
             point.size = point.size, 
             scg = scg,
             sfg = sfg,
             plot.each = T, 
             tile = T, 
             show.LD.obs = F, 
             psmc..ylim = psmc..ylim,
             xlim.D = c(-0.02, .2),
             squeeze = squish,
             layout = "B")
  }
}

sumup <- function(label, sim, obs, precision = 2, as.percentage = F)
{
  cat("\n\n",label,"\n")
  sim <- na.omit(sim)
  #
  mrae <- MRAE(sim=sim, obs=obs)
  mrae <- 100 * mrae
  cat("MRAE: ",round(mean(mrae),1),"%±",round(sd(mrae),1),"%\n", sep="")
  #
  rmse <- RMSE(sim=sim, obs=obs)
  if (as.percentage) {
    rmse <- 100 * rmse
    cat("RMSE: ",round(mean(rmse),precision),"%±",round(sd(rmse),precision),"%\n", sep="")
  } else {
    cat("RMSE: ",round(mean(rmse),precision),"±",round(sd(rmse),precision),"\n", sep="")
  }
  #
  if (ncol(sim)==1) {
    if (as.percentage) {
      cat("Obs: ",round(100*obs,precision),"%\n", sep="")
    } else {
      cat("Obs: ",round(obs,precision),"\n", sep="")
    }
    #
    if (as.percentage) {
      cat("Sim mean: ",round(100*mean(sim),precision),"%\n", sep="")
      cat("Sim sd: ",round(100*sd(sim),precision),"%\n", sep="")
      cat("Sim se: ",round(100*sd(sim)/sqrt(nrow(sim)),precision),"%\n", sep="")
      cat("Sim IP95: ",round(100*quantile(sim, 2.5/100),precision),"%—",round(100*quantile(sim, 97.5/100),precision),"%\n", sep="")
      cat("Sim MinMax: ",round(100*min(sim),precision),"%—",round(100*max(sim),precision),"%\n", sep="")
    } else {
      cat("Sim mean: ",round(mean(sim),precision),"\n", sep="")
      cat("Sim sd: ",round(sd(sim),precision),"\n", sep="")
      cat("Sim se: ",round(sd(sim)/sqrt(nrow(sim)),precision),"\n", sep="")
      cat("Sim IP95: ",round(quantile(sim, 2.5/100),precision),"—",round(quantile(sim, 97.5/100),precision),"\n", sep="")
      cat("Sim MinMax: ",round(min(sim),precision),"—",round(max(sim),precision),"\n", sep="")
    }
  }
}

multidata_plots <- function(data_paths, 
                            labels,
                            output,
                            title,
                            do_tile,
                            plot.legend = FALSE,
                            color.lims = c(0.1, 0.7), 
                            xylims.sprime = c("xmin"=5e3, "xmax"=750e3, "ymin"=0.1e-2, "ymax"=100e-2), 
                            xylims.crf = c("xmin"=1e4, "xmax"=250e3, "ymin"=0.001, "ymax"=0.03),
                            ylim.afs.ceu = c(0, .32),
                            ylim.afs.yri = c(0, .40),
                            ylim.dcfs = c(0, .40),
                            xlim.D = c(-0.02, .15),
                            squeeze = FALSE,
                            show.RMSE = FALSE,
                            log10.crf = c(TRUE, TRUE),
                            log10.sprime = c(TRUE, TRUE),
                            plot.obs.TLD = FALSE )
{
  print.params <- c("ncol" = 2, "width" = 15, height = 27, scale = .7)
  scg <- scale_color_viridis(discrete = F, option = "C", begin = color.lims[1], end = color.lims[2])
  point.size <- 2.0
  line.size <- 0.7
  names(xylims.sprime) <- c("xmin", "xmax", "ymin", "ymax")
  names(xylims.crf) <- c("xmin", "xmax", "ymin", "ymax")
  #___________________________________________________________________________________#
  if (length(labels) != length(data_paths)) stop("-i must have same length as -l")
  invisible(capture.output(L <- abc.import(data_paths,
                                           read.params = FALSE, 
                                           generation.time = tg,
                                           Sprime.min.match.rate = Sprime.min.match.rate)))
  Accepted <- data.frame(index = seq_along(data_paths), name = data_paths, total.distance = 1)
  Accepted$name <- plyr::mapvalues(Accepted$name, data_paths, labels)
  abc.plot( obs.data, L, Accepted,
            output = output,
            title = title,
            tile = do_tile,
            line.size = line.size,
            point.size = point.size,
            scg = scg,
            n = nrow(Accepted),
            show.LD.obs = FALSE,
            show.legend = plot.legend,
            plot.each = FALSE,
            psmc..ylim = psmc..ylim,
            print.params = print.params,
            plot.title = FALSE,
            log10.crf = log10.crf,
            xylims.crf = xylims.crf,
            log10.sprime = log10.sprime,
            xylims.sprime = xylims.sprime,
            show.RMSE = show.RMSE,
            axis.labs = TRUE,
            ylim.afs.ceu = ylim.afs.ceu,
            ylim.afs.yri = ylim.afs.yri,
            xlim.D = xlim.D,
            ylim.dcfs = ylim.dcfs,
            layout = "C",
            squeeze = squeeze,
            plot.obs.TLD = plot.obs.TLD )
}

rmse <- function(sim, obs) 
{
  return( sqrt(mean((sim - obs)**2, na.rm = T)) )
}


#_____________________________________________________________________________________________________________________#
# import
tg <- 25
Sprime.min.match.rate <- 0.0
psmc..run <- 20
psmc..ylim <- c(0, 100e3)
# analyze
psmc..xlim <- c(1e4, 8e6)
obs.data <- abc.import("archives/obs/obs", 
                       read.params = FALSE, 
                       do.AFS.2D = TRUE,
                       generation.time = tg, 
                       Sprime.min.match.rate = Sprime.min.match.rate, 
                       psmc..mut.rate = 1.2e-8, 
                       psmc..run = psmc..run)
# plot
scg <- scale_color_viridis(discrete = FALSE, option = "C", begin = .1, end = .9)
sfg <- scale_fill_viridis(discrete = TRUE, option = "C", begin = .1, end = .9, direction = -1)
sfg <- scale_fill_viridis(discrete = TRUE, option = "A", direction = -1)
point.size = 2.0
line.size = .7
squish = FALSE
#_____________________________________________________________________________________________________________________#

spar <- "param_files/Avvaiyar.spar"
N <- 20
DirOut <- "Final.Blake"

#__________________ANALYSIS___________________# 

if (F) {
  DATA <- abc.import("Final.Blake/data/1M", 
                     do.AFS.2D = T, 
                     read.params = T, 
                     generation.time = tg, 
                     Sprime.min.match.rate = Sprime.min.match.rate)
  # 10361
  # Remove duplicated parameters
  cat("\nSubsetting to non-duplicated parameters:\n")
  P1 <- DATA$PAR
  nori <- nrow(P1)
  NonDups <- (1:nrow(P1))[!duplicated(P1)]
  DATA <- extract(DATA, NonDups)
  cat("  Removed:  ",nori-nrow(DATA$PAR),"\n")
  cat("  Retained: ",nrow(DATA$PAR)," | ~",round(nrow(DATA$PAR)/nori*100,1),"%\n\n")
  # Analyze
  R <- abc.analyze(obs.data, 
                   DATA, 
                   parfile = spar, 
                   N = N, 
                   dispersion.fun = "mad", 
                   psmc..xlim = psmc..xlim)
  #---
  R.accepted <- R$ACCEPTED
  fl <- sapply(R.accepted$name, function(x) rev(strsplit(x, "/")[[1]])[1])
  write.table(fl, paste0(DirOut,"/FINAL.accepted"), row.names=F, col.names=F, quote=F)
  #---
  abc.plot(obs.data, 
           DATA, 
           R.accepted, 
           output = paste0(DirOut,"/FINAL.tiledA.png"), 
           title = spar, 
           line.size = line.size, 
           point.size = point.size, 
           scg = scg,
           sfg = sfg,
           plot.each = F, 
           tile = T, 
           show.LD.obs = F, 
           psmc..ylim = psmc..ylim,
           xlim.D = c(-0.02, .2),
           squeeze = squish,
           layout = "A")
  #---
  abc.plot(obs.data, 
           DATA, 
           R.accepted, 
           output = paste0(DirOut,"/FINAL.tiledB.png"), 
           title = spar, 
           line.size = line.size, 
           point.size = point.size, 
           scg = scg,
           sfg = sfg,
           plot.each = F, 
           tile = T, 
           show.LD.obs = F, 
           psmc..ylim = psmc..ylim,
           xlim.D = c(-0.02, .2),
           squeeze = squish,
           layout = "B")
  #---
  abc.plot(obs.data, 
           DATA, 
           R.accepted, 
           output = paste0(DirOut,"/FINAL.A.png"), 
           title = spar, 
           line.size = line.size, 
           point.size = point.size, 
           scg = scg,
           sfg = sfg,
           plot.each = T, 
           tile = F, 
           show.LD.obs = F, 
           psmc..ylim = psmc..ylim,
           xlim.D = c(-0.02, .2),
           squeeze = squish,
           layout = "A")
  #---
  do.abc.par(DATA, R, 
             paste0(DirOut,"/FINAL"), 
             prior.file = "bo_1k5f.est", 
             dir.priors = "Final.Blake/data/param_prior_sampling")
}

#_________________Model Comparison ::: 20x30Mbp__________________#

if (F) {
  AccNames <- unname(unlist(read.table(paste0(DirOut,"/FINAL.accepted"), header = F)))
  fun.summary <- "min"
  PubModels <- c("durvasula20", "fu14", "gower21", "iasi21" , "jacobs19", "kamm19", "moorjani16" , "ragsdale19", "schaefer21", "skov20" , "yang12")
  DirRT <- "Final.Blake/data/20x30Mbp.50sim/1M.accepted"
  
  #================> var mutrate ::: ScaledL2_CrossCor_2PSMC.spar
  DirPublished <- "Final.Blake/data/20x30Mbp.50sim/"
  Published <- paste0(DirPublished,PubModels)
  names(Published) <- gsub(DirPublished, "", Published, fixed = TRUE)
  Res3 <- compare(DirRT, AccNames, Published, N, spar, fun.summary,
                  comparison.spar = "param_files/ScaledL2_CrossCor_2PSMC.spar", 
                  outprefix = paste0(DirOut,"/FINAL.20x30Mbp.varMut.2PSMC"),
                  tg, Sprime.min.match.rate, psmc..xlim, do.pdf = TRUE )
  
  #================> var mutrate ::: ScaledL2_CrossCor_2PSMC_CRF.spar
  DirPublished <- "Final.Blake/data/20x30Mbp.50sim/"
  Published <- paste0(DirPublished,PubModels)
  names(Published) <- gsub(DirPublished, "", Published, fixed = TRUE)
  Res4 <- compare(DirRT, AccNames, Published, N, spar, fun.summary,
                  comparison.spar = "param_files/ScaledL2_CrossCor_2PSMC_CRF.spar", 
                  outprefix = paste0(DirOut,"/FINAL.20x30Mbp.varMut.2PSMC.CRF"), 
                  tg, Sprime.min.match.rate, psmc..xlim )
  #
  Distances <- Res4$unscaled.distances
  RawStats <- Res4$all.raw.stats
  #
  AllStats <- Res4$split.by.statistics[[1]]
  AllStatsRank <- apply(AllStats, 2, rank)
  write.table(AllStats, paste0(DirOut,"/FINAL.20x30Mbp.varMut.2PSMC.CRF.ScaledDistPerStat.txt"))
  write.table(AllStatsRank, paste0(DirOut,"/FINAL.20x30Mbp.varMut.2PSMC.CRF.ScaledDistPerStat.txt"), append=T, col.names=F)
  #
  RPE <- Res4$relative.percentage.errors
  write.table(RPE, paste0(DirOut,"/FINAL.20x30Mbp.varMut.2PSMC.CRF.RPE.txt"))
  #
  CNAMES <- colnames(RPE)[-c(1,2)]
  CNAMES <- CNAMES[!grepl("FS.", CNAMES, fixed=T)]
  pdf(paste0(DirOut,"/FINAL.20x30Mbp.varMut.2PSMC.CRF.RPE.pdf"), height = 10, width = 3)
  for (cname in CNAMES) {
    df <- RPE[,c("model", "run", cname)]
    names(df) <- c("model", "run", "y")
    rg <- max(abs(df$y), na.rm=T)
    #
    su <- df %>%
      group_by(model) %>% 
      summarize(value = round(mean(y, na.rm = TRUE),1))
    su <- as.data.frame(su)
    su$col <- sapply(su$value, function(a) ifelse(a<0,"-","+"))
    su$value <- paste0(su$value,"%")
    #
    g=ggplot(df) +
      geom_histogram(aes(y, fill=model), bins = 30) +
      theme_bw() +
      facet_wrap(~model, ncol=1) +
      scale_x_continuous(limits = c(-1*rg, rg)) +
      geom_vline(xintercept = 0) +
      theme(legend.position = "none") +
      ggtitle(cname) +
      geom_text(data=su, aes(x=Inf, y=Inf, label=value, col=col), hjust=1, vjust=1, size=6) +
      xlab("") + ylab("")
    print(g)
  }
  dev.off()
  
  sink(paste0(DirOut,"/FINAL.20x30Mbp.COMPARISON.varMut.2PSMC.CRF.Statistics.txt"))
  # 
  cat("\nModels with RPE > 30% on either pi:CEU or pi:YRI:\n")
  A1 <- as.matrix(by(RPE$PI.CEU, RPE$model, mean))
  A2 <- as.matrix(by(RPE$PI.YRI, RPE$model, mean))
  A <- cbind(A1, A2[rownames(A1),1])
  colnames(A) <- c("PI.CEU", "PI.YRI")
  A <- as.data.frame(A)
  A$ok <- apply(A, 1, function(a) ifelse(any(abs(a)>30), "pb", "ok"))
  A$model <- rownames(A)
  A <- subset(A, model != "z. This study")
  A <- A[order(A$ok),]
  print(A)
  print(sum(A$ok=="pb"))
  #
  cat("\nModels with mean m_SPRIME < 60%\n")
  A1 <- as.data.frame(as.matrix(by(RawStats$SPRIME.match.mean, RawStats$model, mean, na.rm=T)))
  A1$model <- rownames(A1)
  A1 <- subset(A1, model != "z. This study")
  A1 <- A1[order(A1$V1),]
  print(A1[A1$V1<0.60,])
  print(nrow(A1[A1$V1<0.60,]))
  # 
  cat("\nOrdered mean values of alpha_CRF\n")
  A1 <- as.data.frame(as.matrix(by(RawStats$CRF.alpha.mean, RawStats$model, mean, na.rm=T)))
  A1$model <- rownames(A1)
  A1 <- subset(A1, model != "z. This study")
  A1 <- A1[order(A1$V1),]
  print(A1)
  #
  cat("\nModels with RPE > 50% on alpha_CRF: \n")
  A1 <- as.matrix(by(RPE$CRF.alpha.mean, RPE$model, mean))
  A1 <- as.data.frame(A1)
  A1$ok <- apply(A1, 1, function(a) ifelse(any(abs(a)>50), "pb", "ok"))
  A1$model <- rownames(A1)
  A1 <- subset(A1, model != "z. This study")
  A1 <- A1[order(A1$ok),]
  print(A1)
  print(sum(A1$ok=="pb"))
  #
  cat("\nMean RPE for Fst: \n")
  A <- as.data.frame(as.matrix(by(RPE$Fst.1, RPE$model, mean)))
  A$model <- rownames(A)
  A <- subset(A, model != "z. This study")
  A <- A[order(A[,1]),]
  print(A)
  #
  cat("\nModels with RPE > 30% on D:\n")
  A <- as.matrix(by(RPE$D.D, RPE$model, mean))
  colnames(A) <- "D.D"
  A <- as.data.frame(A)
  A$ok <- apply(A, 1, function(a) ifelse(any(abs(a)>30), "pb", "ok"))
  A$model <- rownames(A)
  A <- subset(A, model != "z. This study")
  A <- A[order(A$ok, A[,1]),]
  print(A)
  print(sum(A$ok=="pb"))
  #
  cat("\nMean T_LD (kya):\n")
  A <- as.data.frame((1e-3*as.matrix(by(RawStats$LD.mean, RawStats$model, mean))))
  A$model <- rownames(A)
  A <- subset(A, model != "z. This study")
  A <- A[order(A[,1]),]
  print(as.data.frame((A)))
  # 
  cat("\n Mean coef distance for PSMC CEU per model:\n")
  A <- as.data.frame(as.matrix(by(1-Distances$PSMC.CEU, Distances$model, mean)))
  A$model <- rownames(A)
  A <- A[order(A[,1]),]
  print(A)
  # 
  cat("\n Mean coef distance for PSMC YRI per model:\n")
  A <- as.data.frame(as.matrix(by(1-Distances$PSMC.YRI, Distances$model, mean)))
  A$model <- rownames(A)
  A <- A[order(A[,1]),]
  print(A)
  # 
  cat("\n Max coef distance for PSMC CEU per model:\n")
  A <- as.data.frame(as.matrix(by(1-Distances$PSMC.CEU, Distances$model, max)))
  A$model <- rownames(A)
  A <- A[order(A[,1]),]
  print(A)
  # 
  cat("\n Max coef distance for PSMC YRI per model:\n")
  A <- as.data.frame(as.matrix(by(1-Distances$PSMC.YRI, Distances$model, max)))
  A$model <- rownames(A)
  A <- A[order(A[,1]),]
  print(A)
  sink()
    
}

#_________________Model Comparison.HOTSPOTS ::: 20x30Mbp.HOTSPOTS__________________# 

if (F) {
  AccNames <- unname(unlist(read.table(paste0(DirOut,"/FINAL.accepted"), header = F)))
  fun.summary <- "min"
  PubModels <- c("durvasula20", "fu14", "gower21", "iasi21" , "jacobs19", "kamm19", "moorjani16" , "ragsdale19", "schaefer21", "skov20" , "yang12")
  DirRT <- "Final.Blake/data/20x30Mbp.50sim.HOTSPOTS/1M.accepted"
  
  #================> var mutrate ::: ScaledL2_CrossCor_2PSMC.spar
  DirPublished <- "Final.Blake/data/20x30Mbp.50sim.HOTSPOTS/"
  Published <- paste0(DirPublished,PubModels)
  names(Published) <- gsub(DirPublished, "", Published, fixed = TRUE)
  Res3 <- compare(DirRT, AccNames, Published, N, spar, fun.summary,
                  comparison.spar = "param_files/ScaledL2_CrossCor_2PSMC.spar", 
                  outprefix = paste0(DirOut,"/FINAL.20x30Mbp.varMut.2PSMC.HOTSPOTS"), 
                  tg, Sprime.min.match.rate, psmc..xlim )
}

#__________________20x30Mbp___________________#

if (F) plot_set("Final.Blake/data/20x30Mbp", paste0(DirOut,"/FINAL.20x30Mbp"), spar, N, do.pdf = TRUE)

#__________________20x30Mbp.HOTSPOTS___________________#

if (F) plot_set("Final.Blake/data/20x30Mbp.HOTSPOTS", paste0(DirOut,"/FINAL.20x30Mbp.HOTSPOTS"), spar, N)

#__________________Statistics ::: 20x30Mbp___________________# 

if (F) {
  AccNames <- unname(unlist(read.table(paste0(DirOut,"/FINAL.accepted"), header = F)))
  DATA <- abc.import("Final.Blake/data/20x30Mbp", 
                     files = AccNames,
                     do.AFS.2D = F, 
                     read.params = F, 
                     generation.time = tg, 
                     Sprime.min.match.rate = Sprime.min.match.rate)
  R <- abc.analyze(obs.data, 
                   DATA, 
                   parfile = "param_files/ScaledL2_CrossCor_2PSMC.spar", 
                   N = length(AccNames), 
                   dispersion.fun = "mad", 
                   psmc..xlim = psmc..xlim)
  #___
  sink(paste0(DirOut,"/FINAL.Statistics.txt"))
  sumup(label = "D", DATA$SS$D[,"D",drop=F], obs.data$SS$D[,"D",drop=F], precision = 1, as.percentage = T)
  sumup(label = "L_S'", DATA$SS$SPRIME[,"length.mean",drop=F], obs.data$SS$SPRIME[,"length.mean",drop=F], precision = 0, as.percentage = F)
  sumup(label = "m_S'", DATA$SS$SPRIME[,"match.mean",drop=F], obs.data$SS$SPRIME[,"match.mean",drop=F], precision = 1, as.percentage = T)
  sumup(label = "L_CRF", DATA$SS$CRF[,"length.mean",drop=F], obs.data$SS$CRF[,"length.mean",drop=F], precision = 0, as.percentage = F)
  sumup(label = "α_CRF", DATA$SS$CRF[,"alpha.mean",drop=F], obs.data$SS$CRF[,"alpha.mean",drop=F], precision = 1, as.percentage = T)
  sumup(label = "T_LD", DATA$SS$LD[,"mean",drop=F], obs.data$SS$LD[,"mean",drop=F], precision = 0, as.percentage = F)
  sumup(label = "π_CEU", DATA$SS$PI[,"CEU",drop=F], obs.data$SS$PI[,"CEU",drop=F], precision = 3, as.percentage = F)
  sumup(label = "π_YRI", DATA$SS$PI[,"YRI",drop=F], obs.data$SS$PI[,"YRI",drop=F], precision = 3, as.percentage = F)
  cat("\n\n PSMC CEU\nCrosscor=",round(mean(1-R$DISTANCES$raw[,"PSMC.CEU"]),2),"±",round(sd(1-R$DISTANCES$raw[,"PSMC.CEU"]),2),sep="")
  cat("\n\n PSMC YRI\nCrosscor=",round(mean(1-R$DISTANCES$raw[,"PSMC.YRI"]),2),"±",round(sd(1-R$DISTANCES$raw[,"PSMC.YRI"]),2),sep="")
  sink()
  
}
  
#__________________Statistics ::: 20x30Mbp.HOTSPOTS___________________# 

if (F) {
  AccNames <- unname(unlist(read.table(paste0(DirOut,"/FINAL.accepted"), header = F)))
  DATA <- abc.import("Final.Blake/data/20x30Mbp.HOTSPOTS", 
                     files = AccNames,
                     do.AFS.2D = F, 
                     read.params = F, 
                     generation.time = tg, 
                     Sprime.min.match.rate = Sprime.min.match.rate)
  R <- abc.analyze(obs.data, 
                   DATA, 
                   parfile = "param_files/ScaledL2_CrossCor_2PSMC.spar", 
                   N = length(AccNames), 
                   dispersion.fun = "mad", 
                   psmc..xlim = psmc..xlim)
  #___
  sink(paste0(DirOut,"/FINAL.Statistics.HOTSPOTS.txt"))
  sumup(label = "D", DATA$SS$D[,"D",drop=F], obs.data$SS$D[,"D",drop=F], precision = 1, as.percentage = T)
  sumup(label = "L_S'", DATA$SS$SPRIME[,"length.mean",drop=F], obs.data$SS$SPRIME[,"length.mean",drop=F], precision = 0, as.percentage = F)
  sumup(label = "m_S'", DATA$SS$SPRIME[,"match.mean",drop=F], obs.data$SS$SPRIME[,"match.mean",drop=F], precision = 1, as.percentage = T)
  sumup(label = "L_CRF", DATA$SS$CRF[,"length.mean",drop=F], obs.data$SS$CRF[,"length.mean",drop=F], precision = 0, as.percentage = F)
  sumup(label = "α_CRF", DATA$SS$CRF[,"alpha.mean",drop=F], obs.data$SS$CRF[,"alpha.mean",drop=F], precision = 1, as.percentage = T)
  sumup(label = "T_LD", DATA$SS$LD[,"mean",drop=F], obs.data$SS$LD[,"mean",drop=F], precision = 0, as.percentage = F)
  sumup(label = "π_CEU", DATA$SS$PI[,"CEU",drop=F], obs.data$SS$PI[,"CEU",drop=F], precision = 3, as.percentage = F)
  sumup(label = "π_YRI", DATA$SS$PI[,"YRI",drop=F], obs.data$SS$PI[,"YRI",drop=F], precision = 3, as.percentage = F)
  cat("\n\n PSMC CEU\nCrosscor=",round(mean(1-R$DISTANCES$raw[,"PSMC.CEU"]),2),"±",round(sd(1-R$DISTANCES$raw[,"PSMC.CEU"]),2),sep="")
  cat("\n\n PSMC YRI\nCrosscor=",round(mean(1-R$DISTANCES$raw[,"PSMC.YRI"]),2),"±",round(sd(1-R$DISTANCES$raw[,"PSMC.YRI"]),2),sep="")
  sink()
  
}

#__________________aDNA ::: 20x30Mbp ::: stats_aDNA ::: sampling in original pSampleCEU___________________# 

if (F)  {
  Dir <- "Final.Blake/further/20x30Mbp.U/stats_aDNA"; OutFile <- paste0(DirOut,"/FINAL.aDNA.png")
  #___
  DATA <- abc.import(Dir,
                     do.AFS.2D = F, 
                     read.params = F, 
                     generation.time = tg, 
                     Sprime.min.match.rate = Sprime.min.match.rate,
                     ld.exclude.negative.ic95 = FALSE)
  NAMES <- unname(sapply(DATA$NAME, function(x) rev(strsplit(x, "/", fixed=T)[[1]])[1]))
  RUN <- sapply(NAMES, function(x) strsplit(x, "__")[[1]][1])
  DEME <- sapply(NAMES, function(x) strsplit(strsplit(x, "__")[[1]][2], ".", fixed=T)[[1]][3])
  AGE <- sapply(NAMES, function(x) strsplit(strsplit(x, "__")[[1]][2], ".", fixed=T)[[1]][4])
  #___
  # D
  df <- cbind(run=RUN, deme=DEME, age=AGE, y=DATA$SS$D[,"D"], Z=DATA$SS$D[,"Z"])
  df <- na.omit(as.data.frame(apply(df, 2, as.numeric)))
  df.D <- df
  gD=ggplot(df.D) +
    stat_lineribbon(aes(x=age*1e-3, y=y), .width = c(0.5, 0.75, 0.95, 1.0), size = .7, alpha = 1.) +
    theme_ggdist() +
    scale_fill_paletteer_d("rcartocolor::BrwnYl", direction = 1) +
    facet_wrap(~deme, nrow = 1) +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    #scale_y_continuous(limits = c(0, .2), labels = scales::percent_format(accuracy = 1)) +
    xlab("Sampling age of Eurasian individual (kya)") + ylab("D") +
    ##
    geom_point(x=45, y=0.068, col="dodgerblue4") + 
    geom_errorbar(aes(x=45, ymin=0.068-1.96*0.008, ymax=0.068+1.96*0.008), col="dodgerblue4", width=0) + 
    geom_text(x=45, y=0.068, label="Ust\'-Ishim", hjust=0, vjust=0, angle=90, col="darkgray", alpha=.8) +
    ##
    theme(legend.position = "none")
  #___
  # LD  
  df <- cbind(run=RUN, deme=DEME, age=AGE, y=DATA$SS$LD[,"mean"], lb=DATA$SS$LD[,"mean"]-1.96*DATA$SS$LD[,"SE"])
  df <- na.omit(as.data.frame(apply(df, 2, as.numeric)))
  df$yBS <- df$y
  df$y <- df$y + df$age
  df.LD <- df
  gLD=ggplot(df.LD) +
    stat_lineribbon(aes(x=age*1e-3, y=y*1e-3), .width = c(0.5, 0.75, 0.95, 1.0), size = .7, alpha = 1.) +
    theme_ggdist() +
    scale_fill_paletteer_d("rcartocolor::PurpOr", direction = 1) +
    facet_wrap(~deme, nrow = 1) +
    geom_hline(yintercept = 50e3*1e-3, linetype = "dotted") +
    xlab("Sampling age of Eurasian individual (kya)") + ylab("Ancestry LD decay rate (kya)") +
    scale_y_log10(limits = c(1., 1e3)) +
    theme(legend.position = "none")
  #___
  # f4
  df <- cbind(run=RUN, deme=DEME, age=AGE, y=DATA$SS$f4[,"f4"])
  df <- na.omit(as.data.frame(apply(df, 2, as.numeric)))
  df.f4 <- df
  gf4=ggplot(df.f4) +
    stat_lineribbon(aes(x=age*1e-3, y=y), .width = c(0.5, 0.75, 0.95, 1.0), size = .7, alpha = 1.) +
    theme_ggdist() +
    scale_fill_paletteer_d("rcartocolor::BrwnYl", direction = 1) +
    facet_wrap(~deme, nrow = 1) +
    scale_y_continuous(limits = c(0, .1), labels = scales::percent_format(accuracy = 1)) +
    xlab("Sampling age of Eurasian individual (kya)") + ylab("f4") + 
    ##
    geom_point(x=39.405, y=0.064, col="dodgerblue4") + 
    geom_errorbar(aes(x=39.405, ymin=0.057, ymax=0.071), col="dodgerblue4", width=0) + 
    geom_errorbarh(aes(y=0.064, xmin=36.950, xmax=41.860), col="dodgerblue4", height=0) + 
    geom_text(x=39.405, y=0.064, label="Oase1", hjust=0, vjust=0, angle=90, col="darkgray", alpha=.8) +
    ##
    geom_point(x=45, y=0.026, col="dodgerblue4") + 
    geom_errorbar(aes(x=45, ymin=0.026-1.96*0.0033, ymax=0.026+1.96*0.0033), col="dodgerblue4", width=0) + 
    geom_text(x=45, y=0.026, label="Ust\'-Ishim", hjust=0, vjust=0, angle=90, col="darkgray", alpha=.8) +
    ##
    theme(legend.position = "none")
  #___
  G <- cowplot::plot_grid(gD, gLD, gf4, nrow = 3)
  ggsave(G, filename = OutFile, width = 13, height = 10)
  #G <- cowplot::plot_grid(gD, gLD, nrow = 2)
  #ggsave(G, filename = paste0("../plots_080622/",dir,".aDNA.png"), width = 13, height = 8)
  #___
  # Statistics
  sink(paste0(DirOut,"/FINAL.Statistics.aDNA.AllSNPs.txt"))
  cat("\n\n________________________ D ________________________\n\n")
  print(summary(lm(y ~ deme + age, data = df.D)))
  cat("\n\n")
  sumup(label = "D", as.matrix(df.D[,"y",drop=F]), obs.data$SS$D[,"D",drop=F], precision = 1, as.percentage = T)
  sumup(label = "Z", as.matrix(df.D[,"Z",drop=F]), obs.data$SS$D[,"Z",drop=F], precision = 1, as.percentage = 0)
  #
  cat("\n\n________________________ LD yBS ________________________\n\n")
  print(summary(lm(yBS ~ deme + age, data = df.LD)))
  cat("\n\n________________________ LD yBP ________________________\n\n")
  print(summary(lm(y ~ deme + age, data = df.LD)))
  cat("\n\n")
  sumup(label = "LD yBP", as.matrix(df.LD[,"y",drop=F]), obs.data$SS$LD[,"mean",drop=F], precision = 0, as.percentage = F)
  ###
  cat("\n\n\n\n=================================\n")
  cat("__________________ONLY DEME 3______________\n\n\n")
  df.LD <- subset(df.LD, deme == 3)
  cat("\n\n________________________ LD yBS ________________________\n\n")
  print(summary(lm(yBS ~ age, data = df.LD)))
  cat("\n\n________________________ LD yBP ________________________\n\n")
  print(summary(lm(y ~ age, data = df.LD)))
  cat("\n\n")
  sumup(label = "LD yBP", as.matrix(df.LD[,"y",drop=F]), obs.data$SS$LD[,"mean",drop=F], precision = 0, as.percentage = F)
  sink()
  #___
  # LD ::: Deme 3
  df.LD <- subset(df.LD, deme == 3)
  Runs <- unique(df.LD$run)
  PAR <- sapply(Runs, function(x) {
    xx <- read.table(paste0("Final.Blake/data/20x30Mbp/",x,".par"), header=F)
    rownames(xx) <- xx[,1]
    return(xx["T.nea",2])
    })
  PAR <- data.frame(run=Runs, par=PAR)
  df <- merge(df.LD, PAR, by = "run", all.x = T, sort = F)
  g=ggplot(df, aes(y=factor(run))) +
    geom_violin(aes(x=y*1e-3), alpha = .5, fill = "gray", col = "gray70") +
    geom_point(aes(x=y*1e-3, col=age*1e-3), size = 2.5, alpha = .8) +
    geom_point(aes(x=par*1e-3), shape = 4, fill = "black", alpha = .8, size = 3.5) +
    theme_classic() +
    scale_x_continuous(limits = c(0, 1.1e6*1e-3), expand = c(0, 0)) +
    theme(panel.grid.major.y = element_line(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 13)) +
    labs(x = "Time (kya)",
         y = "",
         title = expression(paste("Mean decay constant of the single-sample statistics: ",T[LD][","][1])),
         caption = expression(Crosses~represent~the~time~of~split~between~M[B]~and~M[N]),
         color = "Sample\nAge (kya)") +
    scale_color_viridis(option = "B"); g
  ggsave(g, filename = paste0(DirOut,"/FINAL.Statistics.aDNA.Deme3.LD.png"), width = 10, height = 8)
}

#__________________aDNA ::: 20x30Mbp ::: with observed aDNA___________________# 

if (F) {
  
  DemeIdx <- 3
  Dir <- paste0("Final.Blake/further/20x30Mbp.U.CEU",DemeIdx)
  SubDirs <- c("stats_aDNA.di.All", "stats_aDNA.pseudodi.Archaic", "stats_aDNA.di.1M", "stats_aDNA.pseudodi.1M")
  DirOut <- "Final.Blake"
  #
  obs_ind <- "archives/obs/aDNA/v54.1_1240K_public.ind"
  obs_anno <- "archives/obs/aDNA/v54.1_1240K_public.anno"
  obs_stats <- archives/obs/aDNA/stats_aDNA"
  #
  schaefer21 <- "Final.Blake/further/20x30Mbp.U/stats_aDNA_schaefer21"
  
  
  #______Observed aDNA
  IND <- read.table(obs_ind, sep = "", header=F, stringsAsFactors = F)
  ANNO <- read.delim2(obs_anno, sep = "\t", header = T, na.strings="()()()", comment.char = "", stringsAsFactors = F, quote = "")
  ANNO <- data.frame(GroupID=ANNO$Group.ID, C14.mean=ANNO[,8], C14.sd=ANNO[,9], YearPubli=ANNO[,5], 
                     Loc=ANNO[,13], Country=ANNO[,14], Source=ANNO[,18], Master.ID=ANNO[,2], Genetic.ID=ANNO[,1])
  ANNO$GroupID <- IND$V3
  #
  lf <- list.files(obs_stats, pattern = ".stats")
  lf <- lf[!grepl(".cov.stats", lf)]
  obs <- lapply(lf, function(x) {
    x <- scan(paste0(obs_stats,"/",x), what = "character", sep = "\n", quiet = T)[5]
    x <- strsplit(x, " ")[[1]]
    if (x[1] != "D_diploid") stop("Error")
    return(as.numeric(x[c(3,4)]))
  })
  obs <- do.call(rbind, obs)
  obs <- data.frame(SPECIMEN=gsub(".stats", "", lf), D=obs[,1], SE=obs[,2])
  # merge with ANNO
  obs <- merge(obs, ANNO, by.x="SPECIMEN", by.y="GroupID", all.x=T, all.y=F, sort=F)
  obs$SPECIMEN <- as.character(obs$SPECIMEN)
  y <- rev(sort(table(obs$SPECIMEN)))
  if (any(y[!names(y)%in%"French.DG"]!=1)) stop("PB")
  obs <- obs[!duplicated(obs$SPECIMEN),]
  obs$SPECIMEN <- factor(obs$SPECIMEN, levels = obs$SPECIMEN[order(obs$C14.mean)])
  # remove Romania_Oase_UP
  obs <- subset(obs, SPECIMEN != "Romania_Oase_UP")
  # scale so that French's D = 0.0457
  obs$D.scaled <- obs$D / obs[obs$SPECIMEN=="French.DG","D"] * obs.data$SS$D[,"D"]
  # IC95
  obs$C14.low <- obs$C14.mean - 1.96*obs$C14.sd
  obs$C14.up <- obs$C14.mean + 1.96*obs$C14.sd
  obs$D.low <- obs$D.scaled - 1.96*obs$SE
  obs$D.up <- obs$D.scaled + 1.96*obs$SE
  write.table(obs, paste0(DirOut,"/aDNA.obs.txt"), quote=F, row.names=F, sep="\t")
  
  
  #______Simulated aDNA
  aDNA <- c()
  for (dir in SubDirs) {
    DATA <- abc.import(paste0(Dir,"/",dir),
                       do.AFS.2D = F, 
                       read.params = F, 
                       generation.time = tg, 
                       Sprime.min.match.rate = Sprime.min.match.rate,
                       ld.exclude.negative.ic95 = FALSE)
    NAMES <- unname(sapply(DATA$NAME, function(x) rev(strsplit(x, "/", fixed=T)[[1]])[1]))
    RUN <- sapply(NAMES, function(x) strsplit(x, "__")[[1]][1])
    DEME <- sapply(NAMES, function(x) strsplit(strsplit(x, "__")[[1]][2], ".", fixed=T)[[1]][3])
    AGE <- as.numeric(sapply(NAMES, function(x) strsplit(strsplit(x, "__")[[1]][2], ".", fixed=T)[[1]][4]))
    aDNA <- rbind(aDNA, data.frame(dir=dir, run=RUN, deme=DEME, age=AGE, D=DATA$SS$D[,"D"], Z=DATA$SS$D[,"Z"]))
  }
  
  
  #______Schaefer21 aDNA
    lf <- list.files(schaefer21, pattern = ".stats")
  lf <- lf[!grepl(".cov.stats", lf)]
  S21 <- lapply(lf, function(x) {
    x <- scan(paste0(schaefer21,"/",x), what = "character", sep = "\n", quiet = T)[5]
    x <- strsplit(x, " ")[[1]]
    if (x[1] != "D_diploid") stop("Error")
    return(as.numeric(x[c(3,4)]))
  })
  S21 <- do.call(rbind, S21)
  S21 <- data.frame(SPECIMEN=gsub(".stats", "", lf), D=S21[,1], SE=S21[,2], stringsAsFactors = F)
  S21$age <- as.numeric(sapply(S21$SPECIMEN, function(x) strsplit(x, ".", fixed=T)[[1]][4]))
  
  
  #______Plot!
  rg.ribbon <- 1.0
  rg.ribbon2 <- 0.95
  #
  aDNAs <- subset(aDNA, deme == DemeIdx)
  #
  X=aDNAs %>%
    group_by(dir, age) %>% 
    summarize(ymin = quantile(D, 0+(1-rg.ribbon)/2),
              mean = mean(D),
              ymax = quantile(D, 1-(1-rg.ribbon)/2),
              y2.5  = quantile(D, 0+(1-rg.ribbon2)/2),
              y97.5 = quantile(D, 1-(1-rg.ribbon2)/2) )
  X <- as.data.frame(X)
  
  # A => but small bug with the stat_lineribbon
  g=ggplot(aDNAs) +
    #theme_classic() +
    scale_y_continuous(limits = c(0, .23), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 1, decimal.mark = ','),
                       limits = c(-1000, 47e3)) +
    
    stat_lineribbon(aes(x=age, y=D, group=dir, fill=dir), .width = c(rg.ribbon), size = .7, alpha = .3) +
    stat_lineribbon(data = subset(aDNAs, dir == "stats_aDNA.pseudodi.Archaic" & age>=35e3 & age<50e3), 
                    aes(x=age, y=D, group=dir, fill=dir), .width = c(rg.ribbon), size = .7, alpha = .5) +
    
    #scale_fill_viridis(discrete = T, direction = -1) +
    scale_fill_manual(values = c("#4d749b", "#E8891D", "#05A3A4", "#006373")) +
    
    #scale_color_uchicago() +
    
    labs(x = "Sample 14C age (years BP)", 
         y = "D", 
         caption = paste0("Sampled CEU deme: ",DemeIdx+1)) +
    theme_light(base_size = 15, base_family = "Poppins") +
    theme(#legend.position = "none",
      axis.title = element_text(size = 12),
      axis.text.x = element_text(family = "Roboto Mono", size = 10),
      plot.caption = element_text(size = 9, color = "gray50"),
      panel.grid = element_blank()) +
    
    geom_line(data = S21, aes(x = age, y = D), linetype = "dashed", alpha=.7) +
    #geom_ribbon(data = S21, aes(x = age, ymin = D-1.96*SE, ymax = D+1.96*SE), fill = "lightgray") +
    #geom_line(data = S21, aes(x = age, y = D-1.96*SE), linetype = "dotted", alpha=.7) +
    #geom_line(data = S21, aes(x = age, y = D+1.96*SE), linetype = "dotted", alpha=.7) +
    #geom_errorbar(data = S21, aes(x = age, ymin = D-1.96*SE, ymax = D+1.96*SE), width = 20, alpha=.5) +
    
    geom_errorbarh(data = obs, aes(xmin=C14.low, xmax=C14.up, y=D.scaled), height=0, alpha=.7, size=.9) +
    geom_errorbar(data = obs, aes(x=C14.mean, ymin=D.low, ymax=D.up), width=0, alpha=.7, size=.9) +
    geom_point(data = obs, aes(x=C14.mean, y=D.scaled), size = 4, col = "white") +
    geom_point(data = obs, aes(x=C14.mean, y=D.scaled, col=SPECIMEN), size = 3) +
    geom_point(data = obs, aes(x=C14.mean, y=D), size = 3, alpha = .7) +
    geom_point(data = obs, aes(x=C14.mean, y=D, col=SPECIMEN), size = 1.5, alpha = 1.0)
  
  # B
  g=ggplot(X) +
    geom_ribbon(aes(x = age, ymin = ymin, ymax = ymax, fill = dir), alpha = .05) +
    geom_ribbon(data = subset(X, dir == "stats_aDNA.pseudodi.Archaic" & age>=35e3 & age<50e3),
                aes(x = age, ymin = ymin, ymax = ymax, fill = dir), alpha = .2 ) +
    geom_line(aes(x = age, y = ymin, col = dir), linetype = "solid") +
    geom_line(aes(x = age, y = ymax, col = dir), linetype = "solid") +
    geom_line(aes(x = age, y = mean, col = dir), size = 1.2) +
    #coord_cartesian(ylim = c(0, .23)) +#, labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 1, decimal.mark = ','),
                       limits = c(-1000, 47e3)) +
    scale_color_manual(values = c("#4d749b", "#E8891D", "#05A3A4", "#006373")) +
    scale_fill_manual(values = c("#4d749b", "#E8891D", "#05A3A4", "#006373")) +
    labs(x = "Sample 14C age (years BP)", 
         y = "D", 
         caption = "Sampled CEU deme: 4th",
         color = "Data",
         fill = "Data") +
    theme_light(base_size = 15, base_family = "Poppins") +
    theme(#legend.position = "none",
      axis.title = element_text(size = 12),
      axis.text.x = element_text(family = "Roboto Mono", size = 10),
      plot.caption = element_text(size = 9, color = "gray50"),
      panel.grid = element_blank()) +
    
    geom_line(data = S21, aes(x = age, y = D), linetype = "dashed", alpha=.7) +
    
    new_scale_colour() +
    #scale_color_uchicago() + 
    
    geom_errorbarh(data = obs, aes(xmin=C14.low, xmax=C14.up, y=D.scaled), height=0, alpha=.7, size=.9) +
    geom_errorbar(data = obs, aes(x=C14.mean, ymin=D.low, ymax=D.up), width=0, alpha=.7, size=.9) +
    geom_segment(data = obs, aes(x=C14.mean, xend=C14.mean, y=D, yend=D.scaled), linetype = "dotted") +
    geom_point(data = obs, aes(x=C14.mean, y=D.scaled), size = 4, col = "black") +
    geom_point(data = obs, aes(x=C14.mean, y=D.scaled, col = SPECIMEN), size = 2.2) +
    #
    geom_point(data = obs, aes(x=C14.mean, y=D), size = 2, alpha = .7, col = "black")
  #geom_point(data = obs, aes(x=C14.mean, y=D, shape = SPECIMEN), size = 1.5, alpha = 1.0, col = "white") +
  #scale_shape_manual(values = 0:length(unique(obs$SPECIMEN)))
  
  # D
  XX <- X
  XX$facet <- ""
  XX[XX$dir=="stats_aDNA.pseudodi.Archaic","facet"] <- "Ascertained"
  obs$facet <- ""
  obs[obs$SPECIMEN=="Romania_Oase_UP_enhanced","facet"] <- "Ascertained"
  S21$facet <- ""
  g=ggplot(XX) +
    facet_wrap(~facet) +
    geom_ribbon(aes(x = age, ymin = ymin, ymax = ymax, fill = dir), alpha = .05) +
    geom_ribbon(data = subset(XX, dir == "stats_aDNA.pseudodi.Archaic" & age>=35e3 & age<50e3),
                aes(x = age, ymin = ymin, ymax = ymax, fill = dir), alpha = .2 ) +
    geom_line(aes(x = age, y = ymin, col = dir), linetype = "solid") +
    geom_line(aes(x = age, y = ymax, col = dir), linetype = "solid") +
    geom_line(aes(x = age, y = y2.5, col = dir), linetype = "dotdash") +
    geom_line(aes(x = age, y = y97.5, col = dir), linetype = "dotdash") +
    geom_line(aes(x = age, y = mean, col = dir), size = 1.2) +
    #coord_cartesian(ylim = c(0, .23)) +#, labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 1, decimal.mark = ','),
                       limits = c(-1000, 47e3)) +
    scale_color_manual(values = c("#4d749b", "#E8891D", "#05A3A4", "#006373")) +
    scale_fill_manual(values = c("#4d749b", "#E8891D", "#05A3A4", "#006373")) +
    labs(x = "Sample 14C age (years BP)", 
         y = "D", 
         caption = paste0("Sampled CEU deme: ",DemeIdx+1),
         color = "Data",
         fill = "Data") +
    theme_light(base_size = 15, base_family = "Poppins") +
    theme(#legend.position = "none",
      axis.title = element_text(size = 12),
      axis.text.x = element_text(family = "Roboto Mono", size = 10),
      plot.caption = element_text(size = 9, color = "gray50"),
      panel.grid = element_blank()) +
    
    geom_line(data = S21, aes(x = age, y = D), linetype = "dashed", alpha=.7) +
    
    new_scale_colour() +
    #scale_color_uchicago() + 
    
    geom_errorbarh(data = obs, aes(xmin=C14.low, xmax=C14.up, y=D.scaled), height=0, alpha=.7, size=.9) +
    geom_errorbar(data = obs, aes(x=C14.mean, ymin=D.low, ymax=D.up), width=0, alpha=.7, size=.9) +
    geom_segment(data = obs, aes(x=C14.mean, xend=C14.mean, y=D, yend=D.scaled), linetype = "dotted") +
    geom_point(data = obs, aes(x=C14.mean, y=D.scaled), size = 4, col = "black") +
    geom_point(data = obs, aes(x=C14.mean, y=D.scaled, col = SPECIMEN), size = 2.2) +
    #
    geom_point(data = obs, aes(x=C14.mean, y=D), size = 2, alpha = .7, col = "black")
  
  # Add visual representation of which deme is sampled
  Demes <- data.frame(x=seq(0, 9, 1), col = 0)
  Demes[4,"col"] <- 1
  Demes$col <- as.character(Demes$col)
  map_demes=ggplot(Demes) +
    theme_void() +
    theme(legend.position = "none") +
    geom_segment(y = 0, yend = 0, x = 0, xend = 9, col = "#909090") +
    geom_point(aes(x = x, col = col), y = 0, size = 3) +
    scale_color_manual(values = c("#c3c1c1", "#505050"))
  g2=g +
    annotation_custom(ggplotGrob(map_demes), xmin = 2e3 + 1e3, xmax = 15e3, ymin = 0.20, ymax = 0.25) +
    annotate("text", x = 2e3, y = 0.20, label = "CEU", hjust = 1, vjust = 0)
  
  ggsave(g, filename = paste0(DirOut,"/FINAL.aDNA.Obs_aDNA2.png"), width = 16, height = 8, scale = .9)
  gPDF <- g +
    theme_light(base_size = 15, base_family = NA) +
    theme(axis.title = element_text(size = 12),
          axis.text.x = element_text(family = NA, size = 10),
          plot.caption = element_text(size = 10, color = "gray50"),
          panel.grid = element_blank() )
  ggsave(gPDF, filename = paste0(DirOut,"/FINAL.aDNA.Obs_aDNA2.pdf"), width = 16, height = 8, scale = .9)
  
  # C
  XX <- X
  XX$facet <- ""
  XX[XX$dir=="stats_aDNA.pseudodi.Archaic","facet"] <- "Ascertained"
  obs$facet <- ""
  obs[obs$SPECIMEN=="Romania_Oase_UP_enhanced","facet"] <- "Ascertained"
  S21$facet <- ""
  g=ggplot(XX) +
    facet_wrap(~facet) +
    geom_ribbon(aes(x = age, ymin = ymin, ymax = ymax, fill = dir), alpha = .05) +
    geom_ribbon(data = subset(XX, dir == "stats_aDNA.pseudodi.Archaic" & age>=35e3 & age<50e3),
                aes(x = age, ymin = ymin, ymax = ymax, fill = dir), alpha = .2 ) +
    geom_line(aes(x = age, y = ymin, col = dir), linetype = "solid") +
    geom_line(aes(x = age, y = ymax, col = dir), linetype = "solid") +
    geom_line(aes(x = age, y = mean, col = dir), size = 1.2) +
    #coord_cartesian(ylim = c(0, .23)) +#, labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 1, decimal.mark = ','),
                       limits = c(-1000, 47e3)) +
    scale_color_manual(values = c("#4d749b", "#E8891D", "#05A3A4", "#006373")) +
    scale_fill_manual(values = c("#4d749b", "#E8891D", "#05A3A4", "#006373")) +
    labs(x = "Sample 14C age (years BP)", 
         y = "D", 
         caption = paste0("Sampled CEU deme: ",DemeIdx+1),
         color = "Data",
         fill = "Data") +
    theme_light(base_size = 15, base_family = "Poppins") +
    theme(#legend.position = "none",
      axis.title = element_text(size = 12),
      axis.text.x = element_text(family = "Roboto Mono", size = 10),
      plot.caption = element_text(size = 9, color = "gray50"),
      panel.grid = element_blank()) +
    
    geom_line(data = S21, aes(x = age, y = D), linetype = "dashed", alpha=.7) +
    
    new_scale_colour() +
    #scale_color_uchicago() + 
    
    geom_errorbarh(data = obs, aes(xmin=C14.low, xmax=C14.up, y=D.scaled), height=0, alpha=.7, size=.9) +
    geom_errorbar(data = obs, aes(x=C14.mean, ymin=D.low, ymax=D.up), width=0, alpha=.7, size=.9) +
    geom_segment(data = obs, aes(x=C14.mean, xend=C14.mean, y=D, yend=D.scaled), linetype = "dotted") +
    geom_point(data = obs, aes(x=C14.mean, y=D.scaled), size = 4, col = "black") +
    geom_point(data = obs, aes(x=C14.mean, y=D.scaled, col = SPECIMEN), size = 2.2) +
    #
    geom_point(data = obs, aes(x=C14.mean, y=D), size = 2, alpha = .7, col = "black")
  
  # Add visual representation of which deme is sampled
  Demes <- data.frame(x=seq(0, 9, 1), col = 0)
  Demes[4,"col"] <- 1
  Demes$col <- as.character(Demes$col)
  map_demes=ggplot(Demes) +
    theme_void() +
    theme(legend.position = "none") +
    geom_segment(y = 0, yend = 0, x = 0, xend = 9, col = "#909090") +
    geom_point(aes(x = x, col = col), y = 0, size = 3) +
    scale_color_manual(values = c("#c3c1c1", "#505050"))
  g2=g +
    annotation_custom(ggplotGrob(map_demes), xmin = 2e3 + 1e3, xmax = 15e3, ymin = 0.20, ymax = 0.25) +
    annotate("text", x = 2e3, y = 0.20, label = "CEU", hjust = 1, vjust = 0)
  
  ggsave(g, filename = paste0(DirOut,"/FINAL.aDNA.Obs_aDNA.png"), width = 16, height = 8, scale = .9)
  gPDF <- g +
    theme_light(base_size = 15, base_family = NA) +
    theme(axis.title = element_text(size = 12),
          axis.text.x = element_text(family = NA, size = 10),
          plot.caption = element_text(size = 10, color = "gray50"),
          panel.grid = element_blank() )
  ggsave(gPDF, filename = paste0(DirOut,"/FINAL.aDNA.Obs_aDNA.pdf"), width = 16, height = 8, scale = .9)
  #___
  sink(paste0(DirOut,"/FINAL.Statistics.aDNA.Obs_aDNA.txt"))
  cat("\n\n________________________ D ________________________\n\n")
  All <- subset(aDNA, dir=="stats_aDNA.di.All")
  print(summary(lm(D ~ age, data = All)))
  cat("\n\n")
  sumup(label = "D", as.matrix(All[,"D",drop=F]), obs.data$SS$D[,"D",drop=F], precision = 1, as.percentage = T)
  sumup(label = "Z", as.matrix(All[,"Z",drop=F]), obs.data$SS$D[,"Z",drop=F], precision = 1, as.percentage = 0)
  #
  cat("\n\nD ~ age, for Schaefer 21:\n\n")
  print(summary(lm(D ~ age, data = S21)))
  #
  cat("\n\nANOVA between stats_aDNA.di.All + stats_aDNA.di.1M + stats_aDNA.pseudodi.1M:\n\n")
  print(summary(aov(D ~ dir, data = subset(aDNAs, dir!="stats_aDNA.pseudodi.Archaic"))))
  #
  cat("\n\nRatio stats_aDNA.pseudodi.Archaic over stats_aDNA.di.1M:\n\n")
  M1 <- subset(aDNAs, dir == "stats_aDNA.pseudodi.Archaic")
  rownames(M1) <- apply(M1[,c("run","deme","age")], 1, paste0, collapse=".")
  M2 <- subset(aDNAs, dir == "stats_aDNA.di.1M")
  rownames(M2) <- apply(M2[,c("run","deme","age")], 1, paste0, collapse=".")
  M2 <- M2[rownames(M1),]
  MM <- as.matrix(M1$D / M2$D)
  sumup(label = "D", MM, as.matrix(NA), precision = 1, as.percentage = F)
  sink()
  
}

#__________________Sensitivity Analyses___________________# 

if (F) {
  #_______________
  #  MAC vs no-MAC
  multidata_plots(data_paths = c("Final.Blake/further/sensitivity/rt.1014230.MAC", 
                                 "Final.Blake/further/sensitivity/rt.1014230.noMAC"), 
                labels = c("Filtered", "Non-filtered"),                    
                output = "Final.Blake/further/MAC_noMAC.png",
                title = "",
                do_tile = FALSE,
                plot.legend = TRUE,
                color.lims = c(0.1, 0.7), 
                log10.crf = c(T, T),
                xylims.crf = c("xmin"=5e3, "xmax"=1e6, "ymin"=0.1e-2, "ymax"=10e-2),
                log10.sprime = c(T, F),
                xylims.sprime = c("xmin"=1e4, "xmax"=3e6, "ymin"=0., "ymax"=1.),
                ylim.afs.ceu = c(0, .35),
                ylim.afs.yri = c(0, .45),
                ylim.dcfs = c(0, .45),
                xlim.D = c(-0.02, .20),
                squeeze = squish,
                show.RMSE = FALSE )
  #_______________
  #  Neanderthal sampling ages
  multidata_plots(data_paths = c("Final.Blake/further/sensitivity/rt.1014230.MAC", 
                                 "Final.Blake/further/sensitivity/rt.1014230.Altai"), 
                  labels = c("50 kya", "130 kya"),                    
                  output = "Final.Blake/further/Nea_Ages.png",
                  title = "",
                  do_tile = FALSE,
                  plot.legend = TRUE,
                  color.lims = c(0.1, 0.7), 
                  log10.crf = c(T, T),
                  xylims.crf = c("xmin"=5e3, "xmax"=1e6, "ymin"=0.1e-2, "ymax"=10e-2),
                  log10.sprime = c(T, F),
                  xylims.sprime = c("xmin"=1e4, "xmax"=3e6, "ymin"=0., "ymax"=1.),
                  ylim.afs.ceu = c(0, .35),
                  ylim.afs.yri = c(0, .45),
                  ylim.dcfs = c(0, .45),
                  xlim.D = c(-0.02, .20),
                  squeeze = squish,
                  show.RMSE = FALSE )
  #_______________
  #  Ancestral allele misspecification errors
  multidata_plots(data_paths = c("Final.Blake/further/sensitivity/rt.1014230.MAC", 
                                 "Final.Blake/further/sensitivity/rt.1014230.AncError3p"), 
                  labels = c("0%", "3%"),                    
                  output = "Final.Blake/further/Ancestral_Error.png",
                  title = "",
                  do_tile = FALSE,
                  plot.legend = TRUE,
                  color.lims = c(0.1, 0.7), 
                  log10.crf = c(T, T),
                  xylims.crf = c("xmin"=5e3, "xmax"=1e6, "ymin"=0.1e-2, "ymax"=10e-2),
                  log10.sprime = c(T, F),
                  xylims.sprime = c("xmin"=1e4, "xmax"=3e6, "ymin"=0., "ymax"=1.),
                  ylim.afs.ceu = c(0, .35),
                  ylim.afs.yri = c(0, .45),
                  ylim.dcfs = c(0, .45),
                  xlim.D = c(-0.02, .20),
                  squeeze = squish,
                  show.RMSE = FALSE )
  #_______________
  #  Polyploidy
  multidata_plots(data_paths = c("Final.Blake/further/sensitivity/rt.1014230.noMAC", 
                                 "Final.Blake/further/sensitivity/rt.1014230.noMAC.pseudodiplo"), 
                  labels = c("Diploid", "Pseudodiploid"),                    
                  output = "Final.Blake/further/Pseudodiploidy.png",
                  title = "(no MAC)",
                  do_tile = FALSE,
                  plot.legend = TRUE,
                  color.lims = c(0.1, 0.7), 
                  log10.crf = c(T, T),
                  xylims.crf = c("xmin"=5e3, "xmax"=1e6, "ymin"=0.1e-2, "ymax"=10e-2),
                  log10.sprime = c(T, F),
                  xylims.sprime = c("xmin"=1e4, "xmax"=3e6, "ymin"=0., "ymax"=1.),
                  ylim.afs.ceu = c(0, .35),
                  ylim.afs.yri = c(0, .45),
                  ylim.dcfs = c(0, .45),
                  xlim.D = c(-0.02, .20),
                  squeeze = squish,
                  show.RMSE = FALSE )
  #_______________
  # Empirical recombination maps
  multidata_plots(data_paths = c("Final.Blake/further/sensitivity/rt.1014230.MAC",
                                 "Final.Blake/data/10x7Mbp.50sim.HOTSPOTS/1M.accepted/1014230",
                                 "Final.Blake/data/empirical_gmaps/HapMap", 
                                 "Final.Blake/data/empirical_gmaps/Spence19"), 
                  labels = c("Uniform", "Hotspots", "HapMap", "Spence19"),                    
                  output = "Final.Blake/further/EmpiricalMaps.png",
                  title = "[MAC]",
                  do_tile = FALSE,
                  plot.legend = TRUE,
                  color.lims = c(0.1, 0.7), 
                  log10.crf = c(T, T),
                  xylims.crf = c("xmin"=5e3, "xmax"=1e6, "ymin"=0.1e-2, "ymax"=10e-2),
                  log10.sprime = c(T, F),
                  xylims.sprime = c("xmin"=1e4, "xmax"=3e6, "ymin"=0., "ymax"=1.),
                  ylim.afs.ceu = c(0, .35),
                  ylim.afs.yri = c(0, .45),
                  ylim.dcfs = c(0, .45),
                  xlim.D = c(-0.02, .20),
                  squeeze = squish,
                  show.RMSE = FALSE )
  }

#__________________Collinearity___________________# 
# Collinearity between variables

if (F) {
  require(corrplot)
  DATA <- abc.import("Final.Blake/data/20x30Mbp", 
                     do.AFS.2D = T, 
                     read.params = T, 
                     generation.time = tg, 
                     Sprime.min.match.rate = Sprime.min.match.rate)
  COL <- abc.analyze(obs.data, 
                     DATA, 
                     parfile = "param_files/All.Collinearity.spar", 
                     N = NA, 
                     dispersion.fun = "mad", 
                     psmc..xlim = psmc..xlim)
  MR <- cor(COL$DISTANCES$raw, use = "pairwise.complete.obs")
  MP <- cor.mtest(COL$DISTANCES$raw, conf.level = 0.95)
  png(paste0("Final.Blake/Collinearity_5p.png"), width=1200, height=1200, pointsize = 28)
  corrplot(MR, order = 'hclust', diag = F, cl.pos = "n", number.cex = 0.8, p.mat = MP$p, 
           sig.level = 0.05, addrect = 2, 
           insig = "blank", addCoef.col = "black", tl.col = "black")
  dev.off()
}

#__________________LD Algorithms___________________# 

if (F) {
  
  require(minpack.lm)

  dir.data <- "Final.Blake/further/LD_implementation"
  nsim <- 30
  
  V=c()
  for (i in 1:nsim) {
    par <- read.table(paste0(dir.data,"/",i,".par"), header=F)
    sanka <- read.table(paste0(dir.data,"/",i,".AsSanka12.0.cov.txt"), header=T)
    rt <- read.table(paste0(dir.data,"/",i,".RT.0.cov.txt"), header=T)
    d <- read.table(paste0(dir.data,"/",i,".computed.txt"), header=F)
    moorjani <- read.table(paste0(dir.data,"/",i,".ld1.cov.txt.gz"), header=T)
    #
    M <- merge(subset(sanka, jackknife_run=="all"), subset(rt, jackknife_run=="all"), by="bin_cM")
    M <- merge(M, d, by.x="bin_cM", by.y="V1")
    M <- merge(M, subset(moorjani, jackknife_run=="all"), by="bin_cM")
    M <- data.frame(x=M[,1], d=M$V2, rt=M$exp_fit.y, asSanka=M$exp_fit.x, moorjani=M$average_cov)
    #
    M <- subset(M, x>=0.02 & x<=1.0)
    #
    N <- c(); for (j in 2:ncol(M)) {
      df <- data.frame(x=M[,1], y=M[,j])
      # Levenberg-Marquardt nls algo
      model <- minpack.lm::nlsLM(y ~ a*exp(-1*t*x)+c, data = df, start = list("a" = 0.01, "t"=10, "c"=1e-5))
      coef <- summary(model)$coefficients[2,1]*100
      coef <- ifelse(summary(model)$coefficients[1,1] < 0, NA, coef)
      N <- c(N, summary(model)$coefficients[2,1]*100)
    }
    #
    V <- rbind(V, data.frame(true.age = par[par$V1=="T.admix",2],
                             true.alpha = par[par$V1=="a.admix.EuA__0_from_Nea__0",2],
                             as.data.frame(t(N))))
  }
  V <- as.data.frame(V)
  colnames(V) <- c("true.age", "true.alpha", "computed", "rt.rt", "rt.likesanka", "moorjani")
  rownames(V) <- NULL
  
  V$true.age <- V$true.age / tg
  
  v <- melt(V, id.vars = c("true.age", "true.alpha"))
  v <- subset(v, true.alpha > 0 & variable != "moorjani")
  v$implementation <- sapply(as.character(v$variable), function(x) strsplit(x, ".", fixed = T)[[1]][1])
  v$implementation <- plyr::mapvalues(v$implementation, c("computed", "rt"), c("computed", "This study"))
  v$variable <- plyr::mapvalues(as.character(v$variable), c("computed", "rt.rt", "rt.likesanka"), c("computed", "This study", "Original"))
  lims <- range(c(v$true.age, v$value))
  v$variable <- factor(v$variable, levels = c("Original", "This study", "computed"))
  g1=ggplot(v) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", alpha = .7) +
    geom_point(aes(true.age, value, col = variable, shape = implementation, size = variable), alpha = 1) +
    facet_wrap(~true.alpha) +
    labs(x = "Simulated admixture age (ga)", 
         y = "Estimated admixture age (ga)", 
         caption = paste0("Ascertainment 0\nAdmixture proportions: 2.5% (left) or 5% (right)"),
         color = "Parameters",
         size = "Parameters",
         shape = "Implementation") +
    theme_light(base_size = 15, base_family = "Poppins") +
    theme(panel.grid = element_blank()) +
    coord_equal(xlim = lims) +
    scale_color_manual(values = c("goldenrod2", "black", "dodgerblue4")) +
    scale_shape_manual(values = c(3, 20)) +
    scale_size_manual(values = c(3, 5, 6)) +
    ggtitle("(A) With archaic admixture")
    
  v <- melt(V, id.vars = c("computed", "true.alpha"))
  v <- subset(v, true.alpha == 0 & variable != "moorjani" & variable != "true.age")
  v$implementation <- sapply(as.character(v$variable), function(x) strsplit(x, ".", fixed = T)[[1]][1])
  v$variable <- plyr::mapvalues(as.character(v$variable), c("rt.rt", "rt.likesanka"), c("This study", "Original"))
  lims <- range(c(v$computed, v$value))
  v$variable <- factor(v$variable, levels = c("Original", "This study", "computed"))
  g2=ggplot(v) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", alpha = .7) +
    geom_point(aes(computed, value, col = variable, size = variable), alpha = 1) +
    labs(x = "Estimated with computed (ga)", 
         y = "Estimated with our algo (ga)", 
         caption = paste0("Ascertainment 0\nNo admixture"),
         color = "Parameters",
         size = "Parameters",
         shape = "Implementation") +
    theme_light(base_size = 15, base_family = "Poppins") +
    theme(panel.grid = element_blank()) +
    coord_equal(xlim = lims) +
    scale_color_manual(values = c("goldenrod2", "dodgerblue4")) +
    scale_shape_manual(values = c(3, 1)) +
    scale_size_manual(values = c(3, 6)) +
    ggtitle("(B) No admixture")
  
  G <- cowplot::plot_grid(g1, g2, ncol=1)
  ggsave(G, filename = paste0(DirOut,"/further/LD_implementation.png"), scale = 1.2)
  
  sink(paste0(DirOut,"/further/LD_implementation.txt"))
  cat("*** With admixture\n")
  cat("RMSE computed vs true.age\n")
  print(with(subset(V, true.alpha > 0), rmse(computed, true.age)))
  cat("RMSE RT_Original vs true.age\n")
  print(with(subset(V, true.alpha > 0), rmse(rt.likesanka, true.age)))
  cat("RMSE RT_RT vs true.age\n")
  print(with(subset(V, true.alpha > 0), rmse(rt.rt, true.age)))
  cat("R2 computed vs RT_RT\n")
  print(with(subset(V, true.alpha > 0), cor.test(computed, rt.rt, method = "pearson")))
  
  cat("*** No admixture\n")
  cat("R2 computed vs RT_Orignal\n")
  print(with(subset(V, true.alpha == 0), cor.test(computed, rt.likesanka, method = "pearson")))
  cat("R2 computed vs RT_RT\n")
  print(with(subset(V, true.alpha == 0), cor.test(computed, rt.rt, method = "pearson")))
  
  cat("*** No admixture - computed only\n")
  cat("Nsim\n")
  print(length(subset(V, true.alpha == 0)$computed))
  print(summary(subset(V, true.alpha == 0)$computed))
  cat("SD\n")
  print(sd(subset(V, true.alpha == 0)$computed))
  sink()
  
}
  
#__________________SupMat ::: Published Models___________________# 

if (F) {
  dir.create(paste0(DirOut,"/plots"), showWarnings = FALSE)
  Models <- c("kamm19", "ragsdale19", "jacobs19", "gower21", "iasi21", "durvasula20", "fu14", "yang12", "skov20", "moorjani16", "schaefer21")
  #
  for (Model in Models) {
    cat(Model," ")
    #___ N best runs
    DATA <- abc.import(paste0("Final.Blake/data/20x30Mbp.50sim/",Model), 
                       do.AFS.2D = F, 
                       read.params = T, 
                       generation.time = tg, 
                       Sprime.min.match.rate = Sprime.min.match.rate)
    R <- abc.analyze(obs.data, 
                     DATA, 
                     parfile = spar, 
                     N = N, 
                     dispersion.fun = "mad", 
                     psmc..xlim = psmc..xlim)
    DF <- data.frame(x=R$ACCEPTED$name, mu=DATA$PAR[R$ACCEPTED$index,], lab=paste0(round(1e8*DATA$PAR[R$ACCEPTED$index,],3),"e-8"), stringsAsFactors = F)
    DF <- DF[order(DF$mu),]
    multidata_plots(data_paths = DF[,1], 
                    labels = DF[,3],                    
                    output = paste0(DirOut,"/plots/20x30Mbp.",Model,".png"),
                    title = "",
                    do_tile = FALSE,
                    plot.legend = TRUE,
                    color.lims = c(0.1, 0.9), 
                    log10.crf = c(T, T),
                    xylims.crf = c("xmin"=5e3, "xmax"=1e6, "ymin"=0.1e-2, "ymax"=10e-2),
                    log10.sprime = c(T, F),
                    xylims.sprime = c("xmin"=1e4, "xmax"=3e6, "ymin"=0., "ymax"=1.),
                    ylim.afs.ceu = c(0, .35),
                    ylim.afs.yri = c(0, .45),
                    ylim.dcfs = c(0, .45),
                    xlim.D = c(-0.02, .20),
                    squeeze = squish,
                    show.RMSE = TRUE,
                    plot.obs.TLD = TRUE )
    #___ 1.2e-8
    multidata_plots(data_paths = paste0("Final.Blake/data/20x30Mbp.1_2e-8/",Model,"/",Model,"_",1), 
                    labels = paste0("M",1),                    
                    output = paste0(DirOut,"/plots/20x30Mbp.",Model,".1_2e-8.png"),
                    title = "",
                    do_tile = FALSE,
                    plot.legend = FALSE,
                    color.lims = c(0.7, 0.7), 
                    log10.crf = c(T, T),
                    xylims.crf = c("xmin"=5e3, "xmax"=1e6, "ymin"=0.1e-2, "ymax"=10e-2),
                    log10.sprime = c(T, F),
                    xylims.sprime = c("xmin"=1e4, "xmax"=3e6, "ymin"=0., "ymax"=1.),
                    ylim.afs.ceu = c(0, .35),
                    ylim.afs.yri = c(0, .45),
                    ylim.dcfs = c(0, .45),
                    xlim.D = c(-0.02, .20),
                    squeeze = squish,
                    show.RMSE = TRUE,
                    plot.obs.TLD = TRUE )
  }
  
  #___ All
  X <- c(paste0("Final.Blake/data/20x30Mbp.1_2e-8/",Models,"/",Models,"_",1), "Final.Blake/data/20x30Mbp/1014230")
  Y <- c(Models, "THIS STUDY")
  multidata_plots(data_paths = X, 
                  labels = Y,                    
                  output = paste0(DirOut,"/plots/20x30Mbp.All.png"),
                  title = "",
                  do_tile = FALSE,
                  plot.legend = TRUE,
                  color.lims = c(0.1, 0.9), 
                  log10.crf = c(F, F),
                  xylims.sprime = c("xmin"=5e3, "xmax"=800e3, "ymin"=0, "ymax"=1), 
                  log10.sprime = c(F, F),
                  xylims.crf = c("xmin"=1e4, "xmax"=400e3, "ymin"=.001, "ymax"=.035),
                  ylim.afs.ceu = c(0, .35),
                  ylim.afs.yri = c(0, .45),
                  ylim.dcfs = c(0, .45),
                  xlim.D = c(-0.02, .20),
                  squeeze = squish,
                  plot.obs.TLD = TRUE )
  
}

#__________________Fst trajectories___________________# 

if (F) {
  get_fst <- function(dir, label) {
    ls <- list.files(dir, pattern = "stats", full.names = T)
    fst <- sapply(ls, function(x) as.numeric(strsplit(scan(x, what = "character", sep = "\n", quiet=T)[6], " ")[[1]][3]))
    ages <- sapply(ls, function(x) rev(strsplit(x, "/")[[1]])[1] )
    ages <- as.numeric(gsub(".stats", "", ages))
    fst <- data.frame(type=label, age=ages, Fst=fst)
    return(fst)
  }
  run <- "1014230"
  #___
  pars <- read.table(paste0("Final.Blake/data/20x30Mbp/",run,".par"), header = F)
  pars <- subset(pars, grepl("^T", V1, perl=T))
  pars <- subset(pars, V1!="T.0" & V1!="T.mis5")
  #___
  Fst <- rbind(get_fst("Final.Blake/further/Fst_traj_AfW_AfE", "Between Ma-Mb"),
               get_fst("Final.Blake/further/Fst_traj_AfE_EuA", "Between Mb-Mc"),
               get_fst("Final.Blake/further/Fst_traj_AfW0_AfW9", "Within Ma"),
               get_fst("Final.Blake/further/Fst_traj_AfE0_AfE9", "Within Mb"))
  Fst[Fst$Fst<0,"Fst"] <- 0
  g= ggplot(Fst) + 
    geom_vline(data = pars, aes(xintercept = V2, linetype = V1), size = .7, alpha = .6) +
    geom_line(aes(age, Fst, color = type), size = 1.0) +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    theme_light(base_size = 15, base_family = "Poppins") +
    scale_color_manual(values = viridis(4, end = .8), 
                          labels = c(expression(paste("between ",M[A]-M[B])), 
                                     expression(paste("between ",M[B]-M[C])),
                                     expression(paste("within ",M[A])), 
                                     expression(paste("within ",M[B]))  )) +
    scale_y_sqrt(limits = c(0, 1.), 
                 breaks = c(0.01, 0.05, 0.1, 0.25, .5, .75, 1.)) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted"),
                       labels = c(expression(paste("Founding of ",M[N])),
                                  expression(paste("Founding of ",M[C])),
                                  expression(paste("Founding of ",M[B])) ) ) +
    labs(x = expression(paste("Sampling age (ya) (lo",g[10],"-scale)")),
         y = expression(paste(F[ST])),
         col = "Times",
         linetype = "Pairs") +
    theme(legend.position = "right",
          panel.grid.minor.y = element_blank()); g
  ggsave(g, filename = paste0(DirOut,"/Fst.png"), width = 9.5, height = 6.5)
  #
  sink(paste0(DirOut,"/Fst.txt"))
  sumup("Fst : Ma-Mb", matrix(subset(Fst, type=="Between Ma-Mb")$Fst, ncol=1), matrix(NA, ncol=1))
  sumup("Fst : Mb-Mc", matrix(subset(Fst, type=="Between Mb-Mc")$Fst, ncol=1), matrix(NA, ncol=1))
  sumup("Fst : Ma", matrix(subset(Fst, type=="Within Ma")$Fst, ncol=1), matrix(NA, ncol=1))
  sumup("Fst : Mb", matrix(subset(Fst, type=="Within Mb")$Fst, ncol=1), matrix(NA, ncol=1))
  cat("\n\nFst at 0 ya ::: Between / Within_A ::: ",subset(Fst, type=="Between Ma-Mb" & age==1)$Fst / subset(Fst, type=="Within Ma" & age==1)$Fst)
  cat("\n\nFst at 0 ya ::: Between / Within_B ::: ",subset(Fst, type=="Between Ma-Mb" & age==1)$Fst / subset(Fst, type=="Within Mb" & age==1)$Fst)
  sink()
}

#__________________PSMC ::: 20x30Mbp___________________# 

if (F) {
  AccNames <- unname(unlist(read.table(paste0(DirOut,"/FINAL.accepted"), header = F)))
  fun.summary <- "min"
  PubModels <- c("durvasula20", "fu14", "gower21", "iasi21" , "jacobs19", "kamm19", "moorjani16" , "ragsdale19", "schaefer21", "skov20" , "yang12")
  DirRT <- "Final.Blake/data/20x30Mbp.50sim/1M.accepted"
  #================> var mutrate ::: ScaledL2_CrossCor_2PSMC.spar
  DirPublished <- "Final.Blake/data/20x30Mbp.50sim/"
  Published <- paste0(DirPublished,PubModels)
  names(Published) <- gsub(DirPublished, "", Published, fixed = TRUE)
  Res3 <- compare(DirRT, AccNames, Published, N, spar, fun.summary,
                  comparison.spar = "param_files/ScaledL2_CrossCor_2PSMC.spar", 
                  outprefix = paste0(DirOut,"/temp"),
                  tg, Sprime.min.match.rate, psmc..xlim, do.pdf = TRUE )
  
  # Extract the fittest PSMC for each model
  PSMC <- c()
  for (model in c(PubModels)) { # c(PubModels, "z. This study")
    indir <- Res3$accepted.runs[[spar]][[model]]
    full <- indir
    targets <- unname(sapply(indir, function(x) rev(strsplit(x, "/", fixed=T)[[1]])[1]))
    indir <- paste0(rev(rev(strsplit(indir[1], "/", fixed=T)[[1]])[-1]), collapse="/")
    DATA <- abc.import(indir,
                       target = full,
                       do.AFS.2D = F, 
                       read.params = F, 
                       generation.time = tg, 
                       Sprime.min.match.rate = Sprime.min.match.rate)
    R <- abc.analyze(obs.data, 
                     DATA, 
                     parfile = "param_files/PSMC.CEU_YRI.spar", 
                     N = 1, 
                     dispersion.fun = "mad", 
                     psmc..xlim = psmc..xlim)
    acc <- R$ACCEPTED$index
    psmc <- rbind(
      data.frame(x=DATA$SS$PSMC.CEU[[1]][acc,], 
                 xend=c(DATA$SS$PSMC.CEU[[1]][acc,][-1], Inf),
                 y=DATA$SS$PSMC.CEU[[2]][acc,], 
                 pop="CEU"),
      data.frame(x=DATA$SS$PSMC.YRI[[1]][acc,], 
                 xend=c(DATA$SS$PSMC.YRI[[1]][acc,][-1], Inf),
                 y=DATA$SS$PSMC.YRI[[2]][acc,], 
                 pop="YRI")
    )
    psmc$model <- model
    PSMC <- rbind(PSMC, psmc)
  }
  
  # 3 runs for us:
  model <- "z. This study"
  indir <- Res3$accepted.runs[[spar]][[model]]
  indir <- paste0(rev(rev(strsplit(indir[1], "/", fixed=T)[[1]])[-1]), collapse="/")
  full <- paste0(indir,"/",c("152019", "35483", "344388"))
  DATA <- abc.import(indir,
                       target = full,
                       do.AFS.2D = F, 
                       read.params = F, 
                       generation.time = tg, 
                       Sprime.min.match.rate = Sprime.min.match.rate)
  R <- abc.analyze(obs.data, 
                   DATA, 
                   parfile = "param_files/PSMC.CEU_YRI.spar", 
                   N = 3, 
                   dispersion.fun = "mad", 
                   psmc..xlim = psmc..xlim)
  metacc <- R$ACCEPTED$index
  psmc <- c()
  for (acc in metacc) psmc <- rbind(psmc,
    data.frame(x=DATA$SS$PSMC.CEU[[1]][acc,], 
               xend=c(DATA$SS$PSMC.CEU[[1]][acc,][-1], Inf),
               y=DATA$SS$PSMC.CEU[[2]][acc,], 
               pop="CEU",
               model = paste0(model,".",acc)),
    data.frame(x=DATA$SS$PSMC.YRI[[1]][acc,], 
               xend=c(DATA$SS$PSMC.YRI[[1]][acc,][-1], Inf),
               y=DATA$SS$PSMC.YRI[[2]][acc,], 
               pop="YRI",
               model = paste0(model,".",acc))
  )
  PSMC <- rbind(PSMC, psmc)
  
  # For observed PSMC
  PSMC.obs <- rbind(
    data.frame(x=obs.data$SS$PSMC.CEU[[1]][1,], 
               xend=c(obs.data$SS$PSMC.CEU[[1]][1,][-1], Inf),
               y=obs.data$SS$PSMC.CEU[[2]][1,], 
               pop="CEU", 
               model="obs"),
    data.frame(x=obs.data$SS$PSMC.YRI[[1]][1,], 
               xend=c(obs.data$SS$PSMC.YRI[[1]][1,][-1], Inf),
               y=obs.data$SS$PSMC.YRI[[2]][1,], 
               pop="YRI", 
               model="obs")
  )
  
  # Extractions
  PopsToPlot <- c("CEU", "YRI")
  PSMC.2 <- subset(PSMC, 
                   (model %in% c("schaefer21", "moorjani16", "skov20") | grepl("z. This study", model)) & pop%in%PopsToPlot)
  
  # Plot
  g=ggplot(PSMC.2) +
    geom_rect(data=subset(PSMC.obs, pop%in%PopsToPlot), 
              aes(xmin = x, xmax = xend, 
                  ymin = 0, ymax = y), alpha = 1.0, fill = "gray50", col = "gray50") +
    #geom_step(data=subset(PSMC.obs, pop%in%PopsToPlot), aes(x, y), size = 0.1, alpha = 1.0) +
    geom_step(data = subset(PSMC.2, !grepl("z. This study", model)), 
              aes(x, y, col=model), size = 1.2, alpha = 1.0) +
    geom_step(data = subset(PSMC.2, grepl("z. This study", model)), 
              aes(x, y, col=model), size = 1.2, alpha = 1.0) +
    facet_wrap(~pop) +
    theme_classic() +
    ylim(c(0, 1e5)) +
    scale_color_manual(values = c("#337893",
                                  "#5acdda",
                                  "#c8bb5f",
                                  "#c43b60cc", "#c4613bcc", "#c43b9fcc")) +
    labs(x = "",
         y = "") +
    scale_x_log10(#limits = c(3e3, 1.75e7),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b"); g
  ggsave(g, filename = paste0(DirOut,"/20x30Mbp.COMPARISON.PSMC3.pdf"), width = 12, height = 6)
  
 }

#__________________Parameters___________________#

if (F) {
  DATA <- abc.import("Final.Blake/data/20x30Mbp",
                     do.AFS.2D = F, 
                     read.params = T, 
                     generation.time = tg, 
                     Sprime.min.match.rate = Sprime.min.match.rate)
  Par <- DATA$PAR
  rownames(Par) <- sapply(DATA$NAME, function(x) rev(strsplit(x,"/")[[1]])[1])
  write.table(t(Par), "Final.Blake/Parameters.Accepted.txt", quote=F, sep="\t")
}

#__________________Published Models ::: Stats___________________#

if (F) {
  get_data <- function(Model) {
    DATA <- abc.import(paste0("Final.Blake/data/20x30Mbp.50sim/",Model), 
                       do.AFS.2D = F, 
                       read.params = T, 
                       generation.time = tg, 
                       Sprime.min.match.rate = Sprime.min.match.rate)
    R <- abc.analyze(obs.data, 
                     DATA, 
                     parfile = spar, 
                     N = N, 
                     dispersion.fun = "mad", 
                     psmc..xlim = psmc..xlim)
    IDX <- R$ACCEPTED$index
    SS <- DATA$SS
    MUT <- DATA$PAR
    colnames(MUT) <- "MUT"
    return(list("SS"=SS, "IDX"=IDX, "MUT"=MUT)) 
  }
  
  Models <- c("kamm19", "ragsdale19", "jacobs19", "gower21", "iasi21", "durvasula20", "fu14", "yang12", "skov20", "moorjani16", "schaefer21")
  
  L <- get_data("schaefer21")
  
  # ~MutRate
  X <- data.frame(L$MUT, L$SS$PI); X[order(X$MUT),]
  print("20 retained:")
  X <- data.frame(L$MUT, L$SS$PI); X=X[L$IDX,,drop=F]; X[order(X$MUT),]
  
  # RPE (%)
  X <- L$SS$Fst[,1]; X=X[L$IDX]; x=abs(X-obs.data$SS$Fst[,1])/obs.data$SS$Fst[,1]; mean(x)
  X <- L$SS$DCFS[,1]; X=X[L$IDX]; x=abs(X-obs.data$SS$DCFS[,1])/obs.data$SS$DCFS[,1]; mean(x)
  X <- L$SS$AFS.YRI[,1]; X=X[L$IDX]; x=abs(X-obs.data$SS$AFS.YRI[,1])/obs.data$SS$AFS.YRI[,1]; mean(x)
  X <- L$SS$D[,1]; X=X[L$IDX]; x=abs(X-obs.data$SS$D[,1])/obs.data$SS$D[,1]; mean(x)
  X <- L$SS$CRF[,"length.mean"]; X=X[L$IDX]; x=abs(X-obs.data$SS$CRF[,"length.mean"])/obs.data$SS$CRF[,"length.mean"]; mean(x)
  X <- L$SS$CRF[,"alpha.mean"]; X=X[L$IDX]; x=abs(X-obs.data$SS$CRF[,"alpha.mean"])/obs.data$SS$CRF[,"alpha.mean"]; mean(x)
  X <- L$SS$SPRIME[,"match.mean"]; X=X[L$IDX]; x=abs(X-obs.data$SS$SPRIME[,"match.mean"])/obs.data$SS$SPRIME[,"match.mean"]; mean(x)
  # Fold-diff
  X <- L$SS$D[,1]; X=X[L$IDX]; x=abs(X/obs.data$SS$D[,1]); c(mean(x),1/mean(x))
  X <- L$SS$Fst[,1]; X=X[L$IDX]; x=abs(X/obs.data$SS$Fst[,1]); c(mean(x),1/mean(x))
  X <- L$SS$CRF[,"alpha.mean"]; X=X[L$IDX]; x=abs(X/obs.data$SS$CRF[,"alpha.mean"]); c(mean(x),1/mean(x))
  X <- L$SS$SPRIME[,"match.mean"]; X=X[L$IDX]; x=abs(X/obs.data$SS$SPRIME[,"match.mean"]); c(mean(x),1/mean(x))
  # Descriptive
  X <- L$SS$D[,1]; X=X[L$IDX]; print(sort(X))
  X <- L$SS$LD[,"mean"]; X=X[L$IDX]; mean(X); range(X)
  X <- L$SS$CRF[,"alpha.mean"]*100; X=X[L$IDX]; print(sort(X))

}

#___