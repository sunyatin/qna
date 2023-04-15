
# Version Control:
# 0     10.05.21   Creation
# 1     17.05.21   added abc.par()
# 2     19.05.21   changed plot outline (introduced CRF + changed plot order) + flexible scientific notation in gg.scatter()
# 2.1   22.06.21   added `show.legend` option + labels
# 2.2   28.10.21   trick if ggdist not installed (eg issue on elephant01a)
# 2.3   20.12.21   implemented KS test to tell if posterior distrib is uniform or not
# 2.4   04.05.22   added possibility to scales::squish, or not (arg `squeeze`) the points/lines, etc.
#                  if squeeze=F, we report the number of points effectively plotted as `n=n_points` for plot.scatter() [if squeeze=T, we do not print as we plot all points by default]
#                  note for gg.ld.errorbars : we still squish, as this is the only way to represent the error bars when they go beyong the xlim
# 3     01.07.22   changed gg.scatter to really plot IQR95 for S', CRF instead of an artificially estimated CI95
# 4     12.12.22   added a new layout "B" without CRF but with PSMC.YRI
# 4.1   19.12.22   switched args or RMSE from `RMSE <- function(obs, sim)` to `RMSE_plot <- function(sim, obs)`
# 4.1.1 28.03.23   added option to plot vertical line at 50 kya for T_LD


# LAYOUTS:
# A | pi, AFS:CEU, AFS:YRI, PSMC:CEU, Fst~D, DCFS, S', CRF, LD
# B | pi, AFS:CEU, AFS:YRI, PSMC:CEU, PSMC:YRI, Fst~D, DCFS, S', LD

suppressMessages(require(cowplot))
suppressMessages(require(reshape2))
suppressMessages(require(viridis))
suppressMessages(require(ggplot2))
suppressMessages(require(ggridges))
if ("ggdist" %in% rownames(installed.packages())) {
	suppressMessages(require(ggdist))
	general_theme <- theme_ggdist()
	general_theme <- theme_light()
} else {
	general_theme <- theme_bw()
}

#################################################################################
#################################################################################
#############################          Plotting         #########################
#################################################################################
#################################################################################

split <- function(string)
{
  if (is.na(string)) return(list("s"=NA, "c"=NA))
  string <- strsplit(string, ":", fixed = TRUE)[[1]]
  col <- if (grepl("^\\d", string[2])) as.integer(string[2]) else string[2]
  return(list("s"=string[1], "c"=col))
}

get.sm <- function(data, 
                   accepted.indices, 
                   stat.x,
                   stat.y, 
                   se.x=NA, 
                   se.y=NA, 
                   iqr.x.lb=NA, 
                   iqr.x.ub=NA, 
                   iqr.y.lb=NA, 
                   iqr.y.ub=NA, 
                   labels = NULL)
{
  x <- data$SS[[stat.x$s]][accepted.indices,,drop=FALSE]
  y <- data$SS[[stat.y$s]][accepted.indices,,drop=FALSE]
  sm <- data.frame(id = 1:length(accepted.indices),
                   x = x[,stat.x$c],
                   y = y[,stat.y$c])
  sm$se.x <- if (is.na(se.x$s)) NA else x[,se.x$c]
  sm$se.y <- if (is.na(se.y$s)) NA else y[,se.y$c]
  sm$iqr.x.lb <- if (is.na(iqr.x.lb$s)) NA else x[,iqr.x.lb$c]
  sm$iqr.x.ub <- if (is.na(iqr.x.ub$s)) NA else y[,iqr.x.ub$c]
  sm$iqr.y.lb <- if (is.na(iqr.y.lb$s)) NA else x[,iqr.y.lb$c]
  sm$iqr.y.ub <- if (is.na(iqr.y.ub$s)) NA else y[,iqr.y.ub$c]
  if (!is.null(labels)) sm$labels <- labels
  return(sm)
}

empty_plot <- function(x) return(ggplot() + theme_void() + geom_text(aes(0, 0, label="."), col = "white") + ggtitle(x))

RMSE_plot <- function(sim, obs)
{
  RMSE <- sweep(sim, 2, obs, "-")**2
  RMSE <- sqrt(rowMeans(RMSE))
  RMSE <- ifelse(length(RMSE)==1, paste0(round(RMSE*100,1),"%"), paste0(round(mean(RMSE*100),1),"±",round(sd(RMSE*100),1),"%"))
  RMSE <- bquote(atop(" ", atop(textstyle(RMSE == .(RMSE) ~ .("   ")))))
  return(RMSE)
}

gg.sfs <- function(obs.data,
                   sim.data,
                   accepted,
                   stat,
                   title,
                   ylim,
                   scg = scale_color_viridis(discrete = FALSE, option = "B", end = .9),
                   line.size = 1.2,
                   tile = FALSE,
                   show.RMSE = TRUE,
                   axis.labs = c("", ""),
                   oob = oob )
{
  oob <- ifelse(oob, scales::squish, scales::censor)
  if (is_empty(obs.data, sim.data, stat)) return(empty_plot(title))
  s <- sim.data$SS[[stat]][accepted[,1],,drop=FALSE]
  o <- obs.data$SS[[stat]]
  sm <- reshape2::melt(s)
  om <- reshape2::melt(o)
  #
  g <- ggplot(sm)
  if (tile) {
    g <- g + stat_lineribbon(aes(x=Var2, y=value), .width = c(0.5, 0.8, 0.95, 1.0), size=.7) +
      scale_fill_brewer(palette="YlOrRd")
      #scg
  } else {
    g <- g + geom_line(aes(x=Var2, y=value, col=Var1, group=Var1), size = line.size, alpha = .7) + scg
  }
  g <- g +
    geom_point(data = om, aes(x=Var2, y=value), size = 3, col = "white") +
    geom_point(data = om, aes(x=Var2, y=value), size = 2) +
    # format
    #scale_x_continuous(breaks = seq(min(sm$Var2), max(sm$Var2)), oob = oob) +
    #scale_y_continuous(limits = c(ylim[1], ylim[2]), oob = oob, labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(ylim[1], ylim[2])) +
    general_theme +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = "none") +
    ggtitle(title)
  g <- g + xlab(axis.labs[1]) + ylab(axis.labs[2])
  if (show.RMSE) g <- g + annotate(geom="text", x=Inf, y=Inf, vjust=1, hjust=1, label=RMSE_plot(s, o))
  return(g)
}

gg.scatter <- function(obs.data,
                       sim.data,
                       accepted,
                       stat.x,
                       stat.y,
                       title,
                       xlim,
                       ylim,
                       xlab,
                       ylab,
                       se.x = NA,
                       se.y = NA,
                       iqr.x = c(NA, NA),
                       iqr.y = c(NA, NA),
                       error.width.factor = 0, #.0005,
                       error.height.factor = 0,# .05,
                       scg = scale_color_viridis(discrete = FALSE, option = "B", end = .9),
                       log10 = c(FALSE, FALSE),
                       lab.percent = c(FALSE, FALSE),
                       point.size = 1.2,
                       line.size = 1.2,
                       tile = FALSE,
                       show.legend = FALSE,
                       axis.labs = c("", ""),
                       oob = oob )
{
  # Note: by default, we print the IQR95 except if those were not imported and their values is NA
  # in which case, we plot the CI95 as mean +/- 1.96*SE
  is_squish = oob
  oob <- ifelse(oob, scales::squish, scales::censor)
  # format of stat.x, e.g. D:D or PI:CEU
  stat.x <- split(stat.x)
  stat.y <- split(stat.y)
  se.x <- split(se.x)
  se.y <- split(se.y)
  iqr.x.lb <- split(iqr.x[1])
  iqr.x.ub <- split(iqr.x[2])
  iqr.y.lb <- split(iqr.y[1])
  iqr.y.ub <- split(iqr.y[2])
  #
  if (is_empty(obs.data, sim.data, stat.x$s) || is_empty(obs.data, sim.data, stat.y$s)) return(empty_plot(title))
  sm <- get.sm(sim.data, accepted[,1], stat.x, stat.y, se.x, se.y, iqr.x.lb, iqr.x.ub, iqr.y.lb, iqr.y.ub, labels = accepted[,2])
  om <- get.sm(obs.data, 1, stat.x, stat.y, se.x, se.y, iqr.x.lb, iqr.x.ub, iqr.y.lb, iqr.y.ub)
  # errorbar parameters
  eh <- ifelse(log10[2], 10**diff(log10(ylim)), diff(ylim))*error.height.factor
  ew <- ifelse(log10[1], 10**diff(log10(xlim)), diff(xlim))*error.width.factor
  if (!is.na(iqr.x[1])) {
    sm$x.lb <- sm$iqr.x.lb
    sm$x.ub <- sm$iqr.x.ub
    #
    om$x.lb <- om$iqr.x.lb
    om$x.ub <- om$iqr.x.ub
  } else {
    sm$x.lb <- sm$x-1.96*sm$se.x
    sm$x.ub <- sm$x+1.96*sm$se.x
    #
    om$x.lb <- om$x-1.96*om$se.x
    om$x.ub <- om$x+1.96*om$se.x
  }
  if (!is.na(iqr.y[1])) {
    sm$y.lb <- sm$iqr.y.lb
    sm$y.ub <- sm$iqr.y.ub
    #
    om$y.lb <- om$iqr.y.lb
    om$y.ub <- om$iqr.y.ub
  } else {
    sm$y.lb <- sm$y-1.96*sm$se.y
    sm$y.ub <- sm$y+1.96*sm$se.y
    #
    om$y.lb <- om$y-1.96*om$se.y
    om$y.ub <- om$y+1.96*om$se.y
  }
  # errorbar cutting, if not scales::squish
  if (!is_squish) {
    sm[!is.na(sm$x.lb) & sm$x.lb<xlim[1],"x.lb"] <- xlim[1]
    sm[!is.na(sm$x.ub) & sm$x.ub>xlim[2],"x.ub"] <- xlim[2]
    sm[!is.na(sm$y.lb) & sm$y.lb<ylim[1],"y.lb"] <- ylim[1]
    sm[!is.na(sm$y.ub) & sm$y.ub>ylim[2],"y.ub"] <- ylim[2]
    #
    om[!is.na(om$x.lb) & om$x.lb<xlim[1],"x.lb"] <- xlim[1]
    om[!is.na(om$x.ub) & om$x.ub>xlim[2],"x.ub"] <- xlim[2]
    om[!is.na(om$y.lb) & om$y.lb<ylim[1],"y.lb"] <- ylim[1]
    om[!is.na(om$y.ub) & om$y.ub>ylim[2],"y.ub"] <- ylim[2]
  }
  # target + rugs
  g <- ggplot(sm)
  # errorbars
  if (!tile) {
    if (!is.na(iqr.x.lb)||!is.na(se.x)) g <- g + geom_errorbarh(aes(y=y, xmin=x.lb, xmax=x.ub, col=id), alpha = .7, size = line.size, height = eh)
    if (!is.na(iqr.y.lb)||!is.na(se.y)) g <- g + geom_errorbar(aes(x=x, ymin=y.lb, ymax=y.ub, col=id), alpha = .7, size = line.size, width = ew)
  }
  # points
  if (tile) {
    g <- g + #stat_density_2d_filled(aes(x=x, y=y), alpha=.7)
      stat_density_2d_filled(aes(x=x, y=y), contour_var="count", alpha=.9) #+
      #scale_fill_viridis(discrete=T, option="A", direction=-1)
      #geom_hex(aes(x=x, y=y), bins=20) + scale_fill_gradient(low="pink1", high="darkred")

  } else {
    g <- g +
      geom_point(aes(x=x, y=y, col=id, group=id), size = point.size) +
      geom_rug(aes(x=x, y=y), col=rgb(0, .5, 0, alpha=.7))
  }
  # axis-x
  l <- if (lab.percent[1]) scales::percent_format(accuracy = 1) else if (log10[1]) scales::trans_format("log10", scales::math_format(10^.x)) else if (any(om$x > 1000)) scales::unit_format(unit = "k", scale = 1e-3) else waiver()
  g <- if (log10[1]) {
    g + scale_x_log10(limits = xlim, oob = oob, labels = l) +
      annotation_logticks(sides = "b")
  } else {
    g + scale_x_continuous(limits = xlim, oob = oob, labels = l)
  }
  # axis-y
  l <- if (lab.percent[2]) scales::percent_format(accuracy = 1) else if (log10[2]) scales::trans_format("log10", scales::math_format(10^.x)) else if (any(om$y > 1000)) scales::unit_format(unit = "k", scale = 1e-3) else waiver()
  g <- if (log10[2]) {
    g + scale_y_log10(limits = ylim, oob = oob, labels = l) +
      annotation_logticks(sides = "l")
  } else {
    largeScale <-
    g + scale_y_continuous(limits = ylim, oob = oob, labels = l)
  }
  # observed
  g <- g +
    geom_point(data = om, aes(x=x, y=y), size = 4, col = "white") +
    geom_point(data = om, aes(x=x, y=y), size = 2) +
    geom_vline(data = om, aes(xintercept = x), linetype = "dotted", alpha = .6) +
    geom_hline(data = om, aes(yintercept = y), linetype = "dotted", alpha = .6)
  # errorbars for observed
  if (!is.na(iqr.x.lb)||!is.na(se.x)) g <- g + geom_errorbarh(data = om, aes(y=y, xmin=x.lb, xmax=x.ub), alpha = 1., size = line.size, height = eh)
  if (!is.na(iqr.y.lb)||!is.na(se.y)) g <- g + geom_errorbar(data = om, aes(x=x, ymin=y.lb, ymax=y.ub), alpha = 1., size = line.size, width = ew)
  # theme
  g <- g +
    general_theme +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    ggtitle(title)
  g <- g + xlab(axis.labs[1]) + ylab(axis.labs[2])
  if (show.legend==FALSE) {
    g <- g + scg + theme(legend.position = "none")
  } else {
    scg1 <- scg
    scg1$breaks <- 1:nlevels(sm$labels)
    scg1$labels <- sm$labels
    if (length(sm$labels)==2) scg1$breaks <- c(1., 2.)
    g <- g + scg1 + labs(color="")
    g <- g + theme(legend.position = "right", legend.key.height= unit(1.2, 'cm'))
  }
  # number of effective points
  if (!is_squish) {
    n_points <- nrow(subset(sm, x>=xlim[1] & x<=xlim[2] & y>=ylim[1] & y<=ylim[2]))
    #g <- g + geom_text(x = xlim[1], y = ylim[2]*0.9, label = paste0("n=",n_points), vjust = 0, hjust = 0, size = 5) #, family = "serif")
    g <- g + annotate("text", x = Inf, y =Inf, label = paste0("n=",n_points), vjust = 1, hjust = 1, size = 4.5) #, family = "serif")
    #g <- g + annotate("text", x = -Inf, y = Inf, label = paste0("n=",n_points), vjust = 1, hjust = 0)
  }
  return(g)
}


gg.psmc <- function(obs.data,
                     sim.data,
                     accepted,
                     stat,
                     se.x,
                     se.y,
                     title,
                     xlim,
                     ylim,
                     xlab,
                     ylab,
                     scg = scale_color_viridis(discrete = FALSE, option = "B", end = .9),
                     line.size = 1.2,
                     tile = FALSE,
                     axis.labs = c("", ""),
                     oob = oob )
{
  oob <- ifelse(oob, scales::squish, scales::censor)
  if (is_empty(obs.data, sim.data, stat)) return(empty_plot(title))
  if (length(ylim)==1 && is.na(ylim)) ylim = c(0, 6e4)
  x <- sim.data$SS[[stat]]$x.yBP[accepted[,1],,drop=FALSE]; y <- sim.data$SS[[stat]]$y[accepted[,1],,drop=FALSE]
  sm <- data.frame(id = rep(1:nrow(x), each = ncol(x)), x = c(t(x)), y = c(t(y)))
  x <- obs.data$SS[[stat]]$x.yBP; y <- obs.data$SS[[stat]]$y
  om <- data.frame(id = rep(1:nrow(x), each = ncol(x)), x = c(t(x)), y = c(t(y)))
  #
  g <- ggplot(sm) +
    scale_x_log10(limits = xlim, oob = oob, labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    scale_y_continuous(limits = ylim,  oob = oob) +
    general_theme +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = "none") +
    ggtitle(title)
  g <- g + xlab(axis.labs[1]) + ylab(axis.labs[2])
  if (!tile) {
    g <- g + geom_step(aes(x=x, y=y, col=id, group=id), size = line.size, alpha = .7)  +
      scg
  } else {
    g <- g + geom_step(aes(x=x, y=y, group=id, col=as.numeric(id)), size = line.size, alpha = .7) +
     scg # scale_color_distiller(palette="Purples", direction=1)
  }
  g <- g + geom_step(data = om, aes(x=x, y=y), size = 1.2, col = "black", alpha = .9)
  return(g)
}


gg.ld.errorbars <- function(obs.data,
                        sim.data,
                        accepted,
                        title,
                        scg,
                        xlim = c(NA, NA),
                        show.LD.obs = TRUE,
                        show.legend = FALSE,
                        axis.labs = c("kya", ""),
                        oob = oob,
                        plot.obs.TLD = FALSE )
{
  oob <- scales::squish # ifelse(oob, scales::squish, scales::censor)
  show.legend <- FALSE
  if (is_empty(obs.data, sim.data, "LD", control.obs = FALSE)) return(empty_plot(title))
  sm <- as.data.frame(sim.data$SS$LD[accepted[,1],,drop=FALSE])
  sm$labels <- accepted[,2]
  sm$col <- 1:nrow(sm)
  if (nrow(sm)==1) sm$col <- as.factor(sm$col)
  sm <- sm[rev(order(sm$mean, na.last = FALSE)),,drop=FALSE]
  sm$id <- factor(rownames(sm), levels=rownames(sm))
  om <- c("low"=obs.data$SS$LD[,"mean"]-1.96*obs.data$SS$LD[,"SE"], "up"=obs.data$SS$LD[,"mean"]+1.96*obs.data$SS$LD[,"SE"])
  #
  g <- ggplot(data = sm)
  if (show.LD.obs) {
    g <- g + geom_vline(xintercept = c(om["low.mean"], om["up.mean"])*1e-3, col = rgb(0, .5, 0, alpha=.7))
  }
  g <- g +
    geom_segment(aes(x = (mean-0.674*SE)*1e-3, xend = (mean+0.674*SE)*1e-3, y = id, yend = id, col = col), alpha=1, size=2.2) + #IC50
    geom_segment(aes(x = (mean-1.96*SE)*1e-3, xend = (mean+1.96*SE)*1e-3, y = id, yend = id, col = col), alpha=.5, size=1.5) + #IC95
    geom_point(aes(x = mean*1e-3, y = id), alpha=.85, size=4, shape=15, col="white") +
    geom_point(aes(x = mean*1e-3, y = id), alpha=.85, size=1.5, shape=15, col="black") +
    scale_x_continuous(limits = xlim*1e-3, oob = oob) +
    general_theme +
    theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.text.y = element_text(color = "white")) +
    ggtitle(title)
  g <- g + xlab(axis.labs[1]) + ylab(axis.labs[2])
  if (nrow(sm)==1) {
    g <- g + scale_color_manual(values = "#cc4779") + theme(legend.position = "none")
  } else {
    if (show.legend) {
      g <- g + scg + #scale_color_viridis(discrete = T, option = "A", begin = .1, end = .9, labels = stringr::str_wrap(sm$labels, 7)) +
        theme(legend.position = "right") +
        labs(color="")
    } else {
      g <- g + scg + #scale_color_viridis(discrete = T, option = "A", begin = .1, end = .9) +
        theme(legend.position = "none")
    }
  }
  if (plot.obs.TLD) g <- g + geom_vline(xintercept = 50, linetype = "dashed", alpha = .7)
  return(g)
}


abc.plot.core <- function(obs.data,
                          sim.data,
                          accepted,
                          output,
                          title,
                          line.size,
                          point.size,
                          scg,
                          sfg,
                          n,
                          tile,
                          show.LD.obs,
                          show.legend,
                          psmc..ylim,
                          print.params = c("ncol" = 3, "width" = 18, "height" = 14, "scale" = 0.9),
                          plot.title = TRUE,
                          log10.crf = c(T, T),
                          xylims.crf = c("xmin"=5e3, "xmax"=1e6, "ymin"=0.1e-2, "ymax"=10e-2) ,
                          log10.sprime = c(T, F),
                          xylims.sprime = c("xmin"=1e4, "xmax"=1e6, "ymin"=0., "ymax"=1.),
                          show.RMSE = TRUE,
                          axis.labs = FALSE,
                          ylim.afs.ceu = c(0, .32),
                          ylim.afs.yri = c(0, .4),
                          ylim.dcfs = c(0, .4),
                          xlim.D = c(-0.02, .15),
                          xlim.LD = c(0, 150),
                          transparent.bg = FALSE,
                          squeeze = TRUE,
													layout = "A",
													plot.obs.TLD = FALSE )

{
  CorF <- if (tile) sfg else scg
  oob <- squeeze

	if (layout == "A") {
	  G <- list(

	            gg.scatter(obs.data, sim.data, accepted,
	                       stat.x = "PI:CEU", se.x = NA, stat.y = "PI:YRI", se.y = NA,
	                       title = "(A) Nucleotide diversity", xlim = c(.5, 1.8), ylim = c(.5, 1.8),
	                       xlab = "", ylab = "",
	                       lab.percent = c(FALSE, FALSE),
	                       line.size=line.size, point.size=point.size, scg = CorF, tile = tile, show.legend = show.legend,
	                       axis.labs = dplyr::if_else(rep(axis.labs, 2), c("CEU", "YRI"), c("", "")),
	                       oob = oob ),

	            gg.sfs(obs.data, sim.data, accepted, stat = "AFS.CEU", title ="(B) AFS CEU", ylim = ylim.afs.ceu,
	                   line.size = line.size, scg = CorF, tile = tile, show.RMSE = show.RMSE,
	                   axis.labs = dplyr::if_else(rep(axis.labs, 2), c("", "Frequency"), c("", "")),
	                   oob = oob ),

	            gg.sfs(obs.data, sim.data, accepted, stat = "AFS.YRI", title ="(C) AFS YRI", ylim = ylim.afs.yri,
	                   line.size = line.size, scg = CorF, tile = tile, show.RMSE = show.RMSE,
	                   axis.labs = dplyr::if_else(rep(axis.labs, 2), c("", "Frequency"), c("", "")),
	                   oob = oob ),

	            gg.psmc(obs.data, sim.data, accepted,
	                    stat = "PSMC.CEU",
	                    title = "(D) PSMC CEU", xlim = c(1e4, 1e7), ylim = psmc..ylim,
	                    xlab = "", ylab = "", line.size=line.size, scg = scg, tile = tile,
	                    axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Age (ya, log10-scaled)", "IICR"), c("", "")),
	                    oob = oob ),

	            gg.scatter(obs.data, sim.data, accepted,
	                       stat.x = "D:D", se.x = "D:SE", stat.y = "Fst:1", se.y = NA,
	                       title = "(E) D & Fst", xlim = xlim.D, ylim = c(0, .23),
	                       lab.percent = c(TRUE, FALSE),
	                       xlab = "", ylab = "", line.size=line.size, point.size=point.size, scg = CorF, tile = tile,
	                       axis.labs = dplyr::if_else(rep(axis.labs, 2), c("D", "Fst"), c("", "")),
	                       oob = oob ),

	            gg.sfs(obs.data, sim.data, accepted, stat = "DCFS", title ="(F) DCFS Europeans", ylim = ylim.dcfs,
	                   line.size = line.size, scg = CorF, tile = tile, show.RMSE = show.RMSE,
	                   axis.labs = dplyr::if_else(rep(axis.labs, 2), c("", "Frequency"), c("", "")),
	                   oob = oob ),

	            gg.scatter(obs.data, sim.data, accepted,
	                       stat.x = "SPRIME:length.mean", stat.y = "SPRIME:match.mean",
	                       se.x = NA, se.y = NA, iqr.x = c("SPRIME:length.2.5p", "SPRIME:length.97.5p"), iqr.y = c("SPRIME:match.2.5p", "SPRIME:match.97.5p"),
	                       title = "(G) S\'", xlim = c(xylims.sprime["xmin"], xylims.sprime["xmax"]), ylim = c(xylims.sprime["ymin"], xylims.sprime["ymax"]),
	                       log10 = c(log10.sprime[1], log10.sprime[2]), lab.percent = c(FALSE, TRUE),
	                       xlab = "", ylab = "", line.size=line.size, point.size=point.size, scg = CorF, tile = tile,
	                       axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Tract length (bp)", "Match rate"), c("", "")),
	                       oob = oob ),

	            gg.scatter(obs.data, sim.data, accepted,
	                       stat.x = "CRF:length.mean", stat.y = "CRF:alpha.mean",
	                       se.x = NA, se.y = NA, iqr.x = c("CRF:length.2.5p", "CRF:length.97.5p"), iqr.y = c("CRF:alpha.2.5p", "CRF:alpha.97.5p"),
	                       title = "(H) CRF", xlim = c(xylims.crf["xmin"], xylims.crf["xmax"]), ylim = c(xylims.crf["ymin"], xylims.crf["ymax"]),
	                       log10 = c(log10.crf[1], log10.crf[2]),
	                       error.height.factor = 1e-3, error.width.factor = 3e-4, lab.percent = c(FALSE, TRUE),
	                       xlab = "", ylab = "", line.size=line.size, point.size=point.size, scg = CorF, tile = tile,
	                       axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Tract length (bp)", "Introgression rate"), c("", "")),
	                       oob = oob ),

	            #gg.ld(obs.data, sim.data, accepted, title = "(I) Admixture age", xlim = c(5e3, 1e6), n = n)
	            gg.ld.errorbars(obs.data, sim.data, accepted, title = "(I) Ancestry-LD decay rate",
	                            scg = scg, xlim = xlim.LD*1e3, show.LD.obs, show.legend,
	                            axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Age (kya)", ""), c("kya", "")),
	                            oob = oob, plot.obs.TLD = plot.obs.TLD )
	  )

	} else if (layout == "B") {
		G <- list(

	            gg.scatter(obs.data, sim.data, accepted,
	                       stat.x = "PI:CEU", se.x = NA, stat.y = "PI:YRI", se.y = NA,
	                       title = "(A) Nucleotide diversity", xlim = c(.5, 1.8), ylim = c(.5, 1.8),
	                       xlab = "", ylab = "",
	                       lab.percent = c(FALSE, FALSE),
	                       line.size=line.size, point.size=point.size, scg = CorF, tile = tile, show.legend = show.legend,
	                       axis.labs = dplyr::if_else(rep(axis.labs, 2), c("CEU", "YRI"), c("", "")),
	                       oob = oob ),

	            gg.sfs(obs.data, sim.data, accepted, stat = "AFS.CEU", title ="(B) AFS CEU", ylim = ylim.afs.ceu,
	                   line.size = line.size, scg = CorF, tile = tile, show.RMSE = show.RMSE,
	                   axis.labs = dplyr::if_else(rep(axis.labs, 2), c("", "Frequency"), c("", "")),
	                   oob = oob ),

	            gg.sfs(obs.data, sim.data, accepted, stat = "AFS.YRI", title ="(C) AFS YRI", ylim = ylim.afs.yri,
	                   line.size = line.size, scg = CorF, tile = tile, show.RMSE = show.RMSE,
	                   axis.labs = dplyr::if_else(rep(axis.labs, 2), c("", "Frequency"), c("", "")),
	                   oob = oob ),

	            gg.psmc(obs.data, sim.data, accepted,
	                    stat = "PSMC.CEU",
	                    title = "(D) PSMC CEU", xlim = c(1e4, 1e7), ylim = psmc..ylim,
	                    xlab = "", ylab = "", line.size=line.size, scg = scg, tile = tile,
	                    axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Age (ya, log10-scaled)", "IICR"), c("", "")),
	                    oob = oob ),

							gg.psmc(obs.data, sim.data, accepted,
	                    stat = "PSMC.YRI",
	                    title = "(E) PSMC YRI", xlim = c(1e4, 1e7), ylim = psmc..ylim,
	                    xlab = "", ylab = "", line.size=line.size, scg = scg, tile = tile,
	                    axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Age (ya, log10-scaled)", "IICR"), c("", "")),
	                    oob = oob ),

	            gg.scatter(obs.data, sim.data, accepted,
	                       stat.x = "D:D", se.x = "D:SE", stat.y = "Fst:1", se.y = NA,
	                       title = "(F) D & Fst", xlim = xlim.D, ylim = c(0, .23),
	                       lab.percent = c(TRUE, FALSE),
	                       xlab = "", ylab = "", line.size=line.size, point.size=point.size, scg = CorF, tile = tile,
	                       axis.labs = dplyr::if_else(rep(axis.labs, 2), c("D", "Fst"), c("", "")),
	                       oob = oob ),

	            gg.sfs(obs.data, sim.data, accepted, stat = "DCFS", title ="(G) DCFS Europeans", ylim = ylim.dcfs,
	                   line.size = line.size, scg = CorF, tile = tile, show.RMSE = show.RMSE,
	                   axis.labs = dplyr::if_else(rep(axis.labs, 2), c("", "Frequency"), c("", "")),
	                   oob = oob ),

	            gg.scatter(obs.data, sim.data, accepted,
	                       stat.x = "SPRIME:length.mean", stat.y = "SPRIME:match.mean",
	                       se.x = NA, se.y = NA, iqr.x = c("SPRIME:length.2.5p", "SPRIME:length.97.5p"), iqr.y = c("SPRIME:match.2.5p", "SPRIME:match.97.5p"),
	                       title = "(H) S\'", xlim = c(xylims.sprime["xmin"], xylims.sprime["xmax"]), ylim = c(xylims.sprime["ymin"], xylims.sprime["ymax"]),
	                       log10 = c(log10.sprime[1], log10.sprime[2]), lab.percent = c(FALSE, TRUE),
	                       xlab = "", ylab = "", line.size=line.size, point.size=point.size, scg = CorF, tile = tile,
	                       axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Tract length (bp)", "Match rate"), c("", "")),
	                       oob = oob ),

	            gg.ld.errorbars(obs.data, sim.data, accepted, title = "(I) Ancestry-LD decay rate",
	                            scg = scg, xlim = xlim.LD*1e3, show.LD.obs, show.legend,
	                            axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Age (kya)", ""), c("kya", "")),
	                            oob = oob, plot.obs.TLD = plot.obs.TLD )
	  )
		} else if (layout == "C") {
		  G <- list(
		    
		    gg.scatter(obs.data, sim.data, accepted,
		               stat.x = "PI:CEU", se.x = NA, stat.y = "PI:YRI", se.y = NA,
		               title = "(A) Nucleotide diversity", xlim = c(.5, 1.8), ylim = c(.5, 1.8),
		               xlab = "", ylab = "",
		               lab.percent = c(FALSE, FALSE),
		               line.size=line.size, point.size=point.size, scg = CorF, tile = tile, show.legend = show.legend,
		               axis.labs = dplyr::if_else(rep(axis.labs, 2), c("CEU", "YRI"), c("", "")),
		               oob = oob ),
		    
		    gg.sfs(obs.data, sim.data, accepted, stat = "AFS.CEU", title ="(B) AFS CEU", ylim = ylim.afs.ceu,
		           line.size = line.size, scg = CorF, tile = tile, show.RMSE = show.RMSE,
		           axis.labs = dplyr::if_else(rep(axis.labs, 2), c("", "Frequency"), c("", "")),
		           oob = oob ),
		    
		    gg.sfs(obs.data, sim.data, accepted, stat = "AFS.YRI", title ="(C) AFS YRI", ylim = ylim.afs.yri,
		           line.size = line.size, scg = CorF, tile = tile, show.RMSE = show.RMSE,
		           axis.labs = dplyr::if_else(rep(axis.labs, 2), c("", "Frequency"), c("", "")),
		           oob = oob ),
		    
		    gg.psmc(obs.data, sim.data, accepted,
		            stat = "PSMC.CEU",
		            title = "(D) PSMC CEU", xlim = c(1e4, 1e7), ylim = psmc..ylim,
		            xlab = "", ylab = "", line.size=line.size, scg = scg, tile = tile,
		            axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Age (ya, log10-scaled)", "IICR"), c("", "")),
		            oob = oob ),
		    
		    gg.psmc(obs.data, sim.data, accepted,
		            stat = "PSMC.YRI",
		            title = "(E) PSMC YRI", xlim = c(1e4, 1e7), ylim = psmc..ylim,
		            xlab = "", ylab = "", line.size=line.size, scg = scg, tile = tile,
		            axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Age (ya, log10-scaled)", "IICR"), c("", "")),
		            oob = oob ),
		    
		    gg.scatter(obs.data, sim.data, accepted,
		               stat.x = "D:D", se.x = "D:SE", stat.y = "Fst:1", se.y = NA,
		               title = "(F) D & Fst", xlim = xlim.D, ylim = c(0, .23),
		               lab.percent = c(TRUE, FALSE),
		               xlab = "", ylab = "", line.size=line.size, point.size=point.size, scg = CorF, tile = tile,
		               axis.labs = dplyr::if_else(rep(axis.labs, 2), c("D", "Fst"), c("", "")),
		               oob = oob ),
		    
		    gg.sfs(obs.data, sim.data, accepted, stat = "DCFS", title ="(G) DCFS Europeans", ylim = ylim.dcfs,
		           line.size = line.size, scg = CorF, tile = tile, show.RMSE = show.RMSE,
		           axis.labs = dplyr::if_else(rep(axis.labs, 2), c("", "Frequency"), c("", "")),
		           oob = oob ),
		    
		    gg.scatter(obs.data, sim.data, accepted,
		               stat.x = "SPRIME:length.mean", stat.y = "SPRIME:match.mean",
		               se.x = NA, se.y = NA, iqr.x = c("SPRIME:length.2.5p", "SPRIME:length.97.5p"), iqr.y = c("SPRIME:match.2.5p", "SPRIME:match.97.5p"),
		               title = "(H) S\'", xlim = c(xylims.sprime["xmin"], xylims.sprime["xmax"]), ylim = c(xylims.sprime["ymin"], xylims.sprime["ymax"]),
		               log10 = c(log10.sprime[1], log10.sprime[2]), lab.percent = c(FALSE, TRUE),
		               xlab = "", ylab = "", line.size=line.size, point.size=point.size, scg = CorF, tile = tile,
		               axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Tract length (bp)", "Match rate"), c("", "")),
		               oob = oob ),
		    
		    gg.scatter(obs.data, sim.data, accepted,
		               stat.x = "CRF:length.mean", stat.y = "CRF:alpha.mean",
		               se.x = NA, se.y = NA, iqr.x = c("CRF:length.2.5p", "CRF:length.97.5p"), iqr.y = c("CRF:alpha.2.5p", "CRF:alpha.97.5p"),
		               title = "(I) CRF", xlim = c(xylims.crf["xmin"], xylims.crf["xmax"]), ylim = c(xylims.crf["ymin"], xylims.crf["ymax"]),
		               log10 = c(log10.crf[1], log10.crf[2]),
		               error.height.factor = 1e-3, error.width.factor = 3e-4, lab.percent = c(FALSE, TRUE),
		               xlab = "", ylab = "", line.size=line.size, point.size=point.size, scg = CorF, tile = tile,
		               axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Tract length (bp)", "Introgression rate"), c("", "")),
		               oob = oob ),
		    
		    gg.ld.errorbars(obs.data, sim.data, accepted, title = "(J) Ancestry-LD decay rate",
		                    scg = scg, xlim = xlim.LD*1e3, show.LD.obs, show.legend,
		                    axis.labs = dplyr::if_else(rep(axis.labs, 2), c("Age (kya)", ""), c("kya", "")),
		                    oob = oob, plot.obs.TLD = plot.obs.TLD )
		  )
	} else {
		stop("Layout not implemented")
	}
  
  if (F) G <- lapply(G, function(gg) return(gg + theme(axis.title = element_text(size = 9, family = "Poppins"),
                                                       panel.grid = element_blank())  )  )
  
  if (transparent.bg) {
    G <- lapply(G, function(gg) return(gg + theme(panel.background = element_rect(fill = "transparent"),
                                                  plot.background = element_rect(fill = "transparent", color=NA),
                                                  legend.background = element_rect(fill = "transparent"),
                                                  legend.box.background = element_rect(fill = "transparent"))))
  }

  if (plot.title) {
    title_gg = ggplot() + labs(title=title, subtitle=format(Sys.time(), "%d %b %Y - %H:%M:%S"))
    GA <- suppressWarnings(cowplot::plot_grid(title_gg, cowplot::plot_grid(plotlist=G, labels = "", ncol = print.params["ncol"]), ncol = 1, rel_heights = c(.07, 1.)))
  } else {
    GA <- suppressWarnings(cowplot::plot_grid(plotlist=G, labels = "", ncol = print.params["ncol"]))
  }
  if (!is.null(output)) {
    if (transparent.bg) {
      ggsave(GA, filename = output, width = print.params["width"], height = print.params["height"], scale = print.params["scale"],  bg = "transparent")
    } else {
      ggsave(GA, filename = output, width = print.params["width"], height = print.params["height"], scale = print.params["scale"])
    }
  }
  return(GA)
}


abc.plot <- function(obs.data,
                     sim.data,
                     accepted,
                     output,
                     title,
                     line.size,
                     point.size,
                     scg,
                     sfg = NA,
                     n = 200,
                     plot.each = FALSE,
                     tile = FALSE,
                     show.LD.obs = TRUE,
                     show.legend = FALSE,
                     psmc..ylim = c(0, 60e3),
                     print.params = c("ncol" = 3, "width" = 18, "height" = 14, "scale" = 0.9),
                     plot.title = TRUE,
                     log10.crf = c(T, T),
                     xylims.crf = c("xmin"=5e3, "xmax"=1e6, "ymin"=0.1e-2, "ymax"=10e-2),
                     log10.sprime = c(T, F),
                     xylims.sprime = c("xmin"=1e4, "xmax"=1e6, "ymin"=0., "ymax"=1.),
                     show.RMSE = TRUE,
                     axis.labs = FALSE,
                     ylim.afs.ceu = c(0, .32),
                     ylim.afs.yri = c(0, .4),
                     ylim.dcfs = c(0, .4),
                     xlim.D = c(-0.02, .15),
                     xlim.LD = c(0, 150),
                     transparent.bg = FALSE,
                     squeeze = TRUE,
										 # 12.12.22
										 layout = "A",
										 # 28.03.23
										 plot.obs.TLD = FALSE )
{
  if (class(sfg)=="logical" && is.na(sfg)) sfg <- scg
  if (plot.each) {
    G <- abc.plot.core(obs.data, 
                       sim.data, 
                       accepted, 
                       output, 
                       title, 
                       line.size,
                       point.size, 
                       scg, 
                       sfg, 
                       n, 
                       tile = tile, 
                       show.LD.obs, 
                       show.legend, 
                       psmc..ylim, 
                       xlim.D = xlim.D, 
                       squeeze = squeeze,
											 layout = layout, 
											 plot.obs.TLD = plot.obs.TLD)
    pdf(paste0(output,".pdf"), width = 16, height = 14)
    for (i in 1:nrow(accepted)) print(abc.plot.core(obs.data,
                                                    sim.data,
                                                    accepted[i,,drop=FALSE],
                                                    output = NULL,
                                                    paste0(accepted[i,2],"                 d = ",round(accepted[i,3],2)),
                                                    line.size,
                                                    point.size,
                                                    scg,
                                                    sfg,
                                                    n,
                                                    tile = FALSE,
                                                    show.LD.obs,
                                                    show.legend,
                                                    psmc..ylim,
                                                    print.params,
                                                    plot.title,
                                                    log10.crf,
                                                    xylims.crf,
                                                    log10.sprime,
                                                    xylims.sprime,
                                                    show.RMSE,
                                                    axis.labs,
                                                    ylim.afs.ceu = ylim.afs.ceu,
                                                    ylim.afs.yri = ylim.afs.yri,
                                                    ylim.dcfs = ylim.dcfs,
                                                    xlim.D = xlim.D,
                                                    xlim.LD = xlim.LD,
                                                    transparent.bg = transparent.bg,
                                                    squeeze = squeeze,
																										layout = layout,
																										plot.obs.TLD = plot.obs.TLD))
    dev.off()
  } else {
    G <- abc.plot.core(obs.data,
                       sim.data,
                       accepted,
                       output,
                       title,
                       line.size,
                       point.size,
                       scg,
                       sfg,
                       n,
                       tile,
                       show.LD.obs,
                       show.legend,
                       psmc..ylim,
                       print.params,
                       plot.title,
                       log10.crf,
                       xylims.crf,
                       log10.sprime,
                       xylims.sprime,
                       show.RMSE,
                       axis.labs,
                       ylim.afs.ceu = ylim.afs.ceu,
                       ylim.afs.yri = ylim.afs.yri,
                       ylim.dcfs = ylim.dcfs,
                       xlim.D = xlim.D,
                       xlim.LD = xlim.LD,
                       transparent.bg = transparent.bg,
                       squeeze = squeeze,
											 layout = layout,
											 plot.obs.TLD = plot.obs.TLD )
  }
}

#___
