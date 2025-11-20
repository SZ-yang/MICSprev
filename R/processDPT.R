#' Run full DTP3 analysis (direct, FH, unit, stratified + plots/tables)
#'
#' @param dtp3_data   Data frame created by process_DTP3()
#' @param cluster.info  cluster.info object (from surveyPrev::clusterInfo)
#' @param admin.info1   admin.info for admin1 (surveyPrev::adminInfo)
#' @param admin.info2   admin.info for admin2 (surveyPrev::adminInfo)
#' @param geo           sf POINT layer of cluster locations (with HH1 / cluster IDs)
#' @param geo_admin1    sf polygons used for admin1 maps (e.g. admin2 polygons with NAME_2)
#' @param geo_admin2    sf polygons used for admin2 maps (e.g. admin3 polygons or admin2 with full names)
#' @param figs_dir      Output directory for figures/tables
#' @param country_label Character label used in titles (e.g. "Sierra Leone")
#' @param low_cov_threshold Threshold for "low coverage" exceedance (default 0.80)
#'
#' @return A list with all fitted objects and direct estimates
#' @export
run_dtp3_analysis <- function(
    dtp3_data,
    cluster.info,
    admin.info1,
    admin.info2,
    geo,
    geo_admin1,
    geo_admin2,
    figs_dir        = "figs-DTP3",
    country_label   = "Country",
    low_cov_threshold = 0.80
) {
  dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n=== DTP3 Direct Estimates ===\n")

  # ---- direct estimates ----
  dtp3_res_adm0 <- surveyPrev::directEST(
    data         = dtp3_data,
    cluster.info = cluster.info,
    admin        = 0
  )
  cat("National DTP3 coverage:",
      round(dtp3_res_adm0$res.admin0$direct.est, 4), "\n")

  dtp3_res_adm1 <- surveyPrev::directEST(
    data         = dtp3_data,
    cluster.info = cluster.info,
    admin        = 1
  )

  dtp3_res_adm2 <- surveyPrev::directEST(
    data         = dtp3_data,
    cluster.info = cluster.info,
    admin        = 2
  )

  cat("\n=== DTP3 Model Fitting ===\n")

  # ---- Fay-Herriot admin1 ----
  dtp3_FH_adm1 <- surveyPrev::fhModel(
    data         = dtp3_data,
    cluster.info = cluster.info,
    admin.info   = admin.info1,
    admin        = 1,
    model        = "bym2",
    aggregation  = TRUE
  )

  # ---- Fay-Herriot admin2 ----
  bad_admin2 <- subset(dtp3_res_adm2$res.admin2,
                       direct.var < 1e-30)$admin2.name.full
  bad_clusters <- subset(cluster.info$data,
                         admin2.name.full %in% bad_admin2)$cluster
  cat("Removing", length(bad_clusters),
      "clusters with low variance from admin2 model\n")

  dtp3_FH_adm2 <- surveyPrev::fhModel(
    data         = subset(dtp3_data, !cluster %in% bad_clusters),
    cluster.info = cluster.info,
    admin.info   = admin.info2,
    admin        = 2,
    model        = "bym2",
    aggregation  = TRUE
  )

  # ---- cluster-level models ----
  dtp3_unit_adm1 <- surveyPrev::clusterModel(
    data           = dtp3_data,
    cluster.info   = cluster.info,
    admin.info     = admin.info1,
    stratification = FALSE,
    model          = "bym2",
    admin          = 1,
    aggregation    = TRUE,
    CI             = 0.95
  )

  dtp3_unit_adm2 <- surveyPrev::clusterModel(
    data           = dtp3_data,
    cluster.info   = cluster.info,
    admin.info     = admin.info2,
    stratification = FALSE,
    model          = "bym2",
    admin          = 2,
    aggregation    = TRUE,
    CI             = 0.95
  )

  dtp3_unit_strat_adm1 <- surveyPrev::clusterModel(
    data           = dtp3_data,
    cluster.info   = cluster.info,
    admin.info     = admin.info1,
    stratification = TRUE,
    model          = "bym2",
    admin          = 1,
    aggregation    = TRUE,
    CI             = 0.95
  )

  dtp3_unit_strat_adm2 <- surveyPrev::clusterModel(
    data           = dtp3_data,
    cluster.info   = cluster.info,
    admin.info     = admin.info2,
    stratification = TRUE,
    model          = "bym2",
    admin          = 2,
    aggregation    = TRUE,
    CI             = 0.95
  )

  cat("\n=== Model Fitting Complete ===\n")
  cat("\n=== Generating Visualizations ===\n")

  ## ---------------- Interval plots ----------------
  p1 <- surveyPrev::intervalPlot(
    admin   = 1,
    group   = FALSE,
    compare = TRUE,
    model   = list(
      "Direct Estimates"          = dtp3_res_adm1,
      "Fay-Herriot Model"         = dtp3_FH_adm1,
      "Cluster-Level Model"       = dtp3_unit_adm1,
      "Stratified Cluster Model"  = dtp3_unit_strat_adm1
    )
  ) +
    ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::scale_color_discrete("")

  ggplot2::ggsave(
    filename = file.path(figs_dir, "Interval-admin1.pdf"),
    plot     = p1,
    width    = 8, height = 5
  )

  p2 <- surveyPrev::intervalPlot(
    admin   = 2,
    group   = FALSE,
    compare = TRUE,
    model   = list(
      "Direct Estimates"          = dtp3_res_adm2,
      "Fay-Herriot Model"         = dtp3_FH_adm2,
      "Cluster-Level Model"       = dtp3_unit_adm2,
      "Stratified Cluster Model"  = dtp3_unit_strat_adm2
    )
  ) +
    ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::scale_color_discrete("") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::facet_wrap(~admin1.name, ncol = 5, scale = "free_x")

  ggplot2::ggsave(
    filename = file.path(figs_dir, "Interval-admin2.pdf"),
    plot     = p2,
    width    = 20, height = 12
  )

  ## ---------------- Map plots ----------------
  # Admin1 maps
  out1 <- dtp3_res_adm1$res.admin1[, c("admin1.name", "direct.est", "cv")]
  colnames(out1)[2] <- "mean"
  out1$model <- "Direct Estimates"

  out2 <- dtp3_FH_adm1$res.admin1[, c("admin1.name", "mean", "cv")]
  out2$model <- "Fay-Herriot Model"

  out3 <- dtp3_unit_adm1$res.admin1[, c("admin1.name", "mean", "cv")]
  out3$model <- "Cluster-Level Model"

  out4 <- subset(dtp3_unit_strat_adm1$res.admin1, type == "full")[
    , c("admin1.name", "mean", "cv")]
  out4$model <- "Stratified Cluster-Level Model"

  out <- rbind(out1, out2, out3, out4)
  out$model <- factor(out$model, levels = unique(out$model))

  g1 <- surveyPrev::mapPlot(
    data         = out,
    geo          = geo_admin1,
    by.data      = "admin1.name",
    by.geo       = "NAME_2",
    is.long      = TRUE,
    variable     = "model",
    value        = "mean",
    legend.label = "Coverage",
    direction    = -1,
    ncol         = 4
  )

  g2 <- surveyPrev::mapPlot(
    data         = out,
    geo          = geo_admin1,
    by.data      = "admin1.name",
    by.geo       = "NAME_2",
    is.long      = TRUE,
    variable     = "model",
    value        = "cv",
    legend.label = "CV",
    ncol         = 4
  ) +
    ggplot2::scale_fill_viridis_c("CV", option = "E", direction = -1)

  p3 <- g1 / g2
  ggplot2::ggsave(
    filename = file.path(figs_dir, "map-admin1.pdf"),
    plot     = p3,
    width    = 13, height = 7
  )

  # Admin2 maps
  out1 <- dtp3_res_adm2$res.admin2[, c("admin2.name.full", "direct.est", "cv")]
  colnames(out1)[2] <- "mean"
  out1$model <- "Direct Estimates"

  out2 <- dtp3_FH_adm2$res.admin2[, c("admin2.name.full", "mean", "cv")]
  out2$model <- "Fay-Herriot Model"

  out3 <- dtp3_unit_adm2$res.admin2[, c("admin2.name.full", "mean", "cv")]
  out3$model <- "Cluster-Level Model"

  out4 <- subset(dtp3_unit_strat_adm2$res.admin2, type == "full")[
    , c("admin2.name.full", "mean", "cv")]
  out4$model <- "Stratified Cluster-Level Model"

  outadm2 <- rbind(out1, out2, out3, out4)
  outadm2$model <- factor(outadm2$model, levels = unique(outadm2$model))

  g1 <- surveyPrev::mapPlot(
    data         = outadm2,
    geo          = geo_admin2,
    by.data      = "admin2.name.full",
    by.geo       = "admin2.name.full",
    is.long      = TRUE,
    variable     = "model",
    value        = "mean",
    legend.label = "Coverage",
    direction    = -1,
    ncol         = 4
  )

  g2 <- surveyPrev::mapPlot(
    data         = outadm2,
    geo          = geo_admin2,
    by.data      = "admin2.name.full",
    by.geo       = "admin2.name.full",
    is.long      = TRUE,
    variable     = "model",
    value        = "cv",
    legend.label = "CV",
    ncol         = 4
  ) +
    ggplot2::scale_fill_viridis_c("CV", option = "E", direction = -1)

  p4 <- g1 / g2
  ggplot2::ggsave(
    filename = file.path(figs_dir, "map-admin2.pdf"),
    plot     = p4,
    width    = 13, height = 7
  )

  ## ---------------- Scatter plots ----------------
  # Admin1
  s1 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm1$res.admin1,
    res2    = dtp3_FH_adm1$res.admin1,
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Fay-Herriot Model"
  )
  s2 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm1$res.admin1,
    res2    = dtp3_FH_adm1$res.admin1,
    value1  = "direct.se",
    value2  = "sd",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Standard Errors",
    label1  = "Direct Estimation",
    label2  = "Fay–Herriot Model"
  )
  s3 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm1$res.admin1,
    res2    = dtp3_unit_adm1$res.admin1,
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Cluster-Level Model"
  )
  s4 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm1$res.admin1,
    res2    = dtp3_unit_adm1$res.admin1,
    value1  = "direct.se",
    value2  = "sd",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Standard Error",
    label1  = "Direct Estimation",
    label2  = "Cluster-Level Model"
  )
  s5 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm1$res.admin1,
    res2    = subset(dtp3_unit_strat_adm1$res.admin1, type == "full"),
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Stratified Cluster-Level Model"
  )
  s6 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm1$res.admin1,
    res2    = subset(dtp3_unit_strat_adm1$res.admin1, type == "full"),
    value1  = "direct.se",
    value2  = "sd",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Standard Error",
    label1  = "Direct Estimation",
    label2  = "Stratified Cluster-Level Model"
  )
  p5 <- (s1 + s3 + s5) / (s2 + s4 + s6)
  ggplot2::ggsave(
    filename = file.path(figs_dir, "scatter-admin1.pdf"),
    plot     = p5,
    width    = 9, height = 6
  )

  # Admin2
  s1 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm2$res.admin2,
    res2    = dtp3_FH_adm2$res.admin2,
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Fay-Herriot Model"
  )
  s2 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm2$res.admin2,
    res2    = dtp3_FH_adm2$res.admin2,
    value1  = "direct.se",
    value2  = "sd",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Standard Errors",
    label1  = "Direct Estimation",
    label2  = "Fay–Herriot Model"
  )
  s3 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm2$res.admin2,
    res2    = dtp3_unit_adm2$res.admin2,
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Cluster-Level Model"
  )
  s4 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm2$res.admin2,
    res2    = dtp3_unit_adm2$res.admin2,
    value1  = "direct.se",
    value2  = "sd",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Standard Error",
    label1  = "Direct Estimation",
    label2  = "Cluster-Level Model"
  )
  s5 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm2$res.admin2,
    res2    = subset(dtp3_unit_strat_adm2$res.admin2, type == "full"),
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Stratified Cluster-Level Model"
  )
  s6 <- surveyPrev::scatterPlot(
    res1    = dtp3_res_adm2$res.admin2,
    res2    = subset(dtp3_unit_strat_adm2$res.admin2, type == "full"),
    value1  = "direct.se",
    value2  = "sd",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Standard Error",
    label1  = "Direct Estimation",
    label2  = "Stratified Cluster-Level Model"
  )
  p6 <- (s1 + s3 + s5) / (s2 + s4 + s6)
  ggplot2::ggsave(
    filename = file.path(figs_dir, "scatter-admin2.pdf"),
    plot     = p6,
    width    = 9, height = 6
  )

  ## ---------------- Cluster location maps ----------------
  p7 <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = geo_admin2, fill = NA,
                     color = "gray70", linewidth = 1) +
    ggplot2::geom_sf(data = geo_admin1,
                     ggplot2::aes(fill = NAME_2),
                     color = "black", linewidth = 2, alpha = 0.2) +
    ggplot2::geom_sf(data = geo, color = "red",
                     size = 0.15, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle(paste("Cluster Locations - DTP3 Coverage -", country_label))
  ggplot2::ggsave(
    filename = file.path(figs_dir, "cluster-map.pdf"),
    plot     = p7,
    width    = 5, height = 5
  )

  tmp <- dtp3_data[, c("cluster", "strata")]
  tmp <- tmp[match(unique(tmp$cluster), tmp$cluster), ]
  colnames(tmp)[1] <- "HH1"
  geo2 <- merge(geo, tmp)

  p7b <- ggplot2::ggplot(geo2) +
    ggplot2::geom_sf(data = geo_admin2, fill = NA,
                     color = "gray70", linewidth = 1) +
    ggplot2::geom_sf(data = geo_admin1,
                     ggplot2::aes(fill = NAME_2),
                     color = "black", linewidth = 2, alpha = 0.2) +
    ggplot2::geom_sf(color = "red", size = 0.4, alpha = 0.5) +
    ggplot2::facet_wrap(~strata) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle("Cluster Locations by Urban/Rural Strata - DTP3")
  ggplot2::ggsave(
    filename = file.path(figs_dir, "cluster-map-byUR.pdf"),
    plot     = p7b,
    width    = 10, height = 5
  )

  ## ---------------- Pairwise comparison ----------------
  g1 <- surveyPrev::compareEstimates(
    dtp3_FH_adm1, return.plot = TRUE,
    title = "Fay-Herriot Model"
  )
  g2 <- surveyPrev::compareEstimates(
    dtp3_unit_adm1, return.plot = TRUE,
    title = "Cluster-Level Model"
  )
  g3 <- surveyPrev::compareEstimates(
    dtp3_unit_strat_adm1, return.plot = TRUE,
    title = "Stratified Cluster-Level Model"
  )
  p8 <- g1 + g2 + g3
  ggplot2::ggsave(
    filename = file.path(figs_dir, "pairwise-admin1.pdf"),
    plot     = p8,
    width    = 30, height = 10
  )

  ## ---------------- Ridge plots ----------------
  p9a <- surveyPrev::ridgePlot(x = dtp3_FH_adm1, direction = -1) +
    ggplot2::ggtitle("Fay-Herriot Model - DTP3 Coverage")
  p9b <- surveyPrev::ridgePlot(x = dtp3_unit_adm1, direction = -1) +
    ggplot2::ggtitle("Cluster-Level Model - DTP3 Coverage")
  p9c <- surveyPrev::ridgePlot(x = dtp3_unit_strat_adm1, direction = -1) +
    ggplot2::ggtitle("Stratified Cluster-Level Model - DTP3 Coverage")
  p9 <- p9a + p9b + p9c
  ggplot2::ggsave(
    filename = file.path(figs_dir, "ridge-admin1.pdf"),
    plot     = p9,
    width    = 15, height = 7
  )

  p10 <- surveyPrev::ridgePlot(x = dtp3_FH_adm2, direction = -1)
  p11 <- surveyPrev::ridgePlot(x = dtp3_unit_adm2, direction = -1)
  p12 <- surveyPrev::ridgePlot(x = dtp3_unit_strat_adm2, direction = -1)
  ggplot2::ggsave(file.path(figs_dir, "ridge-FH-admin2.pdf"),          p10, width = 15, height = 10)
  ggplot2::ggsave(file.path(figs_dir, "ridge-cluster-admin2.pdf"),     p11, width = 15, height = 10)
  ggplot2::ggsave(file.path(figs_dir, "ridge-cluster-strat-admin2.pdf"), p12,
                  width = 15, height = 10)

  ## ---------------- Rank plots ----------------
  p13a <- surveyPrev::rankPlot(x = dtp3_FH_adm1, direction = 1) +
    ggplot2::ggtitle("Fay-Herriot Model")
  p13b <- surveyPrev::rankPlot(x = dtp3_unit_adm1, direction = 1) +
    ggplot2::ggtitle("Cluster-Level Model")
  p13c <- surveyPrev::rankPlot(x = dtp3_unit_strat_adm1, direction = 1) +
    ggplot2::ggtitle("Stratified Cluster-Level Model")
  p13 <- p13a + p13b + p13c
  ggplot2::ggsave(
    filename = file.path(figs_dir, "rank-admin1.pdf"),
    plot     = p13,
    width    = 15, height = 7
  )

  p14 <- surveyPrev::rankPlot(x = dtp3_FH_adm2, direction = 1)
  p15 <- surveyPrev::rankPlot(x = dtp3_unit_adm2, direction = 1)
  p16 <- surveyPrev::rankPlot(x = dtp3_unit_strat_adm2, direction = 1)
  ggplot2::ggsave(file.path(figs_dir, "rank-FH-admin2.pdf"),        p14, width = 15, height = 10)
  ggplot2::ggsave(file.path(figs_dir, "rank-cluster-admin2.pdf"),   p15, width = 15, height = 10)
  ggplot2::ggsave(file.path(figs_dir, "rank-cluster-strat-admin2.pdf"), p16,
                  width = 15, height = 10)

  ## ---------------- Exceedance plots (low coverage) ----------------
  thr_label <- paste0(round(low_cov_threshold * 100), "%")

  p17 <- surveyPrev::exceedPlot(
    dtp3_FH_adm1, threshold = low_cov_threshold,
    exceed = FALSE, direction = 1,
    geo = geo_admin1, by.geo = "NAME_2", ylim = c(0, 1)
  ) + ggplot2::ggtitle(
    paste0("Fay-Herriot Model: Prob. Coverage < ", thr_label)
  )

  p18 <- surveyPrev::exceedPlot(
    dtp3_FH_adm2, threshold = low_cov_threshold,
    exceed = FALSE, direction = 1,
    geo = geo_admin2, by.geo = "admin2.name.full", ylim = c(0, 1)
  ) + ggplot2::ggtitle(
    paste0("Fay-Herriot Model: Prob. Coverage < ", thr_label)
  )

  p19 <- surveyPrev::exceedPlot(
    dtp3_unit_adm1, threshold = low_cov_threshold,
    exceed = FALSE, direction = 1,
    geo = geo_admin1, by.geo = "NAME_2", ylim = c(0, 1)
  ) + ggplot2::ggtitle(
    paste0("Cluster-Level Model: Prob. Coverage < ", thr_label)
  )

  p20 <- surveyPrev::exceedPlot(
    dtp3_unit_adm2, threshold = low_cov_threshold,
    exceed = FALSE, direction = 1,
    geo = geo_admin2, by.geo = "admin2.name.full", ylim = c(0, 1)
  ) + ggplot2::ggtitle(
    paste0("Cluster-Level Model: Prob. Coverage < ", thr_label)
  )

  p21 <- surveyPrev::exceedPlot(
    dtp3_unit_strat_adm1, threshold = low_cov_threshold,
    exceed = FALSE, direction = 1,
    geo = geo_admin1, by.geo = "NAME_2", ylim = c(0, 1)
  ) + ggplot2::ggtitle(
    paste0("Stratified Model: Prob. Coverage < ", thr_label)
  )

  p22 <- surveyPrev::exceedPlot(
    dtp3_unit_strat_adm2, threshold = low_cov_threshold,
    exceed = FALSE, direction = 1,
    geo = geo_admin2, by.geo = "admin2.name.full", ylim = c(0, 1)
  ) + ggplot2::ggtitle(
    paste0("Stratified Model: Prob. Coverage < ", thr_label)
  )

  ggplot2::ggsave(file.path(figs_dir, "exceed-FH-admin1.pdf"),                 p17, width = 7, height = 6)
  ggplot2::ggsave(file.path(figs_dir, "exceed-FH-admin2.pdf"),                 p18, width = 7, height = 6)
  ggplot2::ggsave(file.path(figs_dir, "exceed-cluster-admin1.pdf"),            p19, width = 7, height = 6)
  ggplot2::ggsave(file.path(figs_dir, "exceed-cluster-admin2.pdf"),            p20, width = 7, height = 6)
  ggplot2::ggsave(file.path(figs_dir, "exceed-cluster-strat-admin1.pdf"),      p21, width = 7, height = 6)
  ggplot2::ggsave(file.path(figs_dir, "exceed-cluster-strat-admin2.pdf"),      p22, width = 7, height = 6)

  ## ---------------- Summary tables ----------------
  admin1_summary <- data.frame(
    Region          = dtp3_res_adm1$res.admin1$admin1.name,
    Direct_Est      = round(dtp3_res_adm1$res.admin1$direct.est, 4),
    Direct_SE       = round(dtp3_res_adm1$res.admin1$direct.se, 4),
    FH_Est          = round(dtp3_FH_adm1$res.admin1$mean, 4),
    FH_SE           = round(dtp3_FH_adm1$res.admin1$sd, 4),
    Unit_Est        = round(dtp3_unit_adm1$res.admin1$mean, 4),
    Unit_SE         = round(dtp3_unit_adm1$res.admin1$sd, 4),
    Unit_Strat_Est  = round(subset(dtp3_unit_strat_adm1$res.admin1, type == "full")$mean, 4),
    Unit_Strat_SE   = round(subset(dtp3_unit_strat_adm1$res.admin1, type == "full")$sd, 4)
  )
  utils::write.csv(admin1_summary,
                   file = file.path(figs_dir, "admin1_summary.csv"),
                   row.names = FALSE)

  out3 <- dtp3_unit_adm2$res.admin2[, c("admin2.name.full", "mean", "sd")]
  out4 <- dtp3_unit_strat_adm2$res.admin2[, c("admin2.name.full", "mean", "sd")]

  sd_table <- data.frame(
    admin2.name.full = out3$admin2.name.full,
    unstratified_sd  = out3$sd,
    stratified_sd    = out4$sd,
    diff             = out3$sd - out4$sd
  )
  utils::write.csv(sd_table,
                   file = file.path(figs_dir, "sd_table.csv"),
                   row.names = FALSE)

  merged_data <- dplyr::left_join(
    dtp3_data,
    cluster.info$data,
    by = "cluster"
  )

  summary_table <- merged_data %>%
    dplyr::group_by(admin2.name.full) %>%
    dplyr::summarize(
      total_children = dplyr::n(),
      vaccinated     = sum(value == 1, na.rm = TRUE),
      crude_coverage = round(vaccinated / total_children, 4),
      .groups        = "drop"
    )
  utils::write.csv(summary_table,
                   file = file.path(figs_dir, "summary_table.csv"),
                   row.names = FALSE)

  cat("\n=== DTP3 Analysis Complete ===\n")
  cat("All plots saved to:", figs_dir, "\n")

  invisible(list(
    res_adm0          = dtp3_res_adm0,
    res_adm1          = dtp3_res_adm1,
    res_adm2          = dtp3_res_adm2,
    FH_adm1           = dtp3_FH_adm1,
    FH_adm2           = dtp3_FH_adm2,
    unit_adm1         = dtp3_unit_adm1,
    unit_adm2         = dtp3_unit_adm2,
    unit_strat_adm1   = dtp3_unit_strat_adm1,
    unit_strat_adm2   = dtp3_unit_strat_adm2
  ))
}
