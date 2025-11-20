#' Run full NMR analysis (direct, FH, unit, stratified + plots/tables)
#'
#' @param nmr_data   Data frame created by process_NMR()
#' @param cluster.info  cluster.info object (from surveyPrev::clusterInfo)
#' @param admin.info1   admin.info for admin1 (surveyPrev::adminInfo)
#' @param admin.info2   admin.info for admin2 (surveyPrev::adminInfo)
#' @param geo           sf POINT layer of cluster locations (with HH1 / cluster IDs)
#' @param geo_admin1    sf polygons used for admin1 maps (e.g. admin2 polygons with NAME_2)
#' @param geo_admin2    sf polygons used for admin2 maps (e.g. admin3 polygons or admin2 with full names)
#' @param figs_dir      Output directory for figures/tables
#' @param country_label Character label used in titles (e.g. "Sierra Leone")
#' @param report_csv    Optional path to report CSV to compare NMR (can be NULL)
#' @param report_nmr_col Name of NMR column in the report file (default "NMR")
#' @param report_join_col Name of admin1 column in the report file (default "admin1.name")
#' @param nat_benchmark Benchmark NMR in **probability** scale (default 12/1000)
#' @param exceed_threshold Threshold for exceedance plots (default 12/1000)
#'
#' @return A list with all fitted objects and direct estimates
#' @export
run_nmr_analysis <- function(
    nmr_data,
    cluster.info,
    admin.info1,
    admin.info2,
    geo,
    geo_admin1,
    geo_admin2,
    figs_dir       = "figs-NMR",
    country_label  = "Country",
    report_csv     = NULL,
    report_nmr_col = "NMR",
    report_join_col = "admin1.name",
    nat_benchmark  = 12/1000,
    exceed_threshold = 12/1000
) {
  dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n=== NMR Direct Estimates ===\n")

  # ---- direct estimates ----
  nmr_res_adm0 <- surveyPrev::directEST(
    data         = nmr_data,
    cluster.info = cluster.info,
    admin        = 0
  )
  cat("National NMR (per 1,000 live births):",
      round(nmr_res_adm0$res.admin0$direct.est * 1000), "\n")

  nmr_res_adm1 <- surveyPrev::directEST(
    data         = nmr_data,
    cluster.info = cluster.info,
    admin        = 1
  )

  nmr_res_adm2 <- surveyPrev::directEST(
    data         = nmr_data,
    cluster.info = cluster.info,
    admin        = 2
  )

  # ---- optional comparison with report ----
  if (!is.null(report_csv) && file.exists(report_csv)) {
    report <- read.csv(report_csv)
    report <- dplyr::left_join(
      report,
      nmr_res_adm1$res.admin1,
      by = setNames("admin1.name", report_join_col)
    )
    report$direct.est1000 <- round(report$direct.est * 1000)
    report$diff <- report[[report_nmr_col]] - report$direct.est1000
    cat("\nComparison with report:\n")
    print(report[, c(report_join_col, report_nmr_col, "direct.est1000", "diff")])

    g <- ggplot2::ggplot(report,
                         ggplot2::aes(x = .data[[report_nmr_col]],
                                      y = .data$direct.est1000)) +
      ggplot2::geom_point(color = "red", size = 3) +
      ggplot2::geom_label(ggplot2::aes(label = .data[[report_join_col]]),
                          fill = NA) +
      ggplot2::geom_abline(intercept = 0, slope = 1) +
      ggplot2::ggtitle("NMR: Report vs surveyPrev Estimates") +
      ggplot2::xlab("DE in the report") +
      ggplot2::ylab("DE surveyPrev") +
      ggplot2::theme_bw()

    ggplot2::ggsave(
      filename = file.path(figs_dir, "compare_with_report.pdf"),
      plot     = g,
      width    = 7, height = 8
    )
  }

  cat("\n=== NMR Model Fitting ===\n")

  # ---- Fay-Herriot admin1 ----
  nmr_FH_adm1 <- surveyPrev::fhModel(
    data         = nmr_data,
    cluster.info = cluster.info,
    admin.info   = admin.info1,
    admin        = 1,
    model        = "bym2",
    aggregation  = TRUE
  )

  # ---- Fay-Herriot admin2 ----
  bad_admin2 <- subset(nmr_res_adm2$res.admin2,
                       direct.var < 1e-30)$admin2.name.full
  bad_clusters <- subset(cluster.info$data,
                         admin2.name.full %in% bad_admin2)$cluster
  cat("Removing", length(bad_clusters),
      "clusters with low variance from admin2 model\n")

  nmr_FH_adm2 <- surveyPrev::fhModel(
    data         = subset(nmr_data, !cluster %in% bad_clusters),
    cluster.info = cluster.info,
    admin.info   = admin.info2,
    admin        = 2,
    model        = "bym2",
    aggregation  = TRUE
  )

  # ---- cluster-level models ----
  nmr_unit_adm1 <- surveyPrev::clusterModel(
    data          = nmr_data,
    cluster.info  = cluster.info,
    admin.info    = admin.info1,
    stratification = FALSE,
    model         = "bym2",
    admin         = 1,
    aggregation   = TRUE,
    CI            = 0.95
  )

  nmr_unit_adm2 <- surveyPrev::clusterModel(
    data          = nmr_data,
    cluster.info  = cluster.info,
    admin.info    = admin.info2,
    stratification = FALSE,
    model         = "bym2",
    admin         = 2,
    aggregation   = TRUE,
    CI            = 0.95
  )

  nmr_unit_strat_adm1 <- surveyPrev::clusterModel(
    data          = nmr_data,
    cluster.info  = cluster.info,
    admin.info    = admin.info1,
    stratification = TRUE,
    model         = "bym2",
    admin         = 1,
    aggregation   = TRUE,
    CI            = 0.95
  )

  nmr_unit_strat_adm2 <- surveyPrev::clusterModel(
    data          = nmr_data,
    cluster.info  = cluster.info,
    admin.info    = admin.info2,
    stratification = TRUE,
    model         = "bym2",
    admin         = 2,
    aggregation   = TRUE,
    CI            = 0.95
  )

  cat("\n=== Model Fitting Complete ===\n")
  cat("\n=== Generating Visualizations ===\n")

  ## ---------------- Interval plots ----------------
  p1 <- surveyPrev::intervalPlot(
    admin   = 1,
    group   = FALSE,
    compare = TRUE,
    model   = list(
      "Direct Estimates"          = nmr_res_adm1,
      "Fay-Herriot Model"         = nmr_FH_adm1,
      "Cluster-Level Model"       = nmr_unit_adm1,
      "Stratified Cluster Model"  = nmr_unit_strat_adm1
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
      "Direct Estimates"          = nmr_res_adm2,
      "Fay-Herriot Model"         = nmr_FH_adm2,
      "Cluster-Level Model"       = nmr_unit_adm2,
      "Stratified Cluster Model"  = nmr_unit_strat_adm2
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
  out1 <- nmr_res_adm1$res.admin1[, c("admin1.name", "direct.est", "cv")]
  colnames(out1)[2] <- "mean"
  out1$model <- "Direct Estimates"

  out2 <- nmr_FH_adm1$res.admin1[, c("admin1.name", "mean", "cv")]
  out2$model <- "Fay-Herriot Model"

  out3 <- nmr_unit_adm1$res.admin1[, c("admin1.name", "mean", "cv")]
  out3$model <- "Cluster-Level Model"

  out4 <- subset(nmr_unit_strat_adm1$res.admin1, type == "full")[
    , c("admin1.name", "mean", "cv")]
  out4$model <- "Stratified Cluster-Level Model"

  out <- rbind(out1, out2, out3, out4)
  out$model <- factor(out$model, levels = unique(out$model))

  g1 <- surveyPrev::mapPlot(
    data        = out,
    geo         = geo_admin1,
    by.data     = "admin1.name",
    by.geo      = "NAME_2",
    is.long     = TRUE,
    variable    = "model",
    value       = "mean",
    legend.label = "Mean",
    direction   = -1,
    ncol        = 4
  )

  g2 <- surveyPrev::mapPlot(
    data        = out,
    geo         = geo_admin1,
    by.data     = "admin1.name",
    by.geo      = "NAME_2",
    is.long     = TRUE,
    variable    = "model",
    value       = "cv",
    legend.label = "CV",
    ncol        = 4
  ) +
    ggplot2::scale_fill_viridis_c("CV", option = "E", direction = -1)

  p3 <- g1 / g2
  ggplot2::ggsave(
    filename = file.path(figs_dir, "map-admin1.pdf"),
    plot     = p3,
    width    = 13, height = 7
  )

  # Admin2 maps
  out1 <- nmr_res_adm2$res.admin2[, c("admin2.name.full", "direct.est", "cv")]
  colnames(out1)[2] <- "mean"
  out1$model <- "Direct Estimates"

  out2 <- nmr_FH_adm2$res.admin2[, c("admin2.name.full", "mean", "cv")]
  out2$model <- "Fay-Herriot Model"

  out3 <- nmr_unit_adm2$res.admin2[, c("admin2.name.full", "mean", "cv")]
  out3$model <- "Cluster-Level Model"

  out4 <- subset(nmr_unit_strat_adm2$res.admin2, type == "full")[
    , c("admin2.name.full", "mean", "cv")]
  out4$model <- "Stratified Cluster-Level Model"

  outadm2 <- rbind(out1, out2, out3, out4)
  outadm2$model <- factor(outadm2$model, levels = unique(outadm2$model))

  g1 <- surveyPrev::mapPlot(
    data        = outadm2,
    geo         = geo_admin2,
    by.data     = "admin2.name.full",
    by.geo      = "admin2.name.full",
    is.long     = TRUE,
    variable    = "model",
    value       = "mean",
    legend.label = "Mean",
    direction   = -1,
    ncol        = 4
  )

  g2 <- surveyPrev::mapPlot(
    data        = outadm2,
    geo         = geo_admin2,
    by.data     = "admin2.name.full",
    by.geo      = "admin2.name.full",
    is.long     = TRUE,
    variable    = "model",
    value       = "cv",
    legend.label = "CV",
    ncol        = 4
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
    res1    = nmr_res_adm1$res.admin1,
    res2    = nmr_FH_adm1$res.admin1,
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Fay-Herriot Model"
  )
  s2 <- surveyPrev::scatterPlot(
    res1    = nmr_res_adm1$res.admin1,
    res2    = nmr_FH_adm1$res.admin1,
    value1  = "direct.se",
    value2  = "sd",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Standard Errors",
    label1  = "Direct Estimation",
    label2  = "Fay–Herriot Model"
  )
  s3 <- surveyPrev::scatterPlot(
    res1    = nmr_res_adm1$res.admin1,
    res2    = nmr_unit_adm1$res.admin1,
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Cluster-Level Model"
  )
  s4 <- surveyPrev::scatterPlot(
    res1    = nmr_res_adm1$res.admin1,
    res2    = nmr_unit_adm1$res.admin1,
    value1  = "direct.se",
    value2  = "sd",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Standard Error",
    label1  = "Direct Estimation",
    label2  = "Cluster-Level Model"
  )
  s5 <- surveyPrev::scatterPlot(
    res1    = nmr_res_adm1$res.admin1,
    res2    = subset(nmr_unit_strat_adm1$res.admin1, type == "full"),
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin1.name",
    by.res2 = "admin1.name",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Stratified Cluster-Level Model"
  )
  s6 <- surveyPrev::scatterPlot(
    res1    = nmr_res_adm1$res.admin1,
    res2    = subset(nmr_unit_strat_adm1$res.admin1, type == "full"),
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
    res1    = nmr_res_adm2$res.admin2,
    res2    = nmr_FH_adm2$res.admin2,
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Fay-Herriot Model"
  )
  s2 <- surveyPrev::scatterPlot(
    res1    = nmr_res_adm2$res.admin2,
    res2    = nmr_FH_adm2$res.admin2,
    value1  = "direct.se",
    value2  = "sd",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Standard Errors",
    label1  = "Direct Estimation",
    label2  = "Fay–Herriot Model"
  )
  s3 <- surveyPrev::scatterPlot(
    res1    = nmr_res_adm2$res.admin2,
    res2    = nmr_unit_adm2$res.admin2,
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Cluster-Level Model"
  )
  s4 <- surveyPrev::scatterPlot(
    res1    = nmr_res_adm2$res.admin2,
    res2    = nmr_unit_adm2$res.admin2,
    value1  = "direct.se",
    value2  = "sd",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Standard Error",
    label1  = "Direct Estimation",
    label2  = "Cluster-Level Model"
  )
  s5 <- surveyPrev::scatterPlot(
    res1    = nmr_res_adm2$res.admin2,
    res2    = subset(nmr_unit_strat_adm2$res.admin2, type == "full"),
    value1  = "direct.est",
    value2  = "mean",
    by.res1 = "admin2.name.full",
    by.res2 = "admin2.name.full",
    title   = "Estimates",
    label1  = "Direct Estimation",
    label2  = "Stratified Cluster-Level Model"
  )
  s6 <- surveyPrev::scatterPlot(
    res1    = nmr_res_adm2$res.admin2,
    res2    = subset(nmr_unit_strat_adm2$res.admin2, type == "full"),
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
    ggplot2::ggtitle(paste("Cluster Locations -", country_label))
  ggplot2::ggsave(
    filename = file.path(figs_dir, "cluster-map.pdf"),
    plot     = p7,
    width    = 5, height = 5
  )

  tmp <- nmr_data[, c("cluster", "strata")]
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
    ggplot2::ggtitle("Cluster Locations by Urban/Rural Strata")
  ggplot2::ggsave(
    filename = file.path(figs_dir, "cluster-map-byUR.pdf"),
    plot     = p7b,
    width    = 10, height = 5
  )

  ## ---------------- Pairwise comparison ----------------
  g1 <- surveyPrev::compareEstimates(
    nmr_FH_adm1, return.plot = TRUE,
    title = "Fay-Herriot Model"
  )
  g2 <- surveyPrev::compareEstimates(
    nmr_unit_adm1, return.plot = TRUE,
    title = "Cluster-Level Model"
  )
  g3 <- surveyPrev::compareEstimates(
    nmr_unit_strat_adm1, return.plot = TRUE,
    title = "Stratified Cluster-Level Model"
  )
  p8 <- g1 + g2 + g3
  ggplot2::ggsave(
    filename = file.path(figs_dir, "pairwise-admin1.pdf"),
    plot     = p8,
    width    = 30, height = 10
  )

  ## ---------------- Ridge plots ----------------
  p9a <- surveyPrev::ridgePlot(x = nmr_FH_adm1, direction = -1) +
    ggplot2::ggtitle("Fay-Herriot Model") +
    ggplot2::geom_vline(xintercept = nat_benchmark, color = "red",
                        linetype = "dashed", alpha = 0.6, linewidth = 1)
  p9b <- surveyPrev::ridgePlot(x = nmr_unit_adm1, direction = -1) +
    ggplot2::ggtitle("Cluster-Level Model") +
    ggplot2::geom_vline(xintercept = nat_benchmark, color = "red",
                        linetype = "dashed", alpha = 0.6, linewidth = 1)
  p9c <- surveyPrev::ridgePlot(x = nmr_unit_strat_adm1, direction = -1) +
    ggplot2::ggtitle("Stratified Cluster-Level Model") +
    ggplot2::geom_vline(xintercept = nat_benchmark, color = "red",
                        linetype = "dashed", alpha = 0.6, linewidth = 1)
  p9 <- p9a + p9b + p9c
  ggplot2::ggsave(
    filename = file.path(figs_dir, "ridge-admin1.pdf"),
    plot     = p9,
    width    = 15, height = 7
  )

  p10 <- surveyPrev::ridgePlot(x = nmr_FH_adm2, direction = -1) +
    ggplot2::geom_vline(xintercept = nat_benchmark, color = "red",
                        linetype = "dashed", alpha = 0.6, linewidth = 1)
  p11 <- surveyPrev::ridgePlot(x = nmr_unit_adm2, direction = -1) +
    ggplot2::geom_vline(xintercept = nat_benchmark, color = "red",
                        linetype = "dashed", alpha = 0.6, linewidth = 1)
  p12 <- surveyPrev::ridgePlot(x = nmr_unit_strat_adm2, direction = -1) +
    ggplot2::geom_vline(xintercept = nat_benchmark, color = "red",
                        linetype = "dashed", alpha = 0.6, linewidth = 1)

  ggplot2::ggsave(file.path(figs_dir, "ridge-FH-admin2.pdf"),     p10, width = 15, height = 10)
  ggplot2::ggsave(file.path(figs_dir, "ridge-cluster-admin2.pdf"), p11, width = 15, height = 10)
  ggplot2::ggsave(file.path(figs_dir, "ridge-cluster-strat-admin2.pdf"), p12,
                  width = 15, height = 10)

  ## ---------------- Rank plots ----------------
  p12a <- surveyPrev::rankPlot(x = nmr_FH_adm1, direction = 1) +
    ggplot2::ggtitle("Fay-Herriot Model")
  p12b <- surveyPrev::rankPlot(x = nmr_unit_adm1, direction = 1) +
    ggplot2::ggtitle("Cluster-Level Model")
  p12c <- surveyPrev::rankPlot(x = nmr_unit_strat_adm1, direction = 1) +
    ggplot2::ggtitle("Stratified Cluster-Level Model")
  p13 <- p12a + p12b + p12c
  ggplot2::ggsave(
    filename = file.path(figs_dir, "rank-admin1.pdf"),
    plot     = p13,
    width    = 15, height = 7
  )

  p14 <- surveyPrev::rankPlot(x = nmr_FH_adm2, direction = 1)
  p15 <- surveyPrev::rankPlot(x = nmr_unit_adm2, direction = 1)
  p16 <- surveyPrev::rankPlot(x = nmr_unit_strat_adm2, direction = 1)
  ggplot2::ggsave(file.path(figs_dir, "rank-FH-admin2.pdf"),        p14, width = 15, height = 10)
  ggplot2::ggsave(file.path(figs_dir, "rank-cluster-admin2.pdf"),   p15, width = 15, height = 10)
  ggplot2::ggsave(file.path(figs_dir, "rank-cluster-strat-admin2.pdf"), p16,
                  width = 15, height = 10)

  ## ---------------- Exceedance plots ----------------
  p17 <- surveyPrev::exceedPlot(
    nmr_FH_adm1, threshold = exceed_threshold,
    exceed = TRUE, direction = -1,
    geo = geo_admin1, by.geo = "NAME_2", ylim = c(0, 1)
  ) + ggplot2::ggtitle("Fay-Herriot Model: Probability of > threshold")
  p18 <- surveyPrev::exceedPlot(
    nmr_FH_adm2, threshold = exceed_threshold,
    exceed = TRUE, direction = -1,
    geo = geo_admin2, by.geo = "admin2.name.full", ylim = c(0, 1)
  ) + ggplot2::ggtitle("Fay-Herriot Model: Probability of > threshold")
  p19 <- surveyPrev::exceedPlot(
    nmr_unit_adm1, threshold = exceed_threshold,
    exceed = TRUE, direction = -1,
    geo = geo_admin1, by.geo = "NAME_2", ylim = c(0, 1)
  ) + ggplot2::ggtitle("Cluster-Level Model: Probability of > threshold")
  p20 <- surveyPrev::exceedPlot(
    nmr_unit_adm2, threshold = exceed_threshold,
    exceed = TRUE, direction = -1,
    geo = geo_admin2, by.geo = "admin2.name.full", ylim = c(0, 1)
  ) + ggplot2::ggtitle("Cluster-Level Model: Probability of > threshold")
  p21 <- surveyPrev::exceedPlot(
    nmr_unit_strat_adm1, threshold = exceed_threshold,
    exceed = TRUE, direction = -1,
    geo = geo_admin1, by.geo = "NAME_2", ylim = c(0, 1)
  ) + ggplot2::ggtitle("Stratified Model: Probability of > threshold")
  p22 <- surveyPrev::exceedPlot(
    nmr_unit_strat_adm2, threshold = exceed_threshold,
    exceed = TRUE, direction = -1,
    geo = geo_admin2, by.geo = "admin2.name.full", ylim = c(0, 1)
  ) + ggplot2::ggtitle("Stratified Model: Probability of > threshold")

  ggplot2::ggsave(file.path(figs_dir, "exceed-FH-admin1.pdf"),                 p17, width = 7, height = 6)
  ggplot2::ggsave(file.path(figs_dir, "exceed-FH-admin2.pdf"),                 p18, width = 7, height = 6)
  ggplot2::ggsave(file.path(figs_dir, "exceed-cluster-admin1.pdf"),            p19, width = 7, height = 6)
  ggplot2::ggsave(file.path(figs_dir, "exceed-cluster-admin2.pdf"),            p20, width = 7, height = 6)
  ggplot2::ggsave(file.path(figs_dir, "exceed-cluster-strat-admin1.pdf"),      p21, width = 7, height = 6)
  ggplot2::ggsave(file.path(figs_dir, "exceed-cluster-strat-admin2.pdf"),      p22, width = 7, height = 6)

  ## ---------------- Summary tables ----------------
  admin1_summary <- data.frame(
    Region          = nmr_res_adm1$res.admin1$admin1.name,
    Direct_Est      = round(nmr_res_adm1$res.admin1$direct.est * 1000, 2),
    Direct_SE       = round(nmr_res_adm1$res.admin1$direct.se * 1000, 2),
    FH_Est          = round(nmr_FH_adm1$res.admin1$mean * 1000, 2),
    FH_SE           = round(nmr_FH_adm1$res.admin1$sd * 1000, 2),
    Unit_Est        = round(nmr_unit_adm1$res.admin1$mean * 1000, 2),
    Unit_SE         = round(nmr_unit_adm1$res.admin1$sd * 1000, 2),
    Unit_Strat_Est  = round(subset(nmr_unit_strat_adm1$res.admin1, type == "full")$mean * 1000, 2),
    Unit_Strat_SE   = round(subset(nmr_unit_strat_adm1$res.admin1, type == "full")$sd * 1000, 2)
  )
  utils::write.csv(admin1_summary,
                   file = file.path(figs_dir, "admin1_summary.csv"),
                   row.names = FALSE)

  out3 <- nmr_unit_adm2$res.admin2[, c("admin2.name.full", "mean", "sd")]
  out4 <- nmr_unit_strat_adm2$res.admin2[, c("admin2.name.full", "mean", "sd")]

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
    nmr_data,
    cluster.info$data,
    by = "cluster"
  )

  summary_table <- merged_data %>%
    dplyr::group_by(admin2.name.full) %>%
    dplyr::summarize(
      total_births    = dplyr::n(),
      neonatal_deaths = sum(value == 1, na.rm = TRUE),
      crude_nmr       = round((neonatal_deaths / total_births) * 1000, 2),
      .groups         = "drop"
    )
  utils::write.csv(summary_table,
                   file = file.path(figs_dir, "summary_table.csv"),
                   row.names = FALSE)

  cat("\n=== NMR Analysis Complete ===\n")
  cat("All plots saved to:", figs_dir, "\n")

  invisible(list(
    res_adm0          = nmr_res_adm0,
    res_adm1          = nmr_res_adm1,
    res_adm2          = nmr_res_adm2,
    FH_adm1           = nmr_FH_adm1,
    FH_adm2           = nmr_FH_adm2,
    unit_adm1         = nmr_unit_adm1,
    unit_adm2         = nmr_unit_adm2,
    unit_strat_adm1   = nmr_unit_strat_adm1,
    unit_strat_adm2   = nmr_unit_strat_adm2
  ))
}
