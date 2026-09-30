#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

source(file.path("scripts", "utils", "bootstrap.R"))
configure_project_library()
source(file.path("scripts", "utils", "paths.R"))
source(file.path("scripts", "utils", "model_helpers.R"))

paths <- project_paths()
ensure_output_directories(paths)
figure_dir <- file.path(paths$results, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

theme_ei <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "grey94", colour = "grey60"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom", legend.title = element_blank(),
      plot.title = element_blank(), plot.margin = margin(8, 10, 8, 8)
    )
}

save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ragg::agg_tiff(
    file.path(figure_dir, paste0(stem, ".tiff")), width = width,
    height = height, units = "in", res = 600, compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(
    file.path(figure_dir, paste0(stem, ".png")), width = width,
    height = height, units = "in", res = 200
  )
  print(plot)
  grDevices::dev.off()
}

wrap_text <- function(x, width = 34) {
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

# Figure 1: participant flow for the two primary analytic paths.
figure1_font <- "Arial"
figure1_text_size <- 4.0
panel_title_size <- 10
panel_title_margin <- margin(0, 0, 6, 0)
flow <- utils::read.csv(file.path(paths$results, "participant_flow.csv"))
flow <- flow[flow$branch %in% c(
  "Primary PFAS cohort", "Primary PFAS-metals overlap cohort"
), ]
flow$panel <- factor(
  flow$branch,
  levels = c("Primary PFAS cohort", "Primary PFAS-metals overlap cohort"),
  labels = c("A  Primary PFAS analysis", "B  Same-participant comparison")
)
flow$y <- 8 - flow$step_order
flow$label <- paste0(wrap_text(flow$step, 35), "\n", "n = ",
                     format(flow$n_remaining, big.mark = ",", trim = TRUE))
arrows <- do.call(rbind, lapply(split(flow, flow$panel), function(x) {
  data.frame(panel = unique(x$panel), x = 0, xend = 0,
             y = head(x$y, -1) - 0.25, yend = tail(x$y, -1) + 0.25,
             excluded = paste0("Excluded n = ",
                               format(
                                 tail(x$n_excluded_since_previous, -1),
                                 big.mark = ",", trim = TRUE
                               )))
}))
figure1 <- ggplot(flow, aes(x = 0, y = y)) +
  geom_segment(
    data = arrows, aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, colour = "grey35",
    arrow = grid::arrow(length = grid::unit(0.09, "in"))
  ) +
  geom_text(
    data = arrows, aes(x = 0.18, y = (y + yend) / 2, label = excluded),
    inherit.aes = FALSE, hjust = 0, size = figure1_text_size,
    family = figure1_font, colour = "grey35"
  ) +
  geom_label(
    aes(label = label), size = figure1_text_size, family = figure1_font,
    linewidth = 0.35,
    label.padding = grid::unit(0.14, "lines"), lineheight = 0.92,
    fill = "white", colour = "black"
  ) +
  facet_wrap(~panel, nrow = 1) +
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(0.5, 7.5), clip = "off") +
  theme_void(base_size = 12, base_family = figure1_font) +
  theme(
    strip.text = element_text(
      face = "bold", hjust = 0, size = panel_title_size,
      margin = panel_title_margin
    ),
    strip.background = element_blank(),
    panel.spacing = grid::unit(0.5, "in"), plot.margin = margin(8, 10, 8, 8)
  )
save_figure(figure1, "figure1_participant_flow", 10.5, 8.0)

# Figure 2: symmetric model estimates.
cox <- utils::read.csv(file.path(paths$results, "primary_cox_models.csv"))
cox$panel <- with(cox, ifelse(
  cohort == "pfas_primary", "A  PFAS primary cohort",
  ifelse(analyte == "pfas", "B  PFAS overlap cohort",
         ifelse(analyte == "lead", "C  Lead overlap cohort", "D  Cadmium overlap cohort"))
))
cox$model <- factor(cox$model, levels = rev(c("M0", "M1", "M2", "M3")))
model_colours <- c(M0 = "#4D4D4D", M1 = "#0072B2", M2 = "#D55E00", M3 = "#009E73")
figure2 <- ggplot(cox, aes(x = hr, y = model, colour = model)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = 0.18, linewidth = 0.55) +
  geom_point(size = 2.3) +
  facet_wrap(~panel, ncol = 2) +
  scale_x_log10(breaks = c(0.8, 0.9, 1, 1.1, 1.25, 1.4)) +
  scale_colour_manual(values = model_colours) +
  labs(x = "Hazard ratio per exposure doubling (95% CI)", y = NULL) +
  theme_ei(10) + theme(legend.position = "none")
save_figure(figure2, "figure2_symmetric_model_forest", 8.2, 6.5)

# Figure 3: coefficient changes and exposure-albumin concentration structure.
contrasts <- utils::read.csv(file.path(paths$results, "coefficient_contrasts.csv"))
contrasts <- contrasts[contrasts$contrast == "albumin_conditional_egfr", ]
structure <- utils::read.csv(file.path(paths$results, "albumin_exposure_structure.csv"))
analyte_labels <- c(pfas = "PFAS", lead = "Lead", cadmium = "Cadmium")
contrasts$label <- factor(analyte_labels[contrasts$analyte],
                          levels = rev(c("PFAS", "Lead", "Cadmium")))
structure$label <- factor(analyte_labels[structure$analyte],
                          levels = rev(c("PFAS", "Lead", "Cadmium")))
analyte_colours <- c(PFAS = "#D55E00", Lead = "#0072B2", Cadmium = "#009E73")

panel_a <- ggplot(contrasts, aes(x = estimate, y = label, colour = label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = 0.17, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = analyte_colours) +
  labs(x = expression(beta(M3) - beta(M1)~": albumin | eGFR"), y = NULL,
       subtitle = "A  Mortality coefficient change") +
  theme_ei(10) + theme(legend.position = "none", plot.subtitle = element_text(face = "bold"))

structure$r2_label <- sprintf(
  "Deviance-based incremental\nR² = %.4f", structure$deviance_incremental_r2
)
panel_b <- ggplot(structure, aes(x = beta_albumin, y = label, colour = label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = 0.17, linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_text(aes(x = ci_high + 0.025, label = r2_label), hjust = 0,
            size = 2.7, lineheight = 0.9, colour = "black") +
  scale_colour_manual(values = analyte_colours) +
  coord_cartesian(xlim = c(-0.18, 0.72), clip = "off") +
  labs(x = "Adjusted change in log₂ biomarker per 1-g/dL albumin", y = NULL,
       subtitle = "B  Exposure-albumin structure") +
  theme_ei(10) + theme(
    legend.position = "none", plot.subtitle = element_text(face = "bold"),
    plot.margin = margin(8, 40, 8, 8)
  )
figure3 <- patchwork::wrap_plots(panel_a, panel_b, widths = c(1, 1.15))
save_figure(figure3, "figure3_albumin_conditioned_structure", 10.0, 4.4)

# Figure 4: competing structures used to motivate, not identify, conditioning.
nodes <- rbind(
  data.frame(panel = "A  Carrier- and handling-related concentration structure",
             id = c("exposure", "albumin", "kidney", "measured", "health", "mortality"),
             x = c(0, 0, 0, 1.3, 2.4, 2.4), y = c(2.2, 1.3, 0.4, 1.3, 0.4, 2.2),
             half_width = c(0.33, 0.30, 0.30, 0.30, 0.26, 0.20),
             half_height = c(0.11, 0.07, 0.07, 0.11, 0.07, 0.07),
             label = c("Underlying PFAS\nexposure", "Serum albumin", "Kidney function",
                       "Measured total\nserum PFAS", "Health status", "Mortality")),
  data.frame(panel = "B  Physiological confounding or reverse causation",
             id = c("exposure", "albumin", "kidney", "measured", "health", "mortality"),
             x = c(0, 1.2, 1.2, 2.4, 0, 2.4), y = c(2.2, 1.55, 0.85, 2.2, 0.35, 0.35),
             half_width = c(0.33, 0.30, 0.30, 0.30, 0.26, 0.20),
             half_height = c(0.11, 0.07, 0.07, 0.11, 0.07, 0.07),
             label = c("Underlying PFAS\nexposure", "Serum albumin", "Kidney function",
                       "Measured total\nserum PFAS", "Health status", "Mortality"))
)
edge_specs <- list(
  A = data.frame(
    from = c("exposure", "exposure", "albumin", "kidney", "health"),
    to = c("measured", "mortality", "measured", "measured", "mortality")
  ),
  B = data.frame(
    from = c("exposure", "health", "health", "health", "albumin", "kidney"),
    to = c("measured", "albumin", "kidney", "mortality", "measured", "measured")
  )
)
edges <- do.call(rbind, lapply(seq_along(edge_specs), function(i) {
  panel <- unique(nodes$panel)[[i]]
  spec <- edge_specs[[i]]
  from <- nodes[nodes$panel == panel, c("id", "x", "y")]
  to <- nodes[nodes$panel == panel, c(
    "id", "x", "y", "half_width", "half_height"
  )]
  a <- merge(spec, from, by.x = "from", by.y = "id")
  a <- merge(a, to, by.x = "to", by.y = "id", suffixes = c("", "end"))
  dx <- a$xend - a$x
  dy <- a$yend - a$y
  boundary_distance <- pmin(
    a$half_width / pmax(abs(dx), .Machine$double.eps),
    a$half_height / pmax(abs(dy), .Machine$double.eps)
  )
  a$xend <- a$xend - boundary_distance * dx
  a$yend <- a$yend - boundary_distance * dy
  a$panel <- panel
  a
}))
figure4_dag <- ggplot() +
  geom_segment(
    data = edges, aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey35", linewidth = 0.55,
    arrow = grid::arrow(length = grid::unit(0.08, "in"), type = "closed")
  ) +
  geom_label(data = nodes, aes(x = x, y = y, label = label), size = 3,
             family = figure1_font, linewidth = 0.35, fill = "white", lineheight = 0.9) +
  facet_wrap(~panel, nrow = 1) +
  coord_cartesian(xlim = c(-0.35, 2.75), ylim = c(-0.2, 2.5), clip = "off") +
  theme_void(base_size = 12, base_family = figure1_font) + theme(
    strip.text = element_text(
      face = "bold", hjust = 0, size = panel_title_size,
      margin = panel_title_margin
    ),
    strip.background = element_blank(),
    panel.spacing = grid::unit(0.5, "in"), plot.margin = margin(8, 10, 8, 8)
  )
save_figure(figure4_dag, "figure4_competing_dags", 10.0, 4.5)

# Figure S1: survey-weighted cycle-specific PFAS distribution.
data <- readRDS(file.path(paths$processed, "analysis_cohorts.rds"))
pfas_design <- make_domain_design(
  data, "domain_pfas_source", "cohort_pfas_primary", "weight_pfas_primary"
)
cycle_quantiles <- do.call(rbind, lapply(
  levels(droplevels(pfas_design$variables$cycle_factor)), function(cycle) {
    keep <- as.character(pfas_design$variables$cycle_factor) == cycle
    design_cycle <- pfas_design[keep, ]
    q <- survey::svyquantile(
      ~log2_sum4_pfas, design_cycle,
      quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9), ci = FALSE, na.rm = TRUE
    )
    values <- as.numeric(q[[1]])
    data.frame(cycle = cycle, q10 = values[[1]], q25 = values[[2]],
               median = values[[3]], q75 = values[[4]], q90 = values[[5]])
  }
))
cycle_quantiles$cycle <- factor(cycle_quantiles$cycle, levels = cycle_quantiles$cycle)
figure_s1_cycle_distribution <- ggplot(cycle_quantiles, aes(x = cycle, y = median)) +
  geom_linerange(aes(ymin = q10, ymax = q90), linewidth = 0.6, colour = "grey45") +
  geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 3.0, colour = "#0072B2") +
  geom_point(size = 2.0, colour = "white", fill = "white", shape = 21) +
  labs(x = "NHANES cycle", y = "log₂ Σ4PFAS (nmol/L)") +
  theme_ei(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_figure(
  figure_s1_cycle_distribution, "figure_s1_cycle_pfas_distribution", 7.2, 4.2
)

# Figure S2: correctly parameterized restricted cubic splines.
rcs <- utils::read.csv(file.path(paths$results, "rcs_curves.csv"))
rcs <- rcs[rcs$analyte %in% c("pfas_primary", "lead", "cadmium"), ]
rcs$panel <- factor(
  rcs$analyte, levels = c("pfas_primary", "lead", "cadmium"),
  labels = c("PFAS primary cohort", "Lead overlap cohort", "Cadmium overlap cohort")
)
figure_s2_rcs <- ggplot(rcs, aes(x = exposure_value, y = hr)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "grey75", alpha = 0.5) +
  geom_line(colour = "#0072B2", linewidth = 0.75) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  geom_vline(
    data = unique(rcs[c("panel", "reference")]),
    aes(xintercept = reference), linetype = 2, colour = "grey45"
  ) +
  facet_wrap(~panel, scales = "free_x", nrow = 1) +
  labs(x = "log₂ exposure", y = "Hazard ratio (reference: weighted median)") +
  theme_ei(9) + theme(legend.position = "none")
save_figure(figure_s2_rcs, "figure_s2_rcs", 10.0, 3.6)

# Figure S3: interval-specific and leave-one-cycle-out estimates.
interval <- utils::read.csv(file.path(paths$results, "ph_interval_estimates.csv"))
interval$label <- factor(interval$interval, levels = rev(interval$interval))
loo <- utils::read.csv(file.path(paths$results, "leave_one_cycle_out.csv"))
loo <- loo[loo$model == "M3", ]
loo$label <- factor(loo$excluded_cycle, levels = rev(unique(loo$excluded_cycle)))
panel_interval <- ggplot(interval, aes(x = hr, y = label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = 0.14, linewidth = 0.6,
                colour = "#0072B2") +
  geom_point(size = 2.2, colour = "#0072B2") +
  labs(x = "Interval-specific HR", y = NULL, subtitle = "A  Follow-up interval") +
  theme_ei(9) + theme(plot.subtitle = element_text(face = "bold"))
panel_loo <- ggplot(loo, aes(x = hr, y = label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = 0.14, linewidth = 0.6,
                colour = "#0072B2") +
  geom_point(size = 2.2, colour = "#0072B2") +
  labs(x = "M3 HR after excluding cycle", y = NULL, subtitle = "B  Leave-one-cycle-out") +
  theme_ei(9) + theme(plot.subtitle = element_text(face = "bold"))
figure_s3_time_cycle <- patchwork::wrap_plots(panel_interval, panel_loo)
save_figure(
  figure_s3_time_cycle, "figure_s3_time_cycle_sensitivity", 9.0, 4.7
)

# Figure S4: focused PFAS M3 sensitivity results.
primary <- utils::read.csv(file.path(paths$results, "primary_cox_models.csv"))
primary <- primary[primary$cohort == "pfas_primary" & primary$model == "M3", ]
primary$analysis_label <- "Complete case; spline albumin/eGFR"
secondary <- utils::read.csv(file.path(paths$results, "secondary_sensitivity_models.csv"))
secondary <- secondary[secondary$model == "M3" & secondary$analysis %in% c(
  "linear_albumin_egfr", "exclude_baseline_diabetes_cvd_cancer",
  "exclude_known_pregnancy"
), ]
label_map <- c(
  linear_albumin_egfr = "Linear albumin/eGFR",
  exclude_baseline_diabetes_cvd_cancer = "No baseline diabetes/CVD/cancer",
  exclude_known_pregnancy = "Known pregnancy excluded"
)
secondary$analysis_label <- unname(label_map[secondary$analysis])
mi <- utils::read.csv(file.path(paths$results, "mi_pfas_models.csv"))
mi <- mi[mi$model == "M3", ]
mi_row <- data.frame(
  hr = mi$hr, ci_low = mi$hr_ci_low, ci_high = mi$hr_ci_high,
  analysis_label = "Multiple imputation (m = 20)"
)
sensitivity <- rbind(
  primary[, c("hr", "ci_low", "ci_high", "analysis_label")],
  secondary[, c("hr", "ci_low", "ci_high", "analysis_label")], mi_row
)
sensitivity$analysis_label <- factor(
  sensitivity$analysis_label, levels = rev(sensitivity$analysis_label)
)
figure_s4_pfas_sensitivity <- ggplot(sensitivity, aes(x = hr, y = analysis_label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = 0.14, linewidth = 0.6,
                colour = "#0072B2") +
  geom_point(size = 2.2, colour = "#0072B2") +
  labs(x = "PFAS M3 hazard ratio (95% CI)", y = NULL) + theme_ei(9)
save_figure(
  figure_s4_pfas_sensitivity, "figure_s4_pfas_sensitivity", 7.5, 4.3
)
# Figure S5: component-specific PFAS estimates in the overlap cohort.
component_models <- utils::read.csv(file.path(paths$results, "congener_models.csv"))
component_models <- component_models[component_models$model %in% c("M0", "M3"), ]
component_models$analyte_label <- factor(
  toupper(component_models$analyte),
  levels = rev(c("PFOA", "PFOS", "PFHXS", "PFNA"))
)
component_models$model <- factor(component_models$model, levels = c("M0", "M3"))

component_contrasts <- utils::read.csv(file.path(
  paths$results, "congener_contrasts.csv"
))
component_contrasts <- component_contrasts[
  component_contrasts$contrast == "albumin_conditional_egfr",
]
component_contrasts$analyte_label <- factor(
  toupper(component_contrasts$analyte),
  levels = rev(c("PFOA", "PFOS", "PFHXS", "PFNA"))
)

component_panel_a <- ggplot(
  component_models,
  aes(x = hr, y = analyte_label, colour = model, shape = model)
) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high), orientation = "y",
    position = position_dodge(width = 0.45), width = 0.14, linewidth = 0.55
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 2.4) +
  scale_colour_manual(values = model_colours[c("M0", "M3")]) +
  scale_shape_manual(values = c(M0 = 16, M3 = 17)) +
  labs(
    x = "Hazard ratio per component doubling (95% CI)", y = NULL,
    subtitle = "A  Separate component models"
  ) +
  theme_ei(10) +
  theme(plot.subtitle = element_text(face = "bold"))

component_panel_b <- ggplot(
  component_contrasts,
  aes(x = estimate, y = analyte_label)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y", width = 0.14, linewidth = 0.6,
    colour = analyte_colours[["PFAS"]]
  ) +
  geom_point(size = 2.5, colour = analyte_colours[["PFAS"]]) +
  labs(
    x = expression(beta(M3) - beta(M1)~": albumin | eGFR"), y = NULL,
    subtitle = "B  Albumin-conditional coefficient change"
  ) +
  theme_ei(10) +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(face = "bold")
  )

figure_s5_congener <- patchwork::wrap_plots(
  component_panel_a, component_panel_b, widths = c(1.12, 1)
)
save_figure(figure_s5_congener, "figure_s5_congener_analysis", 9.2, 4.5)

message("Main and supplementary figures complete: ", figure_dir)
