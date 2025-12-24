#Load libraries ----
library(tidyverse)
library(viridis)
library(patchwork)
library(scales)
library(cowplot)
library(rlang)
library(grid)

#Load data ----
flux_data <- read.csv("data/output_data/flux_data.csv")
model_params <- read.csv("data/output_data/muthusamy_model_parameters.csv")
mass_cumulative <- read.csv("data/output_data/mass_cumulative.csv")
q50_mass_flux <- read.csv("data/output_data/quartile_mass_flux.csv")
metadata <- read.csv("data/input_data/metadata.csv")
runoff_data <- read.csv("data/output_data/runoff_data_vol.csv")
PSA_velocity <- read.csv("data/output_data/PSA_velocity.csv")
PSA_perc_mobil <- read.csv("data/output_data/PSA_percent_mobilized.csv")
vid_data <- read.csv("data/output_data/video_data.csv")

#Wash off model plots ----
#Facet order + labels
params_surface <- model_params %>%
  mutate(
    rainfall_f = factor(rainfall, levels = c("high","med","low")),
    slope_f    = factor(slope,    levels = c(5,10,15)),
    surface_f  = factor(surface, levels = c("sand","concrete"))
  )

facet_labels <- list(
  slope_f    = c(`5`="5º (9%)", `10`="10º (18%)", `15`="15º (27%)"),
  rainfall_f = c(high="89 mm/hr", med="63 mm/hr", low="36 mm/hr"),
  surface_f = c("High Roughness","Low Roughness")
)

#Make curves from model parameters
model_curves <- model_params %>%
  rowwise() %>%
  mutate(curve = list(tibble(
    t = seq(0, 30, length.out = 200),
    pred_frac = f_k * (1 - exp(-(k_prime * i_mean) * t))
  ))) %>%
  unnest(curve) %>%
  select(condition, t, pred_frac)

runoff_prepped <- runoff_data %>%
  mutate(interval = end_time_min - start_time_min,
         volume_mL = if_else(interval == 0.5, volume_mL * 2, volume_mL)) %>%
  group_by(condition, end_time_min) %>%
  summarise(Q = mean(volume_mL, na.rm = TRUE), .groups = "drop") %>%
  arrange(condition, end_time_min) %>%
  group_by(condition) %>%
  #Add (0,0) flow anchor
  group_modify(~ bind_rows(tibble(end_time_min = 0, Q = 0), .x)) %>%
  ungroup()

#Replicate points to overlay
mass_points <- mass_cumulative %>%
  select(condition, end_time_min, cumulative_mass_g)

#Plot functions
plot_surface <- function(surface_type = "sand") {
  
  #Conditions for this surface
  conds_this_surface <- params_surface %>%
    filter(surface == surface_type) %>%
    pull(condition)
  
  #Join labels/params
  curves_surface <- model_curves %>%
    filter(condition %in% conds_this_surface) %>%
    left_join(params_surface %>% select(condition, surface, rainfall_f, slope_f, f_k, rse), by = "condition")
  
  runoff_surface <- runoff_prepped %>%
    filter(condition %in% conds_this_surface) %>%
    left_join(params_surface %>% select(condition, surface, rainfall_f, slope_f), by = "condition")
  
  mass_surface <- mass_points %>%
    filter(condition %in% conds_this_surface) %>%
    left_join(params_surface %>% select(condition, surface, rainfall_f, slope_f), by = "condition")
  
  #Dynamic annotation placement per facet
  annos <- mass_surface %>%
    group_by(condition, rainfall_f, slope_f) %>%
    summarise(
      x_anno = 0.35 * max(end_time_min, na.rm = TRUE),
      y_anno = 0.5 * max(c(cumulative_mass_g,
                            curves_surface %>% filter(condition == first(condition)) %>% pull(pred_frac)), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(params_surface %>% select(condition, k_prime, f_k, rse), by = "condition")
  
  annos <- annos %>%
    mutate(y_anno = if_else(condition == 4 | condition == 5 | condition == 6, y_anno + 0.5, y_anno),
           x_anno = if_else(condition == 4 | condition == 5 | condition == 6, x_anno - 10, x_anno)) %>%
    mutate(y_anno = if_else(condition == 1 | condition == 2 | condition == 3, y_anno - 0.15, y_anno)) %>%
    mutate(y_anno = if_else(condition == 16 | condition == 17 | condition == 18, y_anno + 0.5, y_anno),
           x_anno = if_else(condition == 16 | condition == 17 | condition == 18, x_anno - 10, x_anno))
  
  if (surface_type == "sand") {
    title = "High Roughness Wash-off Models"
  } else {
    title = "Low Roughness Wash-off Models" } 
  
  ggplot() +
    
    #Flow (secondary axis)
    geom_line(data = runoff_surface,
              aes(x = end_time_min, y = Q/1000),
              color = "lightblue", linewidth = 0.7) +
    #Observed points
    geom_point(data = mass_surface,
               aes(x = end_time_min, y = cumulative_mass_g),
               alpha = 0.6) +
    #Model curve
    geom_line(data = curves_surface,
              aes(x = t, y = pred_frac),
              color = "black", linewidth = 1) +
    #Annotations
    geom_text(data = annos,
              aes(x = x_anno, y = y_anno,
                  label = paste0("f(k) = ", format(f_k, digits = 3, nsmall = 3),
                                 "\nk' = ", round(k_prime, 3),
                                 "\nRSE = ", format(rse, digits = 3, nsmall = 3))),
              hjust = 0,
              size = 3.3) +
    facet_grid(
      rows = vars(rainfall_f),
      cols = vars(slope_f),
      labeller = labeller(.rows = facet_labels$rainfall_f, .cols = facet_labels$slope_f)
    ) +
    scale_y_continuous(
      name = "Wash-off fraction (" ~ W[t] ~ "/" ~ W[0] ~ ")",
      sec.axis = sec_axis(~ . * 1000, name = "Q (mL/min)")
    ) +
    scale_x_continuous(name = "Time (min)", limits = c(0, 30)) +
    ggtitle(title) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
}

plot_sand      <- plot_surface("sand")
plot_concrete  <- plot_surface("concrete")

#View
plot_sand
plot_concrete

#Save
ggsave("figures/washoff_model_fit_sand.png", plot_sand, width = 7, height = 5, dpi = 600)
ggsave("figures/washoff_model_fit_concrete.png", plot_concrete, width = 7, height = 5, dpi = 600)
#ggsave("figures/washoff_model_fit_sand.pdf", plot_sand, width = 7, height = 5, dpi = 600)
#ggsave("figures/washoff_model_fit_concrete.pdf", plot_concrete, width = 7, height = 5, dpi = 600)

##Legent for wash-off model plots ----
#Dummy data
legend_df <- tibble(
  x    = 1, y    = 1,
  type = factor(
    c("measured mass flux",
      "wash-off model fit",
      "measured runoff volume"),
    levels = c("measured mass flux",
               "wash-off model fit",
               "measured runoff volume")
  )
)

legend_plot <- ggplot(legend_df, aes(x, y, colour = type)) +
  geom_line(aes(linetype = type)) +
  geom_point(aes(shape = type)) +
  scale_colour_manual(
    values = c(
      "measured mass flux"   = "black",
      "wash-off model fit"   = "black",
      "measured runoff volume" = "lightblue"
    ),
    breaks = c("measured mass flux",
               "wash-off model fit",
               "measured runoff volume"),
    labels = c("measured mass flux",
               "wash-off model fit",
               "measured runoff volume"),
    name = NULL   # no legend title
  ) +
  # Shapes: only the first is a dot, others are "no point"
  scale_shape_manual(
    values = c(
      "measured mass flux"   = 16,  # filled circle
      "wash-off model fit"   = NA,
      "measured runoff volume" = NA
    )
  ) +
  # Linetypes: only the two lines get a solid line
  scale_linetype_manual(
    values = c(
      "measured mass flux"   = "blank",  # no line
      "wash-off model fit"   = "solid",
      "measured runoff volume" = "solid"
    )
  ) +
  guides(
    colour = guide_legend(
      direction = "horizontal",
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        shape    = c(16, NA, NA),
        linetype = c("blank", "solid", "solid")
      )
    ),
    shape = "none",
    linetype = "none"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 10)
  )

#Extract legend
legend_only <- cowplot::get_legend(legend_plot)

#Turn legend into a ggplot object for saving
legend_gg <- cowplot::ggdraw(legend_only)

#Save
ggsave("figures/washoff_legend_only.png",
       legend_gg,
       width = 6, height = 1, dpi = 600, bg = "white")



#Heatmaps for f(k), k', total wash-off, and Q50 ----
#Prep factors and labels
params <- model_params %>%
  mutate(
    slope_f    = factor(slope, levels = c(5, 10, 15)),
    rainfall_f = factor(rainfall, levels = c("low", "med", "high")),
    surface_f  = factor(surface, levels = c("sand", "concrete"))
  )

q50_mass_flux <- q50_mass_flux %>%
  mutate(
    slope_f    = factor(slope, levels = c(5, 10, 15)),
    rainfall_f = factor(rainfall, levels = c("low", "med", "high")),
    surface_f  = factor(surface, levels = c("sand", "concrete"))
  )

vid_data <- vid_data %>%
  mutate(
    slope_f    = factor(slope, levels = c(5, 10, 15)),
    rainfall_f = factor(rainfall, levels = c("low", "med", "high")),
    surface_f  = factor(surface, levels = c("sand", "concrete"))
  )

slope_tick_labels <- c(`5` = "5º", `10` = "10º", `15` = "15º")
rain_tick_labels  <- c(low = "Low", med = "Medium", high = "High")

#Function to make one 3×3 heatmap figure (faceted by surface rows)
make_heatmap <- function(df, value_col, title, limits, palette = "inferno", direction = 1) {
  ggplot(df, aes(x = slope_f, y = rainfall_f, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.7) +
    facet_grid(
      rows = vars(surface_f),
      labeller = labeller(surface_f = c(sand = "High", concrete = "Low"))
    ) +
    scale_x_discrete(name = "Slope", labels = slope_tick_labels) +
    scale_y_discrete(name = "Rainfall", labels = rain_tick_labels) +
    scale_fill_viridis_c(
      option = palette,
      limits = limits,
      name = title,
      direction = direction
    ) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text.y = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.title = element_text(face = "bold")
    )
}

#Set limits
kp_lim <- range(params$k_prime, na.rm = TRUE)
fk_lim <- range(params$f_k, na.rm = TRUE)
q50_lim <- range(q50_mass_flux$Q50_time_min_mean, na.rm = TRUE)
max_frac_lim <- range(q50_mass_flux$max_frac, na.rm = TRUE)

#Make individual plots
p_kp <- make_heatmap(params, "k_prime","k′",   kp_lim, palette = "magma")
p_fk <- make_heatmap(params, "f_k","f(k)",   fk_lim, palette = "magma")
p_q50 <- make_heatmap(q50_mass_flux, "Q50_time_min_mean","Q50 (min)",   q50_lim, palette = "magma", direction = -1)
p_max_frac <- make_heatmap(q50_mass_flux, "max_frac","Fw",   max_frac_lim, palette = "magma")

#Arrange side-by-side
heatmaps <- (p_q50 + theme(legend.position = "right")) +
  (p_max_frac + theme(legend.position = "right")) +
  (p_kp + theme(legend.position = "right")) +
  (p_fk + theme(legend.position = "right")) +
  plot_layout(ncol = 2, widths = c(1,1), guides = "keep")

#View
heatmaps

#Save figs
ggsave("figures/heatmaps.png", heatmaps, width = 9, height = 7, dpi = 600)
#ggsave("figures/heatmaps.pdf",  heatmaps, width = 10, height = 8)

##Make dissipation time heatmap for SI ----
dis_lim <- range(vid_data$dis_time_min, na.rm = TRUE)
p_dis <- make_heatmap(vid_data, "dis_time_min", paste0("Dissipation time",
                                                    "\n(mins)"),   dis_lim, palette = "magma", direction = -1)

ggsave("figures/dis_heatmap.png", p_dis, width = 4, height = 5, dpi = 600, bg = "white")





#Heatmap for flow depth ----
flow_data <- flux_data %>% group_by(slope, rainfall) %>%
  summarize(v_cm_s = mean(V) * 100,
            depth_mm = mean(depth_mm),
            tau_avg = mean(tau)) %>%
  mutate(
    slope_f    = factor(slope, levels = c(5, 10, 15)),
    rainfall_f = factor(rainfall, levels = c("low", "med", "high"))
  )

#Function to make one 3×3 heatmap figure (faceted by surface rows)
single_heatmap <- function(df, value_col, title, limits, palette = "inferno", direction = 1) {
  ggplot(df, aes(x = slope_f, y = rainfall_f, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.7) +
    scale_x_discrete(name = "Slope", labels = c(`5`="5º", `10`="10º", `15`="15º")) +
    scale_y_discrete(name = "Rainfall", labels = c(low="Low", med="Medium", high="High")) +
    scale_fill_viridis_c(option = palette, limits = limits, name = NULL, direction = direction) + 
    coord_equal() +
    ggtitle(title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text.y = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 6))
    )
}


#Set limits
flow_lim <- range(flow_data$depth_mm, na.rm = TRUE)
vel_lim <- range(flow_data$v_cm_s, na.rm = TRUE)
tau_lim <- range(flow_data$tau_avg, na.rm = TRUE)

#Make individual plots
p_flow <- single_heatmap(flow_data, "depth_mm","Flow Depth (mm)", flow_lim, palette = "magma")
p_vel <- single_heatmap(flow_data, "v_cm_s","Flow Velocity (cm/s)", vel_lim, palette = "magma")
p_tau <- single_heatmap(flow_data, "tau_avg","Shear Stress (Pa)", tau_lim, palette = "magma")

#Arrange side-by-side
heatmaps_flow <- (p_flow + theme(legend.position = "right")) +
  (p_vel + theme(legend.position = "right")) +
  (p_tau + theme(legend.position = "right")) +
  plot_layout(ncol = 3, widths = c(1,1,1), guides = "keep")

#View
heatmaps_flow

#Save fig
ggsave("figures/heatmaps_flow.png", heatmaps_flow, width = 12, height = 3, dpi = 600)

#Arrange stacked
heatmaps_flow <- (p_flow + theme(legend.position = "right")) +
  (p_vel + theme(legend.position = "right")) +
  (p_tau + theme(legend.position = "right")) +
  plot_layout(ncol = 1, widths = c(1), guides = "keep")

#View
heatmaps_flow

#Save fig
ggsave("figures/heatmaps_flow_stacked.png", heatmaps_flow, width = 5, height = 7, dpi = 600)



#Parameters vs flow conditions----
#Prep data
flux_by_cond <- flux_data %>% 
  group_by(condition, rainfall, slope, surface) %>%
  summarize(depth = mean(depth_mm),
            velocity = mean(V) * 100,
            shear = mean(tau)) %>%
  left_join(params %>% select(f_k, k_prime, condition, max_frac), by = "condition") %>%
  left_join(q50_mass_flux %>% select(Q50_time_min_mean, condition), by = "condition") %>%
  mutate(surface = factor(surface, levels = c("sand","concrete")))

#Set colors + shapes
mcols <- viridisLite::magma(100)
surf_cols   <- c(concrete = mcols[20], sand = mcols[60])
surf_shapes <- c(concrete = 16, sand = 17)  #16 = filled circle, 17 = filled triangle

#Labels for each variable
label_map <- c(
  k_prime = "k'",
  f_k = "f(k)",
  max_frac = "Fw",
  Q50_time_min_mean = "Q50 (min)",
  depth = "Flow depth (mm)",
  velocity = "Flow velocity (cm/s)",
  shear = "Shear stress (Pa)"
)

surface_labels <- c(
  sand = "High",
  concrete = "Low"
)

#Plot function
make_xy_plot <- function(df, xvar, yvar) {
  ggplot(df, aes_string(x = xvar, y = yvar, color = "surface", shape = "surface")) +
    geom_point(size = 3) +
    scale_color_manual(
      values = surf_cols,             
      breaks = c("sand","concrete"),  
      labels = c(sand = "High", concrete = "Low"),
      name   = "Surface\nRoughness"
    ) +
    scale_shape_manual(
      values = surf_shapes,            
      breaks = c("sand","concrete"),
      labels = c(sand = "High", concrete = "Low"),
      name   = "Surface\nRoughness"
    ) +
    labs(
      x = label_map[[xvar]],
      y = label_map[[yvar]],
      title = paste0(label_map[[yvar]], " vs ", label_map[[xvar]])
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
}

#Build all 12 combinations
x_vars <- c("depth", "velocity", "shear")
y_vars <- c("k_prime", "f_k", "Q50_time_min_mean", "max_frac")

plots_12 <- cross2(x_vars, y_vars) |>
  imap(function(vars, i) {
    x <- vars[[1]]; y <- vars[[2]]
    p <- make_xy_plot(flux_by_cond, x, y)
    attr(p, "plot_name") <- paste0("p_", y, "_vs_", x)
    p
  })

#Print one-by-one
for (p in plots_12) print(p)

#Correlation tests
cor.test(flux_by_cond$depth, flux_by_cond$k_prime, method = "spearman") #rho = 0.3477812
cor.test(flux_by_cond$velocity, flux_by_cond$k_prime, method = "spearman") #rho = 0.1455108
cor.test(flux_by_cond$shear, flux_by_cond$k_prime, method = "spearman") #rho = -0.07327141
cor.test(flux_by_cond$depth, flux_by_cond$f_k, method = "spearman") #rho = 0.6697626
cor.test(flux_by_cond$velocity, flux_by_cond$f_k, method = "spearman") #rho = 0.2899897
cor.test(flux_by_cond$shear, flux_by_cond$f_k, method = "spearman") #rho = -0.1021672
cor.test(flux_by_cond$depth, flux_by_cond$Q50_time_min_mean, method = "spearman") #rho = -0.7358101
cor.test(flux_by_cond$velocity, flux_by_cond$Q50_time_min_mean, method = "spearman") #rho = -0.620227
cor.test(flux_by_cond$shear, flux_by_cond$Q50_time_min_mean, method = "spearman") #rho = -0.1331269
cor.test(flux_by_cond$depth, flux_by_cond$max_frac, method = "spearman") #rho = 0.5892673
cor.test(flux_by_cond$velocity, flux_by_cond$max_frac, method = "spearman") #rho = 0.4365325
cor.test(flux_by_cond$shear, flux_by_cond$max_frac, method = "spearman") #rho = 0.07327141

#Function to add Spearman rho and p to plots (corner = "br" or "tr")
corr_annot <- function(df, xvar, yvar, corner = c("br","tr"),
                       digits_rho = 2, digits_p = 2, size = 4.5) {
  corner <- match.arg(corner)
  d <- df[, c(xvar, yvar)]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  
  ct   <- suppressWarnings(stats::cor.test(d[[xvar]], d[[yvar]],
                                           method = "spearman", exact = FALSE))
  rho  <- unname(ct$estimate)
  pval <- ct$p.value
  plab <- if (pval < 1e-3) {
    formatC(pval, format = "e", digits = 1)    #Scientific notation for very small vals
  } else {
    formatC(pval, format = "f", digits = 3)    #3 decimal places otherwise
  }
  
  xr <- range(d[[xvar]], na.rm = TRUE);    yr <- range(d[[yvar]], na.rm = TRUE)
  x  <- xr[2] - 0.02 * diff(xr) 
  y  <- if (corner == "br") yr[1] + 0.02 * diff(yr) else yr[2] - 0.02 * diff(yr)
  
  ggplot2::geom_text(
    data = data.frame(x = x, y = y,
                      lab = sprintf("\u03c1 = %.*f\n p = %s", digits_rho, rho, plab)),
    ggplot2::aes(x = x, y = y, label = lab),
    inherit.aes = FALSE, hjust = 1, vjust = if (corner == "br") 0 else 1,
    fontface = "bold", size = size
  )
}

#Name plots
names(plots_12) <- sapply(plots_12, function(p) attr(p, "plot_name"))

#Make grid of significant plots (rho > abs(0.5))
#Extract a single legend to place as a right column
legend_g <- cowplot::get_legend(
  plots_12$p_f_k_vs_depth +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12),
      legend.key.size = unit(1, "lines"),
      legend.spacing.y = unit(0.7, "lines")
    )
)
legend_gg <- cowplot::ggdraw(legend_g)

plots_12$p_f_k_vs_depth_n <- plots_12$p_f_k_vs_depth + theme(legend.position = "none") + corr_annot(flux_by_cond, "depth", "f_k", "br")
plots_12$p_Q50_time_min_mean_vs_depth_n <- plots_12$p_Q50_time_min_mean_vs_depth + theme(legend.position = "none") + corr_annot(flux_by_cond, "depth", "Q50_time_min_mean", "tr")
plots_12$p_max_frac_vs_depth_n <- plots_12$p_max_frac_vs_depth + theme(legend.position = "none") + corr_annot(flux_by_cond, "depth", "max_frac", "br")
plots_12$p_Q50_time_min_mean_vs_velocity_n <- plots_12$p_Q50_time_min_mean_vs_velocity + theme(legend.position = "none") + corr_annot(flux_by_cond, "velocity", "Q50_time_min_mean", "tr")

#Stitch plots
panel_grid <- (plots_12$p_f_k_vs_depth_n | plots_12$p_max_frac_vs_depth_n) / (plots_12$p_Q50_time_min_mean_vs_depth_n | plots_12$p_Q50_time_min_mean_vs_velocity_n)

final <- (panel_grid | legend_gg) + patchwork::plot_layout(widths = c(1, 0.2))

ggsave("figures/correlation_grid.png", final, width = 9, height = 7, dpi = 600, bg = "white")


#PSA Velocity Quartiles----
size_labels <- c("<125 µm","125-250 µm","250-500 µm","500-1000 µm",">1000 µm")
size_cols   <- setNames(
  magma(length(size_labels), begin = 0.05, end = 0.85), #Avoid near-black & near-white
  size_labels
)

PSA_clean <- PSA_velocity %>%
  mutate(
    size_class = str_replace_all(size_class, "–", "-"),
    size_f     = factor(size_class, levels = size_labels),
    rainfall_f = factor(rainfall,
                        levels = c("low","med","high"),
                        labels = c("Low","Medium","High")),
    slope_f    = factor(slope, levels = c(5,10,15)),
    surface_f  = factor(surface, levels = c("sand","concrete"),
                        labels = c("High Roughness","Low Roughness"))
  ) %>%
  mutate( #Convert from m/s to cm/s
    Q25_count_v = Q25_count_v * 100,
    Q50_count_v = Q50_count_v * 100,
    Q75_count_v = Q75_count_v * 100,
    Q25_vol_v = Q25_vol_v * 100,
    Q50_vol_v = Q50_vol_v * 100,
    Q75_vol_v = Q75_vol_v * 100
  )

plot_psa <- function(df, kind = c("count","volume")) {
  kind <- match.arg(kind)
  
  if (kind == "count") {
    q25 <- rlang::sym("Q25_count_v"); q50 <- rlang::sym("Q50_count_v"); q75 <- rlang::sym("Q75_count_v")
    title_lab <- expression("Quartile velocities by size")
  } else {
    q25 <- rlang::sym("Q25_vol_v");   q50 <- rlang::sym("Q50_vol_v");   q75 <- rlang::sym("Q75_vol_v")
    title_lab <- expression("Quartile velocities by size, calculated by particle volume")
  }
  
  #5 size classes per slope bin → dodge by size (order from size_labels)
  pos <- position_dodge(width = 0.85)
  tick_w <- 0.28 #Tick marks
  line_w <- 0.7
  
  ggplot(df, aes(x = slope_f, y = !!q50, color = size_f, group = size_f)) +
    #Vertical line from Q25 to Q75
    geom_linerange(aes(ymin = !!q25, ymax = !!q75),
                   position = pos, linewidth = line_w) +
    #Small tick at Q25
    geom_errorbar(aes(ymin = !!q25, ymax = !!q25),
                  position = pos, width = tick_w, linewidth = line_w, show.legend = FALSE) +
    #Small tick at Q75
    geom_errorbar(aes(ymin = !!q75, ymax = !!q75),
                  position = pos, width = tick_w, linewidth = line_w, show.legend = FALSE) +
    #Q50 point
    geom_point(position = pos, size = 2.5) +
    facet_grid(rows = vars(surface_f), cols = vars(rainfall_f)) +
    scale_color_manual(values = size_cols, name = "Size") +
    labs(x = "Slope (°)", y = "Velocity (cm/s)", title = title_lab, subtitle = "Rainfall Intensity") +
    coord_cartesian(ylim = c(0, NA)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.subtitle= element_text(hjust = 0.5,
                                  margin = margin(b = 6))
    )
}

#Build figures
p_psa_count  <- plot_psa(PSA_clean, "count")
p_psa_volume <- plot_psa(PSA_clean, "volume")

#View
p_psa_count
p_psa_volume

#Save
ggsave("figures/psa_quartile_count.png", p_psa_count, width = 9, height = 6, dpi = 600, bg = "white")
#ggsave("figures/psa_quartile_count.pdf",  p_psa_count, width = 9, height = 6)

ggsave("figures/psa_quartile_vol.png", p_psa_volume, width = 9, height = 6, dpi = 600, bg = "white")
#ggsave("figures/psa_quartile_vol.pdf",  p_psa_volume, width = 9, height = 6)



#Percent mobilized by size ----
#Palette & ordering
size_labels <- c("<125 µm","125-250 µm","250-500 µm","500-1000 µm",">1000 µm")
size_cols   <- setNames(magma(length(size_labels), begin = 0.05, end = 0.85), size_labels)

PSA_perc_clean <- PSA_perc_mobil %>%
  left_join(metadata %>% select(-c(run, replicate)) %>% distinct(), by = "condition") %>%
  mutate(
    size_class = str_replace_all(size_class, "–", "-"),
    size_f     = factor(size_class, levels = size_labels),
    slope_f    = factor(slope,    levels = c(5, 10, 15)),
    surface_f  = factor(surface, levels = c("sand","concrete"),
                        labels = c("High Roughness","Low Roughness")),
    rainfall_f = factor(rainfall,
                        levels = c("low","med","high"),
                        labels = c("Low","Medium","High")),
    fw_count_mobilized = perc_count_mobilized/100,
    fw_vol_mobilized = perc_volume_mobilized/100
  )


plot_perc_mobil_grouped <- function(df, kind = c("count","volume")) {
  kind <- match.arg(kind)
  value_sym <- if (kind == "count") rlang::sym("fw_count_mobilized") else rlang::sym("fw_vol_mobilized")
  title_lab <- if (kind == "count") expression("Fraction wash-off ("~F[w]~") by size") else expression("Fraction wash-off ("~F[w]~") by size, calculated by particle volume")
  
  #5 size classes per slope bin → dodge by size
  bar_width <- 0.75
  pos <- position_dodge2(width = 0.90, padding = 0.25, preserve = "single")
  
  ggplot(df, aes(x = slope_f, y = !!value_sym,
                 fill = size_f,            #Color by size class
                 group = size_f)) +        #Dodge by size within each slope
    geom_col(width = bar_width, position = pos, color = "white", linewidth = 0.3) +
    facet_grid(rows = vars(surface_f), cols = vars(rainfall_f)) +
    scale_fill_manual(values = size_cols, name = "Size") +
    scale_color_manual(values = size_cols, guide = "none") +
    labs(x = "Slope (°)", y = expression(F[w]), title = title_lab, subtitle = "Rainfall Intensity") +
    coord_cartesian(ylim = c(0.20, 1.00)) +
    scale_y_continuous(breaks = seq(0.20, 1.00, 0.20),
                       expand = expansion(mult = c(0, 0), add = 0)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing.y = unit(1.4, "lines"),
      plot.subtitle= element_text(hjust = 0.5,
                                  margin = margin(b = 6))
    )
}

#Build figures
p_mobil_count  <- plot_perc_mobil_grouped(PSA_perc_clean, "count")
p_mobil_volume <- plot_perc_mobil_grouped(PSA_perc_clean, "volume")

#View
p_mobil_count
p_mobil_volume

#Save
ggsave("figures/psa_mobilized_count.png", p_mobil_count, width = 9, height = 6, dpi = 600, bg = "white")
#ggsave("figures/psa_quartile_count.pdf",  p_mobil_count, width = 9, height = 6)

ggsave("figures/psa_mobilized_vol.png", p_mobil_volume, width = 9, height = 6, dpi = 600, bg = "white")
#ggsave("figures/psa_quartile_vol.pdf",  p_mobil_volume, width = 9, height = 6)



#CDFs: particle load by size ----
#Read in data (warning, takes a couple minutes)
PSA_data <- read.csv("data/input_data/PSA_data.csv") 

PSA_mobilized <- PSA_data %>%
  filter(interval < 7, interval > 0) %>%
  select(FLength, Volume, run, end_time_min) %>%
  left_join(metadata %>% select(run, condition) %>% distinct(),
            by = "run") %>%
  select(-run)

#Assign size classes
size_breaks <- c(-Inf, 125, 250, 500, 1000, Inf)
size_labels <- c("<125 µm", "125–250 µm", "250–500 µm", "500–1000 µm", ">1000 µm")
mobilized_size_class <- PSA_mobilized %>%
  mutate(size_class = cut(FLength, breaks = size_breaks, labels = size_labels))

#Set factors
size_cols   <- setNames(magma(length(size_labels), begin = 0.10, end = 0.85), size_labels)
cond_meta_f <- metadata %>%
  distinct(condition, surface, rainfall, slope) %>%
  mutate(
    surface_f  = factor(surface,  levels = c("sand","concrete"),
                        labels  = c("High Roughness","Low Roughness")),
    rainfall_f = factor(rainfall, levels = c("high","med","low"),
                        labels  = c("High","Medium","Low")),
    slope_f    = factor(slope,    levels = c(5,10,15))
  )

#Time scale
t_grid <- seq(0, 30, by = 5)

#Prep from particle data
cdf_by_condition <- mobilized_size_class %>%
  transmute(
    condition,
    size_f = factor(size_class, levels = size_labels),
    t      = as.numeric(end_time_min),         
    w      = as.numeric(Volume)
  ) %>%
  filter(!is.na(size_f), is.finite(t)) %>%
  group_by(condition, size_f, t) %>%
  summarise(n = n(), vol = sum(w, na.rm = TRUE), .groups = "drop") %>%
  group_by(condition, size_f) %>%
  complete(t = t_grid, fill = list(n = 0, vol = 0)) %>%
  arrange(t, .by_group = TRUE) %>%
  mutate(
    cum_n  = cumsum(n),
    cum_v  = cumsum(vol),
    N_tot  = sum(n),
    V_tot  = sum(vol),
    frac_count  = ifelse(N_tot > 0, cum_n / N_tot, NA_real_),
    frac_volume = ifelse(V_tot > 0, cum_v / V_tot, NA_real_)
  ) %>%
  ungroup() %>%
  left_join(cond_meta_f, by = "condition")


#Figure function
plot_cdf_surface <- function(cdf_df, surface_name = c("High Roughness","Low Roughness"),
                             kind = c("count","volume"), line_size = 0.9, pt_size = 1.2) {
  surface_name <- match.arg(surface_name)
  kind         <- match.arg(kind)
  
  y_col  <- if (kind == "count") rlang::sym("frac_count") else rlang::sym("frac_volume")
  title_ <- paste0("CDF over time by size class, ", surface_name,
                   " (by ", kind, ")")
  
  df <- cdf_df %>% dplyr::filter(surface_f == surface_name)
  
  ggplot(df, aes(x = t, y = !!y_col, color = size_f, group = size_f)) +
    geom_line(linewidth = line_size) +
    geom_point(size = pt_size) +
    geom_point(data = ~ dplyr::filter(.x, t == 0),
               size = pt_size + 0.4, show.legend = FALSE) +
    facet_grid(rows = vars(rainfall_f), cols = vars(slope_f)) +
    scale_color_manual(values = size_cols, name = "Size (µm)") +
    scale_x_continuous(name = "Time (min)", breaks = seq(0, 30, 5), limits = c(0, 30),
                       expand = expansion(add = c(0.5, 0))) +
    scale_y_continuous(name = "Cumulative fraction (within size class)",
                       labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    labs(title = title_) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

#Make figures
p_cdf_sand_count     <- plot_cdf_surface(cdf_by_condition, "High Roughness", "count")
p_cdf_sand_volume    <- plot_cdf_surface(cdf_by_condition, "High Roughness", "volume")
p_cdf_concrete_count <- plot_cdf_surface(cdf_by_condition, "Low Roughness", "count")
p_cdf_concrete_vol   <- plot_cdf_surface(cdf_by_condition, "Low Roughness", "volume")

#View
p_cdf_sand_count
p_cdf_sand_volume
p_cdf_concrete_count
p_cdf_concrete_vol

#Save (white background)
ggsave("figures/cdf_sand_count.png", p_cdf_sand_count, width = 7, height = 5, dpi = 600, bg = "white")
ggsave("figures/cdf_sand_volume.png", p_cdf_sand_volume, width = 7, height = 5, dpi = 600, bg = "white")
ggsave("figures/cdf_concrete_count.png", p_cdf_concrete_count, width = 7, height = 5, dpi = 600, bg = "white")
ggsave("figures/cdf_concrete_volume.png",p_cdf_concrete_vol, width = 7, height = 5, dpi = 600, bg = "white")





#PSD (All) ----
PSD_data <- PSA_data %>% select(FLength, Volume)

#KDE Function
kde_all <- function(df, weight_col = NULL, n = 1024, bw_mult = 1.8, x_range = NULL) {
  df <- df %>% transmute(size = as.numeric(FLength),
                         w    = if (is.null(weight_col)) 1 else pmax(0, as.numeric(.data[[weight_col]])))
  df <- df %>% filter(is.finite(size), is.finite(w), w > 0)
  
  if (nrow(df) < 10) return(tibble(size_um = numeric(), density = numeric()))
  
  if (is.null(x_range)) x_range <- range(df$size, na.rm = TRUE)
  
  bw  <- stats::bw.nrd0(df$size) * bw_mult
  den <- stats::density(df$size,
                        weights = df$w / sum(df$w),
                        bw      = bw,
                        n       = n,
                        from    = x_range[1],
                        to      = x_range[2],
                        cut     = 0)
  
  tibble(size_um = den$x, density = den$y)
}

#Plot builder function
psd_all_plot <- function(df, kind = c("count","volume"),
                         bw_mult = 1.8, fill_alpha = 0.22, line_size = 1.1) {
  kind <- match.arg(kind)
  
  #Choose weights & colors
  weight_col <- if (kind == "count") NULL else "Volume"
  title_lab  <- if (kind == "count") "PSD (all particles), by count" else "PSD (all particles), by volume"
  pal <- magma(11)
  col <- if (kind == "count") pal[7] else pal[9]
  
  dens <- kde_all(df, weight_col = weight_col, bw_mult = bw_mult, x_range = c(0, 2000))
  
  ggplot(dens, aes(size_um, density)) +
    geom_area(fill = alpha(col, fill_alpha), color = NA) +
    geom_line(color = col, linewidth = line_size) +
    labs(title = title_lab, x = "Particle size (µm)", y = "Relative frequency") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
}

#Build plots
p_psd_count  <- psd_all_plot(PSD_data, "count",  bw_mult = 3)
p_psd_volume <- psd_all_plot(PSD_data, "volume", bw_mult = 3)

#View
p_psd_count
p_psd_volume

#Side-by-side with patchwork
#p_psd_count | p_psd_volume

#Save with white background
ggsave("figures/psd_all_count.png",  p_psd_count,  width = 8, height = 5, dpi = 300, bg = "white")
ggsave("figures/psd_all_volume.png", p_psd_volume, width = 8, height = 5, dpi = 300, bg = "white")




#Shape data correlations ----
#Read in data (warning, takes a couple minutes)
shape_data <- read.csv("data/output_data/shape_data.csv")

#Choose colors
pal_full <- magma(11)
cols_7_9 <- pal_full[c(7, 9)]

#Set factors
surf_lvls <- levels(factor(shape_data$surface))
if (length(surf_lvls) == 2) {
  surf_cols <- setNames(cols_7_9, surf_lvls)
} else {
  surf_cols <- setNames(magma(length(surf_lvls)), surf_lvls)
}

#pal_full  <- viridisLite::magma(11)
#surf_cols <- c(concrete = pal_full[7], sand = pal_full[9])
#Reorder so sand comes first
order_levels <- c("sand","concrete")
shape_data$surface <- factor(shape_data$surface, levels = order_levels)

plot_shape_overlay <- function(df, x,
                               x_lab   = deparse(substitute(x)),
                               title   = NULL,
                               y_range = NULL,
                               rho_mode   = c("auto","manual","none"),
                               rho_vals   = NULL,
                               rho_method = c("spearman","pearson"),
                               rho_digits = 2,
                               beta_mode  = c("auto","manual","none"),
                               beta_vals  = NULL,
                               beta_digits = 2,
                               rho_size   = 5,
                               rho_gap    = 0.18) {
  
  rho_mode   <- match.arg(rho_mode)
  rho_method <- match.arg(rho_method)
  beta_mode  <- match.arg(beta_mode)
  
  legend_order <- c("sand","concrete")
  
  lvls <- intersect(legend_order, levels(factor(df$surface)))
  cols <- surf_cols[lvls]; names(cols) <- lvls
  lbls <- c(sand = "High Roughness", concrete = "Low Roughness")[lvls]
  
  #Compute stacking index so first legend item goes on top
  rank_idx  <- match(lvls, legend_order)
  stack_idx <- length(lvls) - rank_idx
  
  y_rng <- if (is.null(y_range)) range(df$velocity_norm, na.rm = TRUE) else y_range
  x_sym <- ensym(x)
  
  ##ρ (correlation) ----
  rhos <- NULL
  if (rho_mode == "auto") {
    rhos <- sapply(lvls, function(s) {
      d  <- df[df$surface == s, , drop = FALSE]
      xv <- dplyr::pull(d, !!x_sym)
      if (length(unique(xv)) < 2 || length(unique(d$velocity_norm)) < 2) return(NA_real_)
      suppressWarnings(cor(xv, d$velocity_norm, method = rho_method, use = "complete.obs"))
    })
    names(rhos) <- lvls
  } else if (rho_mode == "manual") {
    rhos <- rho_vals[lvls]
  }
  
  ##β (slope from lm) ----
  betas <- NULL
  if (beta_mode == "auto") {
    betas <- sapply(lvls, function(s) {
      d <- df[df$surface == s, , drop = FALSE]
      if (nrow(d) < 2) return(NA_real_)
      fit <- stats::lm(stats::reformulate(deparse(x_sym), response = "velocity_norm"), data = d)
      unname(stats::coef(fit)[2])  #slope for x
    })
    names(betas) <- lvls
  } else if (beta_mode == "manual") {
    betas <- beta_vals[lvls]
  }
  
  xr <- range(dplyr::pull(df, !!x_sym), na.rm = TRUE)
  yr <- y_rng
  
  #Annotation (bottom-right, stacked)
  ann_df <- NULL
  if (!is.null(rhos) || !is.null(betas)) {
    xr <- range(dplyr::pull(df, !!x_sym), na.rm = TRUE)
    yr <- y_rng
    ann_df <- tibble::tibble(
      surface = lvls,
      label = paste0("\u03b2 = ",
                     ifelse(is.null(betas) | is.na(betas), "NA",
                            sprintf(paste0("%.", beta_digits, "f"), betas)),
                     ", \u03c1 = ",
                     ifelse(is.null(rhos)  | is.na(rhos),  "NA",
                            sprintf(paste0("%.", rho_digits,  "f"), rhos))),
      x = xr[2] - 0.02 * diff(xr),
      y = yr[1] + (stack_idx) * rho_gap * diff(yr) + 0.02 * diff(yr)
    )
  }
  
  ##Plot
  p <- ggplot(df, aes(x = !!x_sym, y = velocity_norm)) +
    geom_smooth(aes(color = surface, fill = surface),
                method = "lm", se = TRUE, linewidth = 1) +
    scale_color_manual(values = cols, breaks = lvls, labels = lbls, name = "Surface") +
    scale_fill_manual(values  = scales::alpha(cols, 0.25),
                      breaks = lvls, labels = lbls, name = "Surface") +
    coord_cartesian(ylim = y_rng) +
    labs(x = x_lab, y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  
  if (!is.null(ann_df)) {
    p <- p + geom_text(
      data = ann_df,
      aes(x = x, y = y, label = label, color = surface),
      inherit.aes = FALSE, hjust = 1, vjust = 0, fontface = "bold",
      size = rho_size, show.legend = FALSE
    )
  }
  
  p
}


#Legend extraction and shared-axis label trick

common_y <- range(shape_data$velocity_norm, na.rm = TRUE)

pA <- plot_shape_overlay(shape_data, elongation,  "Elongation (Width/Length)",
                         y_range = common_y, title = "A",
                         rho_mode = "auto", beta_mode = "auto")
pB <- plot_shape_overlay(shape_data, flatness,    "Flatness (Height/Width)",
                         y_range = common_y, title = "B",
                         rho_mode = "auto", beta_mode = "auto")
pC <- plot_shape_overlay(shape_data, Sphericity,  "Sphericity",
                         y_range = common_y, title = "C",
                         rho_mode = "auto", beta_mode = "auto")
pD <- plot_shape_overlay(shape_data, Convexity,   "Convexity",
                         y_range = common_y, title = "D",
                         rho_mode = "auto", beta_mode = "auto")

#Extract a single legend to place as a right column
legend_g <- cowplot::get_legend(pA)

#Remove individual legends and y-axis title
strip_y <- theme(axis.title.y = element_blank(),
                 axis.title.y.right = element_blank())

pA_n <- pA + strip_y + theme(legend.position = "none")
pB_n <- pB + strip_y + theme(legend.position = "none")
pC_n <- pC + strip_y + theme(legend.position = "none")
pD_n <- pD + strip_y + theme(legend.position = "none")

#Place new single y-axis title
y_title <- ggplot() +
  annotate("text", x = 0, y = 0.5,
           label = expression(V[norm]~" ("~V[particle]~"/"~bar(V)[size~class]~")"),
           angle = 90, vjust = 0.5, fontface = "bold") +
  theme_void() +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))

#Stitch plots
panel_grid <- (pA_n | pB_n) / (pC_n | pD_n)

final <- (y_title | panel_grid | legend_g) +
  patchwork::plot_layout(widths = c(0.08, 1, 0.30))  #Adjust to add/remove gap

#View
final

#Save
ggsave("figures/shape_factors.png", final, width = 7, height = 5, dpi = 600, bg = "white")





