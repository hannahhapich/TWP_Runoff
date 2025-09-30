#Load libraries ----
library(tidyverse)
library(viridis)
library(patchwork)

#Load data ----
flux_data <- read.csv("data/output_data/flux_data.csv")
model_params <- read.csv("data/output_data/muthumasy_model_parameters.csv")
mass_cumulative <- read.csv("data/output_data/mass_cumulative.csv")
q50_mass_flux <- read.csv("data/output_data/quartile_mass_flux.csv")
metadata <- read.csv("data/input_data/metadata.csv")
runoff_data <- read.csv("data/output_data/runoff_data_vol.csv")
PSA_velocity <- read.csv("data/output_data/PSA_velocity.csv")
PSA_perc_mobil <- read.csv("data/output_data/PSA_percent_mobilized.csv")

#Wash off model plots ----
#Facet order + labels
params_surface <- model_params %>%
  mutate(
    rainfall_f = factor(rainfall, levels = c("high","med","low")),
    slope_f    = factor(slope,    levels = c(5,10,15))
  )

facet_labels <- list(
  slope_f    = c(`5`="Slope : 5º", `10`="Slope : 10º", `15`="Slope : 15º"),
  rainfall_f = c(high="Rainfall : High", med="Rainfall : Medium", low="Rainfall : Low")
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
    left_join(params_surface %>% select(condition, k_prime, rse), by = "condition")
  
  annos <- annos %>%
    mutate(y_anno = if_else(condition == 4 | condition == 5 | condition == 6, y_anno + 0.5, y_anno),
           x_anno = if_else(condition == 4 | condition == 5 | condition == 6, x_anno - 10, x_anno)) %>%
    mutate(y_anno = if_else(condition == 1 | condition == 2 | condition == 3, y_anno - 0.15, y_anno)) %>%
    mutate(y_anno = if_else(condition == 16 | condition == 17 | condition == 18, y_anno + 0.5, y_anno),
           x_anno = if_else(condition == 16 | condition == 17 | condition == 18, x_anno - 10, x_anno))
  
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
                  label = paste0("k' = ", round(k_prime, 3),
                                 "\nRSE = ", format(rse, digits = 3, nsmall = 3))),
              hjust = 0,
              size = 4) +
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
    ggtitle(paste(surface_type, "surface")) +
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

slope_tick_labels <- c(`5` = "5º", `10` = "10º", `15` = "15º")
rain_tick_labels  <- c(low = "Low", med = "Medium", high = "High")

#Function to make one 3×3 heatmap figure (faceted by surface rows)
make_heatmap <- function(df, value_col, title, limits, palette = "inferno", direction = 1) {
  ggplot(df, aes(x = slope_f, y = rainfall_f, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.7) +
    facet_grid(rows = vars(surface_f)) +
    scale_x_discrete(name = "Slope",    labels = c(`5`="5º", `10`="10º", `15`="15º")) +
    scale_y_discrete(name = "Rainfall", labels = c(low="Low", med="Medium", high="High")) +
    scale_fill_viridis_c(option = palette, limits = limits, name = title, direction = direction) +  # <- many colors
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          strip.text.y = element_text(face = "bold"),
          axis.title.x = element_text(margin = margin(t = 6)),
          axis.title.y = element_text(margin = margin(r = 8)),
          legend.title = element_text(face = "bold"))
}

#Set limits
kp_lim <- range(params$k_prime, na.rm = TRUE)
q50_lim <- range(q50_mass_flux$Q50_time_min_mean, na.rm = TRUE)
max_frac_lim <- range(q50_mass_flux$max_frac, na.rm = TRUE)

#Make individual plots
p_kp <- make_heatmap(params, "k_prime","k′",   kp_lim, palette = "magma")
p_q50 <- make_heatmap(q50_mass_flux, "Q50_time_min_mean","Q50 (min)",   q50_lim, palette = "magma", direction = -1)
p_max_frac <- make_heatmap(q50_mass_flux, "max_frac","Fw",   max_frac_lim, palette = "magma")

#Arrange side-by-side
heatmaps <- (p_kp + theme(legend.position = "right")) +
  (p_q50 + theme(legend.position = "right")) +
  (p_max_frac + theme(legend.position = "right")) +
  plot_layout(ncol = 3, widths = c(1,1,1), guides = "keep")

#View
heatmaps

#Save figs
ggsave("figures/heatmaps.png", heatmaps, width = 9, height = 3, dpi = 600)
#ggsave("figures/heatmaps.pdf",  heatmaps, width = 10, height = 8)



#PSA Velocity Quartiles----
size_levels <- c("<125 µm","125-250 µm","250-500 µm","500-1000 µm",">1000 µm")
size_cols   <- setNames(
  magma(length(size_levels), begin = 0.05, end = 0.85), #Avoid near-black & near-white
  size_levels
)

PSA_clean <- PSA_velocity %>%
  mutate(
    size_class = str_replace_all(size_class, "–", "-"),
    size_f     = factor(size_class, levels = size_levels),
    rainfall_f = factor(rainfall, levels = c("low","med","high")),
    slope_f    = factor(slope, levels = c(5,10,15)),
    surface_f  = factor(surface, levels = c("sand","concrete"))
  )

#Plot function
plot_psa <- function(df, kind = c("count","volume")) {
  kind <- match.arg(kind)
  
  if (kind == "count") {
    q25 <- sym("Q25_count_time"); q50 <- sym("Q50_count_time"); q75 <- sym("Q75_count_time")
    title_lab <- "Particle half load by size ("~Q[25]~", "~Q[50]~", and"~Q[75]~", by count)"
  } else {
    q25 <- sym("Q25_vol_time");   q50 <- sym("Q50_vol_time");   q75 <- sym("Q75_vol_time")
    title_lab <- "Particle half load by size ("~Q[25]~", "~Q[50]~", and"~Q[75]~", by volume)"
  }
  
  pos <- position_dodge(width = 0.6)  #Controls left-right separation within each size bin
  
  ggplot(df, aes(x = size_f, y = !!q50,
                 color = size_f,         #Same color for the three rainfall points within a size bin
                 shape = rainfall_f,     #Diamond / circle / square
                 group = rainfall_f)) +
    geom_errorbar(aes(ymin = !!q25, ymax = !!q75),
                  position = pos, width = 0.2, linewidth = 0.7) +
    geom_point(position = pos, size = 3) +
    facet_grid(rows = vars(surface_f), cols = vars(slope_f)) +
    scale_color_manual(values = size_cols, name = "Size") +
    scale_shape_manual(values = c(low = 18, med = 16, high = 15), name = "Rainfall") +
    labs(x = "Particle size (µm)", y = "Time (min)", title = title_lab) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1)
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
size_levels <- c("<125 µm","125-250 µm","250-500 µm","500-1000 µm",">1000 µm")
size_cols   <- setNames(magma(length(size_levels), begin = 0.05, end = 0.85), size_levels)

PSA_perc_clean <- PSA_perc_mobil %>%
  left_join(metadata %>% select(-c(run, replicate)) %>% distinct(), by = "condition") %>%
  mutate(
    size_class = str_replace_all(size_class, "–", "-"),
    size_f     = factor(size_class, levels = size_levels),
    slope_f    = factor(slope,    levels = c(5, 10, 15)),
    surface_f  = factor(surface,  levels = c("sand", "concrete")),
    rainfall_f = factor(rainfall, levels = c("low", "med", "high"))  # L→R inside bins
  )

plot_perc_mobil_grouped <- function(df, kind = c("count","volume"), diamond_gap = 1) {
  kind <- match.arg(kind)
  value_sym <- if (kind == "count") rlang::sym("perc_count_mobilized") else rlang::sym("perc_volume_mobilized")
  title_lab <- if (kind == "count") "Percent mobilized (by count)" else "Percent mobilized (by volume)"
  
  bar_width <- 0.7
  pos <- position_dodge2(width = 0.7, padding = 0.2, preserve = "single")  #Dodge by rainfall (3 bars per size bin)
  
  ggplot(df, aes(x = size_f, y = !!value_sym, fill = size_f, group = rainfall_f)) +
    geom_col(width = bar_width, position = pos, color = "white", linewidth = 0.3) +
    #Shape above each bar (same color as bar, shape encodes rainfall)
    geom_point(aes(y = !!value_sym + diamond_gap, color = size_f, shape = rainfall_f),
               position = pos, size = 2.25, show.legend = TRUE) +
    facet_grid(rows = vars(surface_f), cols = vars(slope_f)) +
    scale_fill_manual(values = size_cols, name = "Size") +
    scale_color_manual(values = size_cols, guide = "none") +
    scale_shape_manual(values = c(low = 18, med = 16, high = 15), name = "Rainfall") +
    labs(x = "Particle size (µm)", y = "Percent mobilized (%)", title = title_lab) +
    coord_cartesian(ylim = c(0, 100 + diamond_gap + 1)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

#Build figures
p_mobil_count  <- plot_perc_mobil_grouped(PSA_perc_clean, "count",  diamond_gap = 7)
p_mobil_volume <- plot_perc_mobil_grouped(PSA_perc_clean, "volume", diamond_gap = 7)

#View
p_mobil_count
p_mobil_volume

#Save
ggsave("figures/psa_mobilized_count.png", p_mobil_count, width = 9, height = 6, dpi = 600, bg = "white")
#ggsave("figures/psa_quartile_count.pdf",  p_mobil_count, width = 9, height = 6)

ggsave("figures/psa_mobilized_vol.png", p_mobil_volume, width = 9, height = 6, dpi = 600, bg = "white")
#ggsave("figures/psa_quartile_vol.pdf",  p_mobil_volume, width = 9, height = 6)









