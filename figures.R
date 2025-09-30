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

##Plot functions ----
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






