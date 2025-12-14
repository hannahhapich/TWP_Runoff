#Load libraries ----
library(tidyverse)
library(broom)
library(minpack.lm)

#Load data ----
runoff_data <- read.csv("data/input_data/runoff_data.csv")
mass_data <- read.csv("data/input_data/mass_data.csv") #NOTE: sample_mass_g is already blank subtracted
metadata <- read.csv("data/input_data/metadata.csv")
rain_data <- read.csv("data/input_data/rain_data.csv")

#Convert runoff mass to volume -----
#Jones (1992) function: density of air-saturated water
jones_density <- function(temp_C) {
  density_kg_m3 <- 999.84847 +
    6.337563e-2 * temp_C -
    8.523829e-4 * temp_C^2 +
    6.943248e-6 * temp_C^3 -
    3.821216e-8 * temp_C^4
  
  #Convert to g/mL
  density_g_mL <- density_kg_m3 / 1000
  return(density_g_mL)
}

#Calculate density
runoff_data$density_g_mL <- jones_density(runoff_data$water_temp_C)

#Replace missing density rows with mean of available values
mean_density <- mean(runoff_data$density_g_mL, na.rm = TRUE)
runoff_data$density_g_mL[is.na(runoff_data$density_g_mL)] <- mean_density

#Calculate volume (mL)
runoff_data$volume_mL <- runoff_data$water_mass_g / runoff_data$density_g_mL

#Join with condition
runoff_data <- runoff_data %>% 
  left_join(metadata %>% select(condition, run) %>% distinct(), by = "run")

##Approximate runoff steady state ----
tol   <- 50 #mL tolerance
t_min <- 1 #search window (min)
t_max <- 5

runoff_data2 <- runoff_data %>%
  mutate(end_time_min = as.numeric(end_time_min)) %>%
  arrange(condition, run, end_time_min)

find_ss_time <- function(df, tol = 50, t_min = 2, t_max = 4) {
  cands <- df %>% filter(end_time_min >= t_min, end_time_min <= t_max) %>% pull(end_time_min)
  for (t0 in cands) {
    tail_vals <- df %>% filter(end_time_min >= t0) %>% pull(volume_mL)
    mu <- mean(tail_vals, na.rm = TRUE)
    if (all(abs(tail_vals - mu) <= tol)) return(tibble(steady_state_min = t0))
  }
  tibble(steady_state_min = NA_real_)
}

ss_by_run <- runoff_data2 %>%
  group_by(condition, run) %>%
  group_modify(~find_ss_time(.x, tol, t_min, t_max)) %>%
  ungroup()

ss_by_condition <- ss_by_run %>%
  group_by(condition) %>%
  summarise(steady_state_min_mean = mean(steady_state_min, na.rm = TRUE),
            n_runs_detected     = sum(!is.na(steady_state_min)),
            .groups = "drop")

ss_by_condition #Max = 2.5 min


#Blank calculations ----
#Filter data
blank_mass <- mass_data %>% filter(end_time_min == "Blank")
sample_mass <- mass_data %>% filter(!start_time_min == "N/A") %>%
  mutate(raw_mass_g = jar_and_sample_mass_g - jar_mass_g) #Correction for % calculation because sample_mass_g is already blank subtracted

#Mean and SD blank mass, % total sample mass
print(mean(blank_mass$sample_mass_g)) #mean = 0.005259625 g
print(sd(blank_mass$sample_mass_g)) #SD = 0.006978204 g
print(mean(blank_mass$sample_mass_g)/mean(sample_mass$raw_mass_g)) #blank % of avg sample mass = 4.06%
print(sd(blank_mass$sample_mass_g)/mean(sample_mass$raw_mass_g)) #blank SD % of avg sample mass = 5.39%
print(min(blank_mass$sample_mass_g)) #min blank mass = 0.00069
print(max(blank_mass$sample_mass_g)) #max blank mass = 0.05441
print(min(sample_mass$raw_mass_g)) #min sample mass = 0.00364
print(max(sample_mass$raw_mass_g)) #max sample mass = 0.78407

#Variability between replicates ----
##Mass-based variability between individual samples ----

#Mean, SD, and CV
mass_variability <- mass_data %>%
  group_by(end_time_min, slope, surface, rainfall) %>%
  summarize(
    n = n(),
    mean = mean(sample_mass_g, na.rm = TRUE),
    sd   = sd(sample_mass_g, na.rm = TRUE),
    cv   = (sd / mean) * 100,
    .groups = "drop"
  ) %>%
  arrange(slope, surface, rainfall, end_time_min)

print(mean(mass_variability$cv)) #Mean CV = 18.3279%

#Filtered variability
mass_variability_filtered <- mass_variability %>% 
  filter (end_time_min != "30", end_time_min != "25", end_time_min != "Blank", end_time_min != "Remainder")

print(mean(mass_variability_filtered$cv)) #Mean CV = 10.672%

###KS tests for mass flux distributions ----
#Mass based KS test function
ks_mass <- function(df) {
  #Get all pairwise combinations of replicates
  replicates <- unique(df$replicate)
  combs <- combn(replicates, 2, simplify = FALSE)
  
  #Run ks.test for each replicate pair
  map_dfr(combs, function(pair) {
    x <- df$sample_mass_g[df$replicate == pair[1]]
    y <- df$sample_mass_g[df$replicate == pair[2]]
    ks <- ks.test(x, y)
    tibble(
      replicate1 = pair[1],
      replicate2 = pair[2],
      D = ks$statistic,
      p_value = ks$p.value
    )
  })
}

#Apply function
mass_ks <- mass_data %>%
  group_by(slope, surface, rainfall) %>%
  group_modify(~ ks_mass(.x)) %>%
  ungroup()

print(mean(mass_ks$p_value)) #Mean p_value = 0.9898544
print(sd(mass_ks$p_value)) #SD p_value = 0.03340815

print(mean(mass_ks$D)) #Mean D = 0.1712963
print(sd(mass_ks$D)) #SD D = 0.06300634




#Average precipitation rates----
#Rain gauge data
average_rain <- rain_data %>%
  group_by(rainfall) %>%
  summarize(
    n = n(),
    mean_mm_min = mean(gauge_avg_mm_min, na.rm = TRUE),
    sd   = sd(gauge_avg_mm_min, na.rm = TRUE)
  ) %>%
  mutate(mean_mm_hr = mean_mm_min * 60,
         sd_hr = sd * 60)

print(average_rain) #Average rain gauge data

#Runoff calculated water flux
water_flux <- runoff_data %>%
  filter(end_time_min > 5) %>%
  group_by(run) %>%
  summarize(mean_mL_min = (mean(volume_mL))) %>%
  left_join(metadata %>% select(run, slope, rainfall), by = "run") %>%
  mutate(mean_mL_min = ifelse(rainfall == "high", mean_mL_min * 2, mean_mL_min)) %>% #Correcting for high rainfall runoff collections being every 0.5 mins
  mutate(slope_rad = slope * (pi / 180)) %>%
  mutate(projected_area = (7200)*(cos(slope_rad))) %>% #Choi (2017) equation for projected (rain) area, 120 x 60 cm surface = 7200 cm2
  mutate(avg_mm_min = (mean_mL_min/projected_area)*10) %>%
  left_join(metadata %>% select(condition, run) %>% distinct(), by = "run")

average_water_flux <- water_flux %>%
  group_by(rainfall) %>%
  summarize(
    n = n(),
    mean_mm_min = mean(avg_mm_min, na.rm = TRUE),
    sd   = sd(avg_mm_min, na.rm = TRUE)
  ) %>%
  mutate(mean_mm_hr = mean_mm_min * 60,
         sd_hr = sd * 60)

print(average_water_flux) #Average rain data from water flux


#Estimating flow depth ----
#Calculate Re to determine laminar vs turbulent flow
v = 1.0*10^(-6) #Kinematic viscosity (1.0×10^−6 m2/s)
water_flux <- water_flux %>% 
  mutate(q = ((mean_mL_min * 10^(-6))/60)/0.6) %>% #Unit discharge: convert mL/min -> m3/s, divide by flow width (0.6m)
  mutate(Re = q/v) #Re = unit discharge (m2/s) over kinematic viscosity (1.0×10^−6 m2/s)

print(min(water_flux$Re)) #Min Re
print(max(water_flux$Re)) #Max Re

#Calculate mean density
density_mean_run <- runoff_data %>% group_by(run) %>% 
  summarize(density_mean_g_mL = mean(runoff_data$density_g_mL))

rho <- mean(density_mean_run$density_mean_g_mL) *1000 #kg/m^3 (average water density)

#Calculate flow depth with laminar overland flow equation
g = 9.81 #Acceleration due to gravity (9.81 m/s2)
water_flux <- water_flux %>%
  mutate(S = tan(slope_rad)) %>% #Dimensionless slope
  mutate(y = ((3*v*q)/(g*S))^(1/3)) %>% #Flow depth in m
  mutate(depth_mm = y*1000, #Flow depth converted to mm
         V = q/y) %>% #Flow velocity in m/s
  left_join(metadata %>% select(run, surface) %>% distinct(), by = "run") %>%
  mutate(tau = g * rho * y * S)
  



#Fitting wash-off models ----
##Data cleaning ----
#Calculate cumulative sample mass
sample_mass_cumulative <- sample_mass %>%
  group_by(run) %>%
  mutate(start_time_min = as.numeric (start_time_min),
         end_time_min = as.numeric(end_time_min)) %>%
  arrange(end_time_min, .by_group = TRUE) %>%
  left_join(metadata %>% select(run, condition), by = "run") %>% #Add condition ID
  mutate(cumulative_mass_g = cumsum(sample_mass_g)) %>%
  ungroup()

w0 = 1 #Initial load, 1 g
wash_off_mass <- sample_mass_cumulative %>% 
  select(-c(jar_mass_g, jar_and_sample_mass_g, raw_mass_g)) %>%
  left_join(water_flux %>% select(run, avg_mm_min), by = "run") %>% 
  rename(t = end_time_min, #t = time of wash-off measurement (min)
         i = avg_mm_min, #i = rainfall intensity (mm/min)
         wt = cumulative_mass_g) %>% #Wt = total wash-off mass (g) at time t
  mutate(frac = wt/w0) #Wash off fraction here equals Wt because W0 = 1 g, but here for future calculations

#Integrate across replicates
wash_off_avg <- wash_off_mass %>%
  group_by(condition, surface, t) %>%
  summarize(
    frac_mean = mean(frac),
    sd_frac = sd(frac),
    i_mean = mean(i),
    .groups = "drop"
  )

#Uncoupled 2-factor model (capacity factor and k) (Egodawatta et al, 2007) ----
##Check model for single condition ----
model_check <- nls( #Nonlinear least squares
  frac_mean ~ cf * (1 - exp(-k * i_mean * t)),
  data = wash_off_avg %>% filter(condition == 1),
  start = list(cf = 0.5, k = 0.01)  #Starting guesses
)
summary(model_check) #Capacity factor = 0.82, k = 0.091, RSE = 0.03 (good fit, within 3% predicted curve)

#Plot data points vs model fit for visual confirmation
df_cond1 <- wash_off_avg %>% filter(condition == 1)  # or however you're labeling condition 1
i_val <- unique(df_cond1$i_mean)

df_pred <- data.frame(
  t = seq(0, max(df_cond1$t), length.out = 100)
) %>%
  mutate(
    pred_frac = coef(model_check)["cf"] * (1 - exp(-coef(model_check)["k"] * i_val * t))
  )

ggplot(df_cond1, aes(x = t, y = frac_mean)) +
  geom_point(color = "blue") +
  geom_line(data = df_pred, aes(x = t, y = pred_frac), color = "red") +
  labs(x = "Time (min)", y = "Wash-off fraction",
       title = "Condition 1: observed vs fitted curve")


##Batch model fit ----
egodawatta_model_fit <- function(df_cond) {
  df_cond <- df_cond %>% arrange(t)
  i_val <- unique(df_cond$i_mean)[1]
  
  cf0 <- pmin(0.98, max(df_cond$frac_mean, na.rm = TRUE)) #Starting guess for cf: just under the observed max
  slope <- (df_cond$frac_mean[2] - df_cond$frac_mean[1]) /
    (df_cond$t[2] - df_cond$t[1])
  k0 <- max(1e-4, slope / (cf0 * i_val)) #Starting guess for k: slope ≈ cf*k*i  =>  k ≈ slope/(cf*i)
  
  nls(
    frac_mean ~ cf * (1 - exp(-k * i_val * t)),
    data = df_cond,
    start = list(cf = cf0, k = k0),
    control = nls.control(maxiter = 500, warnOnly = TRUE)
  )
}


#Apply function for all conditions
egodawatta_model_results <- map(unique(wash_off_avg$condition), function(cond) {
  df_cond <- filter(wash_off_avg, condition == cond)
  mod <- try(egodawatta_model_fit(df_cond), silent = TRUE)
  if (inherits(mod, "try-error")) return(NULL)
  list(condition = cond, model = mod)
})

#Extract parameter table
param_table <- map_dfr(egodawatta_model_results, function(x) {
  if (is.null(x)) return(NULL)
  df_cond <- dplyr::filter(wash_off_avg, condition == x$condition)
  mod <- x$model
  
  rss <- sum(residuals(mod)^2) #Residual sum of squares
  n   <- nrow(df_cond)
  p   <- length(stats::coef(mod)) #Capacity factor and k
  rse <- sqrt(rss / pmax(1, n - p)) #Residual standard error
  max_frac <- max(df_cond$frac_mean, na.rm = TRUE) #Mass wash-off fraction measured
  
  broom::tidy(mod) %>%
    select(term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    mutate(
      condition = x$condition,
      n = n,
      rss = rss,
      rse = rse,
      max_frac = max_frac
    )
})

#Model predictions
predictions <- map_dfr(egodawatta_model_results, function(x) {
  if (is.null(x)) return(NULL)
  df_cond <- filter(wash_off_avg, condition == x$condition)
  i_val <- unique(df_cond$i_mean)[1]
  t_seq <- seq(0, max(df_cond$t), length.out = 100)
  tibble(
    condition = x$condition,
    t = t_seq,
    pred_frac = coef(x$model)["cf"] * (1 - exp(-coef(x$model)["k"] * i_val * t_seq))
  )
})

#Plot data points vs model predictions for visual confirmation
ggplot() +
  geom_point(data = wash_off_avg, aes(x = t, y = frac_mean), alpha = 0.7, color = "blue") +
  geom_line(data = predictions, aes(x = t, y = pred_frac), color = "red") +
  facet_wrap(~ condition, scales = "free_x") +
  labs(x = "Time (min)", y = "Wash-off fraction",
       title = "Uncoupled 2-factor model: observed vs fitted by condition") +
  theme_minimal()


#Coupled model (f(k) and kprime) (Muthusamy et al, 2018) ----
#Build functions for iterative solving of c* and k'
rss_for_condition <- function(df_cond, cval, kprime) {
  i_val <- unique(df_cond$i_mean)[1]
  fcap  <- cval * kprime
  if (is.na(i_val) || cval <= 0 || kprime <= 0 || fcap > 1) return(Inf)
  pred <- fcap * (1 - exp(-kprime * i_val * df_cond$t))
  sum((df_cond$frac_mean - pred)^2, na.rm = TRUE)
}

fit_kprime_for_c <- function(df_cond, cval) {
  upper <- (1 / cval) - 1e-12 #c*k < 1 (small buffer added so it doesn't equal 1)
  if (upper <= 0) return(NULL)
  opt <- optimize(function(kp) rss_for_condition(df_cond, cval, kp),
                  interval = c(1e-8, upper))
  list(kprime = opt$minimum,
       rss    = opt$objective,
       fcap   = cval * opt$minimum)
}

#Solve for c and find min value (c*) for each surface
c_grid <- seq(1, 20, by = 0.2)

#For one surface's data frame, sweep c and sum RSS across its conditions
fit_for_c_one_surface <- function(df_surface, cval) {
  fits <- df_surface %>%
    dplyr::group_by(condition) %>%
    dplyr::group_map(~{
      f <- fit_kprime_for_c(.x, cval)
      if (is.null(f)) return(NULL)
      tibble::tibble(condition = .y$condition, kprime = f$kprime,
                     fcap = f$fcap, rss = f$rss)
    }) %>% dplyr::bind_rows()
  tibble::tibble(c = cval,
                 total_rss = sum(fits$rss),
                 n_cond = nrow(fits))
}

#Sweep c for each surface and pick the minimizing c* for that surface, save all data
rss_sweeps_by_surface <- wash_off_avg %>%
  dplyr::group_by(surface) %>%
  dplyr::group_map(~{
    df_surf <- .x
    purrr::map_dfr(c_grid, ~fit_for_c_one_surface(df_surf, .x)) %>%
      dplyr::mutate(surface = .y$surface[[1]])
  }) %>%
  dplyr::bind_rows()

#Extract the minimizing point (c* and its RSS) for each surface
c_star_by_surface <- rss_sweeps_by_surface %>%
  dplyr::group_by(surface) %>%
  dplyr::slice_min(total_rss, n = 1, with_ties = FALSE) %>%
  dplyr::rename(c_star = c, rss_min = total_rss)

print(c_star_by_surface) #c*(concrete) = 3.8; c*(sand) = 4.8

#Plot total RSS vs c, faceted by surface, with c* shown
p_rss_by_surface <- ggplot(rss_sweeps_by_surface,
                           aes(x = c, y = total_rss)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ surface, scales = "free_y") +
  geom_vline(data = c_star_by_surface,
             aes(xintercept = c_star),
             linetype = 2) +
  geom_point(data = c_star_by_surface,
             aes(x = c_star, y = rss_min),
             size = 2) +
  ggrepel::geom_label_repel(
    data = c_star_by_surface,
    aes(x = c_star, y = rss_min,
        label = paste0("c* = ", signif(c_star, 3))),
    min.segment.length = 0
  ) +
  labs(title = "Total RSS vs c (per surface)",
       x = "c (mm)",
       y = "Total residual sum-of-squares") +
  theme_bw()

#Print and save figure
print(p_rss_by_surface)
ggsave("figures/c_star_RSS.png", p_rss_by_surface, width = 7, height = 5, dpi = 600)


#Function to fit all conditions with surface specific c*
refit_at_cstar <- function(df_cond, cval) {
  f <- fit_kprime_for_c(df_cond, cval)
  if (is.null(f)) return(tibble::tibble())
  rss <- f$rss
  n   <- nrow(df_cond)
  p   <- 1L
  rse <- sqrt(rss / pmax(1, n - p))
  
  tibble::tibble(
    c_mm      = cval,
    k_prime   = f$kprime,
    f_k       = f$fcap,
    rss       = rss,
    n_points  = n,
    rse       = rse,
    max_frac  = max(df_cond$frac_mean, na.rm = TRUE)
  )
}

#Apply function and extract parameters
param_table_coupled <- wash_off_avg %>%
  dplyr::left_join(c_star_by_surface, by = "surface") %>%
  dplyr::group_by(condition, surface, c_star) %>%
  dplyr::group_modify(~ refit_at_cstar(.x, unique(.y$c_star))) %>%
  dplyr::ungroup()

#Add rain intensity per condition
param_table_coupled <- param_table_coupled %>%
  dplyr::left_join(
    wash_off_avg %>% filter(t == 30) %>% dplyr::select(condition, i_mean, sd_frac) %>% dplyr::distinct(),
    by = "condition"
  )

#Per-condition intensity & time horizon
cond_meta <- wash_off_avg %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(
    i_val = unique(i_mean)[1],
    t_max = max(t, na.rm = TRUE),
    .groups = "drop"
  )

#Extract model predictions
params <- param_table_coupled %>%
  dplyr::inner_join(cond_meta, by = "condition")

predictions_coupled <- params %>%
  purrr::pmap_dfr(function(condition, c_mm, k_prime, f_k, rss, n_points, rse,
                           max_frac, i_mean, t_max, surface, c_star, ...) {
    t_seq <- seq(0, t_max, length.out = 100)
    tibble::tibble(
      condition = condition,
      surface   = surface,
      t         = t_seq,
      pred_frac = f_k * (1 - exp(-(k_prime * i_mean) * t_seq))
    )
  })

#Plot data points vs model predictions for visual confirmation
ggplot() +
  geom_point(data = wash_off_avg, aes(x = t, y = frac_mean), alpha = 0.7, color = "blue") +
  geom_line(data = predictions_coupled, aes(x = t, y = pred_frac), color = "red") +
  facet_wrap(~ condition, scales = "free_x") +
  labs(x = "Time (min)", y = "Wash-off fraction",
       title = "Coupled model (surface-specific c*): observed vs fitted by condition") +
  theme_minimal()



#Mass-based Q50 (from linear interpolation) ----
sample_mass_Q50 <- sample_mass_cumulative %>% 
  select(-c(jar_mass_g, jar_and_sample_mass_g, sample_mass_g, raw_mass_g))

#Calculate Q50 for each run
##Function
qtime_mass <- function(dat, time_col = "end_time_min",
                                cum_col = "cumulative_mass_g", p = 0.5) {
  t  <- dat[[time_col]]
  cM <- dat[[cum_col]]
  
  #Keep rows with finite time and cumulative
  ok <- is.finite(t) & is.finite(cM)
  t  <- t[ok]; cM <- cM[ok]
  if (!length(t)) return(NA_real_)
  o  <- order(t); t <- t[o]; cM <- cM[o]
  
  #Explicit (0,0) anchor
  if (t[1] > 0) { t <- c(0, t); cM <- c(0, cM) } else if (t[1] == 0) cM[1] <- 0
  
  #Final cumulative = value at the last time point in that run (30 if present, otherwise the max time)
  final_cum <- cM[length(cM)]; if (!is.finite(final_cum) || final_cum <= 0) return(NA_real_)
  
  #Define target mass
  target <- p * final_cum
  #First index where cumulative >= target
  if (max(cM, na.rm = TRUE) < target) return(NA_real_)
  
  #Linear interpolation between points
  k <- which(cM >= target)[1]; if (k == 1) return(t[1])
  t[k-1] + (target - cM[k-1])/(cM[k]-cM[k-1]) * (t[k]-t[k-1])
}

#Apply to all runs
mass_qtile <- sample_mass_Q50 %>%
  group_by(surface, slope, rainfall, run, condition) %>%
  summarise(
    Q25_time_min = qtime_mass(cur_data(), p = 0.25),
    Q50_time_min = qtime_mass(cur_data(), p = 0.50),
    Q75_time_min = qtime_mass(cur_data(), p = 0.75),
    .groups = "drop"
  )

#Summarize across groups
mass_qtile_summary <- mass_qtile %>%
  group_by(surface, slope, rainfall, condition) %>%
  summarize(
    n_runs = dplyr::n(),
    across(starts_with("Q"),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd   = ~sd(.x,   na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) 

#Max wash-off fraction
max_frac <- sample_mass_Q50 %>% filter(end_time_min == 30) %>% 
  select(cumulative_mass_g, run, condition) %>%
  group_by(condition) %>%
  summarize(
    max_frac = mean(cumulative_mass_g)
    )

mass_qtile_summary <- mass_qtile_summary %>% 
  left_join(max_frac, by = "condition")



#Summary stats ----
##Wash-off factors vs experimental conditions ----
df <- param_table_coupled %>% 
  left_join(mass_qtile_summary %>% select(condition, Q50_time_min_mean, Q50_time_min_sd), by = "condition") %>%
  left_join(metadata %>% select(condition, slope, rainfall) %>% distinct(), by = "condition") %>%
  select(-c(c_mm, n_points, i_mean))

#Compare high rain intensity Fw per surface for discussion section
df %>% select(max_frac, sd_frac, rainfall, surface) %>%
  filter(rainfall == "high") %>%
  group_by(surface) %>%
  summarize(avg_max_frac = mean(max_frac),
            avg_sd_frac = mean(sd_frac))

vars  <- c("k_prime","f_k","max_frac","Q50_time_min_mean","rse","sd_frac","Q50_time_min_sd")

#Single-factor groupings (summarized separately)
group_vars <- c("surface","rainfall","slope")

#Metrics to pull min/max groups for
minmax_targets <- c("k_prime","f_k","max_frac","Q50_time_min_mean")

#Averages by each individual factor
averages_by_factor <- map(group_vars, function(gv) {
  df %>%
    group_by(.data[[gv]]) %>%
    summarise(
      across(all_of(vars), ~ mean(.x, na.rm = TRUE), .names = "avg_{.col}"),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(.data[[gv]])
}) %>% set_names(group_vars)

#View individual tables
averages_by_factor$surface
averages_by_factor$rainfall
averages_by_factor$slope

#Min/max for each variable
pull_level_chr <- function(df, col) {
  as.character(df %>% dplyr::pull(dplyr::all_of(col)))
}

minmax_by_factor <- purrr::map_dfr(group_vars, function(gv) {
  avg_tbl <- averages_by_factor[[gv]]
  val_for <- function(metric) paste0("avg_", metric)
  
  purrr::map_dfr(minmax_targets, function(metric) {
    valcol <- val_for(metric)
    
    min_row <- avg_tbl %>% dplyr::slice_min(.data[[valcol]], n = 1, with_ties = FALSE)
    max_row <- avg_tbl %>% dplyr::slice_max(.data[[valcol]], n = 1, with_ties = FALSE)
    
    tibble::tibble(
      group     = gv,
      metric    = metric,
      min_level = pull_level_chr(min_row, gv),
      min_value = min_row %>% dplyr::pull(dplyr::all_of(valcol)),
      max_level = pull_level_chr(max_row, gv),
      max_value = max_row %>% dplyr::pull(dplyr::all_of(valcol))
    )
  })
})

#View
minmax_by_factor



#Significance testing ----
##Mass-based significance testing ----
###Model parameters ----
#Data prep
param_table_coupled <- param_table_coupled %>% 
  left_join(metadata %>% select(-c(run, replicate, surface)) %>% distinct(), by = "condition")

#Run ANOVA
#Ensure factors are coded properly
param_table_coupled <- param_table_coupled %>%
  mutate(
    slope = factor(slope),
    rainfall = factor(rainfall),
    surface = factor(surface)
  )

#ANOVA tests
anova_kprime <- aov(k_prime ~ surface + rainfall + slope +
                          surface:rainfall + surface:slope + rainfall:slope,
                        data = param_table_coupled)
anova_fk <- aov(f_k ~ surface + rainfall + slope +
                      surface:rainfall + surface:slope + rainfall:slope,
                    data = param_table_coupled)

summary(anova_kprime)
summary(anova_fk)

###Q50 and total wash-off fraction ----
#Data prep
mass_qtile_summary <- mass_qtile_summary %>% 
  mutate(
    slope = factor(slope),
    rainfall = factor(rainfall),
    surface = factor(surface)
  )

#ANOVA tests
anova_Q50 <- aov(Q50_time_min_mean ~ surface + rainfall + slope +
                      surface:rainfall + surface:slope + rainfall:slope,
                    data = mass_qtile_summary)
anova_Fw <- aov(max_frac ~ surface + rainfall + slope +
                  surface:rainfall + surface:slope + rainfall:slope,
                data = mass_qtile_summary)

summary(anova_Q50)
summary(anova_Fw)



#Data export ----
write.csv(water_flux, "data/output_data/flux_data.csv", row.names = F)
write.csv(param_table_coupled, "data/output_data/muthusamy_model_parameters.csv" , row.names = F)
write.csv(mass_qtile_summary, "data/output_data/quartile_mass_flux.csv" , row.names = F)
write.csv(sample_mass_cumulative, "data/output_data/mass_cumulative.csv" , row.names = F)
write.csv(runoff_data, "data/output_data/runoff_data_vol.csv" , row.names = F)









#PSA Calculations ----
#Note: PSA calculations separated due to data size
##Load data ----
#WARNING, will take a few mins...
PSA_data <- read.csv("data/input_data/PSA_data.csv")

##Blank calculations ----
#Filter data
blank_PSA <- PSA_data %>% filter(interval == "0")

#Mean and SD blank mass, % total sample mass
print(median(blank_PSA$FLength)) #D50 = 113.229
print(sd(blank_PSA$FLength)) #SD = 85.56953

#Average mobilized particle size
mobilized_PSA <- PSA_data %>% filter(interval < 7, interval > 0)
print(median(mobilized_PSA$FLength)) #D50 = 161.014
print(sd(mobilized_PSA$FLength)) #SD = 177.1187

##Variability between replicates ----
###KS tests for PSA distributions ----
#Define size classes
size_breaks <- c(-Inf, 125, 250, 500, 1000, Inf)
size_labels <- c("<125 µm", "125–250 µm", "250–500 µm", "500–1000 µm", ">1000 µm")

mobilized_PSA <- PSA_data %>% filter(interval > 0, interval < 7) %>%
  left_join(metadata %>% select(c(condition, run)) %>% distinct(), by = "run")

#Assign size classes from above
mobilized_size_class <- mobilized_PSA %>%
  mutate(size_class = cut(FLength, breaks = size_breaks, labels = size_labels))

#Summarize mobilized
mob_sum <- mobilized_size_class %>%
  group_by(run, size_class) %>%
  summarize(mobilized_count  = n(),
            mobilized_volume = sum(Volume, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(metadata %>% select(condition, run), by = "run")

#PSA KS test functions
ks_psa_count <- function(df) {
  replicates <- unique(df$run)
  combs <- combn(replicates, 2, simplify = FALSE) #Get pairwise combinations of replicates
  
  map_dfr(combs, function(pair) { #Run ks.test for each replicate pair
    x <- df$mobilized_count[df$run == pair[1]]
    y <- df$mobilized_count[df$run == pair[2]]
    ks <- ks.test(x, y)
    tibble(
      replicate1 = pair[1],
      replicate2 = pair[2],
      D = ks$statistic,
      p_value = ks$p.value)
  })
}

ks_psa_vol <- function(df) {
  replicates <- unique(df$run)
  combs <- combn(replicates, 2, simplify = FALSE) #Get pairwise combinations of replicates
  
  map_dfr(combs, function(pair) { #Run ks.test for each replicate pair
    x <- df$mobilized_volume[df$run == pair[1]]
    y <- df$mobilized_volume[df$run == pair[2]]
    ks <- ks.test(x, y)
    tibble(
      replicate1 = pair[1],
      replicate2 = pair[2],
      D = ks$statistic,
      p_value = ks$p.value)
  })
}

#Apply functions
count_ks <- mob_sum %>%
  group_by(condition) %>%
  group_modify(~ ks_psa_count(.x)) %>%
  ungroup()

print(mean(count_ks$p_value)) #Mean p_value = 0.9247501
print(sd(count_ks$p_value)) #SD p_value = 0.06297975

print(mean(count_ks$D)) #Mean D = 0.3185185
print(sd(count_ks$D)) #SD D = 0.09919311

vol_ks <- mob_sum %>%
  group_by(condition) %>%
  group_modify(~ ks_psa_vol(.x)) %>%
  ungroup()

print(mean(vol_ks$p_value)) #Mean p_value = 0.9929453
print(sd(vol_ks$p_value)) #SD p_value = 0.02936029

print(mean(vol_ks$D)) #Mean D = 0.2111111
print(sd(vol_ks$D)) #SD D = 0.04624246


##Particle transport speed by size ----
#Function to derive quantile times
qtime <- function(dat, weight_col = c("count","volume"), p = 0.5) {
  weight_col <- rlang::arg_match(weight_col)
  w <- dat[[weight_col]]
  tot <- sum(w, na.rm = TRUE)
  if (length(w) == 0 || tot == 0) return(NA_real_)
  
  #Cumulative BEFORE interval
  cum_prev <- dplyr::lag(cumsum(w), default = 0)
  target   <- p * tot
  
  #Identify first interval where target quantile is crossed
  idx <- which(cum_prev < target & (cum_prev + w) >= target)
  if (!length(idx)) return(max(dat$end_time_min, na.rm = TRUE))
  i <- idx[1]
  
  t0 <- dat$start_time_min[i]
  t1 <- dat$end_time_min[i]
  wi <- w[i]
  if (wi <= 0) return(t1)
  
  #Linear extrapolation between bins
  frac <- (target - cum_prev[i]) / wi
  t0 + frac * (t1 - t0)
}

###Calculate count and volume based quantile times ----
size_breaks <- c(-Inf, 125, 250, 500, 1000, Inf)
size_labels <- c("<125 µm", "125–250 µm", "250–500 µm", "500–1000 µm", ">1000 µm")

#Calculate V50 based off of pooled particle data across all runs (per condition)
pooled_qtiles <- mobilized_PSA %>%
  mutate(size_class = cut(FLength, breaks = size_breaks, labels = size_labels, right = TRUE)) %>%
  filter(!is.na(size_class)) %>%
  group_by(condition, size_class, end_time_min) %>%
  summarize(
    count  = n(),
    volume = sum(Volume, na.rm = TRUE),
    .groups = "drop_last"   #Leaves grouping by condition & size_class
  ) %>%
  arrange(condition, size_class, end_time_min) %>%
  group_by(condition, size_class) %>%
  mutate(
    start_time_min = dplyr::lag(end_time_min, default = 0)
  ) %>%
  
  summarize( #Note: all quantile times here are in minutes
    Q25_count_time = qtime(cur_data_all(), "count", 0.25),
    Q50_count_time = qtime(cur_data_all(), "count", 0.50),
    Q75_count_time = qtime(cur_data_all(), "count", 0.75),
    Q25_vol_time   = qtime(cur_data_all(), "volume", 0.25),
    Q50_vol_time   = qtime(cur_data_all(), "volume", 0.50),
    Q75_vol_time   = qtime(cur_data_all(), "volume", 0.75),
    .groups = "drop"
  ) %>%
  arrange(condition, factor(size_class, levels = size_labels))

#Calculate per-run to obtain SD
run_qtiles <- mobilized_PSA %>%
  mutate(size_class = cut(FLength, breaks = size_breaks, labels = size_labels, right = TRUE)) %>%
  filter(!is.na(size_class)) %>%
  group_by(condition, run, size_class, end_time_min) %>%
  summarise(
    count  = n(),
    volume = sum(Volume, na.rm = TRUE),
    .groups = "drop_last"      #Leaves condition, run, size_class
  ) %>%
  arrange(condition, run, size_class, end_time_min) %>%
  group_by(condition, run, size_class) %>%
  mutate(start_time_min = dplyr::lag(end_time_min, default = 0)) %>%
  summarise(  #Per-run qtimes
    Q25_count_time = qtime(cur_data_all(), "count",  0.25),
    Q50_count_time = qtime(cur_data_all(), "count",  0.50),
    Q75_count_time = qtime(cur_data_all(), "count",  0.75),
    Q25_vol_time   = qtime(cur_data_all(), "volume", 0.25),
    Q50_vol_time   = qtime(cur_data_all(), "volume", 0.50),
    Q75_vol_time   = qtime(cur_data_all(), "volume", 0.75),
    .groups = "drop"
  )

run_qtile_summary <- run_qtiles %>%
  group_by(condition, size_class) %>%
  summarise(
    n_runs = dplyr::n(),
    across(
      starts_with("Q"),
      list(mean = ~ mean(.x, na.rm = TRUE),
           sd   = ~ sd(.x,   na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>% select(condition, size_class,
               Q25_count_time_sd, Q50_count_time_sd, Q75_count_time_sd,
               Q25_vol_time_sd, Q50_vol_time_sd, Q75_vol_time_sd)

#Join pooled V50 calculations with SD
PSA_qtile_times <- pooled_qtiles %>%
  left_join(run_qtile_summary, by = c("condition", "size_class"))

#Avoid divide by zero
safe_cv <- function(sd, mean) ifelse(is.finite(sd) & is.finite(mean) & mean > 0, sd/mean, NA_real_)

###Calculate quantile effective velocity ----
PSA_qtile_speeds <- PSA_qtile_times %>% 
  #Calculate quantile speeds (m/s)
  mutate(Q25_count_v = 1.1/(Q25_count_time*60), #Converting from minutes to seconds
         Q50_count_v = 1.1/(Q50_count_time*60),
         Q75_count_v = 1.1/(Q75_count_time*60),
         Q25_vol_v   = 1.1/(Q25_vol_time*60),
         Q50_vol_v   = 1.1/(Q50_vol_time*60),
         Q75_vol_v   = 1.1/(Q75_vol_time*60)) %>%
  #Calculate ratio of mean to SD for quantile times and apply to quantile velocities 
  #(solution to scaling issue where small sd values grow exponentially when multiplying the same as mean)
  mutate(cv_Q25_count_time = safe_cv(Q25_count_time_sd, Q25_count_time),
         cv_Q50_count_time = safe_cv(Q50_count_time_sd, Q50_count_time),
         cv_Q75_count_time = safe_cv(Q75_count_time_sd, Q75_count_time),
         cv_Q25_vol_time   = safe_cv(Q25_vol_time_sd,   Q25_vol_time),
         cv_Q50_vol_time   = safe_cv(Q50_vol_time_sd,   Q50_vol_time),
         cv_Q75_vol_time   = safe_cv(Q75_vol_time_sd,   Q75_vol_time),
         
         Q25_count_v_sd = cv_Q25_count_time * Q25_count_v,
         Q50_count_v_sd = cv_Q50_count_time * Q50_count_v,
         Q75_count_v_sd = cv_Q75_count_time * Q75_count_v,
         Q25_vol_v_sd   = cv_Q25_vol_time   * Q25_vol_v,
         Q50_vol_v_sd   = cv_Q50_vol_time   * Q50_vol_v,
         Q75_vol_v_sd   = cv_Q75_vol_time   * Q75_vol_v) %>%
  left_join(metadata %>% select(-c(replicate, run)) %>% distinct(), by = "condition")



##Particle transport effectiveness by size ----
immobile_PSA <- PSA_data %>% filter(interval == 7) %>%
  left_join(metadata %>% select(c(condition, run)) %>% distinct(), by = "run")

#Assign size classes from above
mobilized_size_class <- mobilized_PSA %>%
  mutate(size_class = cut(FLength, breaks = size_breaks, labels = size_labels))
immobile_size_class <- immobile_PSA %>%
  mutate(size_class = cut(FLength, breaks = size_breaks, labels = size_labels))

#Check the number of replicates is constant before comparing aggregate data
replicate_check_mob <- mobilized_size_class %>%
  group_by(condition) %>%
  summarize(n_replicates = n_distinct(run), .groups = "drop")
replicate_check_immob <- immobile_size_class %>%
  group_by(condition) %>%
  summarize(n_replicates = n_distinct(run), .groups = "drop")

#Summarize mobilized
mob_sum <- mobilized_size_class %>%
  group_by(condition, size_class) %>%
  summarize(
    mobilized_count  = n(),
    mobilized_volume = sum(Volume, na.rm = TRUE),
    .groups = "drop"
  )

#Summarize immobile
immob_sum <- immobile_size_class %>%
  group_by(condition, size_class) %>%
  summarize(
    immobile_count  = n(),
    immobile_volume = sum(Volume, na.rm = TRUE),
    .groups = "drop"
  )

#Combine & calculate % mobilized
percent_mobilized_pooled <- full_join(mob_sum, immob_sum,
                               by = c("condition","size_class")) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(
    total_count  = mobilized_count + immobile_count,
    total_volume = mobilized_volume + immobile_volume,
    perc_count_mobilized  = ifelse(total_count > 0, mobilized_count / total_count * 100, NA),
    perc_volume_mobilized = ifelse(total_volume > 0, mobilized_volume / total_volume * 100, NA)
  ) %>%
  arrange(condition, size_class)


#Per-run % mobilized (to obtain SD for each condition)
mob_sum_run <- mobilized_size_class %>%
  group_by(condition, run, size_class) %>%
  summarize(
    mobilized_count  = n(),
    mobilized_volume = sum(Volume, na.rm = TRUE),
    .groups = "drop"
  )

immob_sum_run <- immobile_size_class %>%
  group_by(condition, run, size_class) %>%
  summarize(
    immobile_count  = n(),
    immobile_volume = sum(Volume, na.rm = TRUE),
    .groups = "drop"
  )

#Per-run % mobilized
percent_mobilized_runs <- full_join(mob_sum_run, immob_sum_run,
                                    by = c("condition","run","size_class")) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  mutate(
    total_count  = mobilized_count + immobile_count,
    total_volume = mobilized_volume + immobile_volume,
    perc_count_mobilized_run  = ifelse(total_count  > 0, mobilized_count  / total_count  * 100, NA_real_),
    perc_volume_mobilized_run = ifelse(total_volume > 0, mobilized_volume / total_volume * 100, NA_real_)
  )

#Summarize SD across conditions
percent_mobilized_runs_summary <- percent_mobilized_runs %>%
  group_by(condition, size_class) %>%
  summarize(
    n_runs = dplyr::n(),
    perc_count_mobilized_sd  = sd(perc_count_mobilized_run,  na.rm = TRUE),
    perc_volume_mobilized_sd = sd(perc_volume_mobilized_run, na.rm = TRUE),
    .groups = "drop"
  )

#Join
percent_mobilized <- percent_mobilized_pooled %>%
  left_join(percent_mobilized_runs_summary,
            by = c("condition","size_class"))



##Particle shape analysis ----
shape_data <- mobilized_size_class %>% 
  filter(size_class != "<125 µm", size_class != "125–250 µm") %>%
  mutate(elongation = FWidth/FLength,
         flatness = FThickness/FWidth,
         velocity = 1.10 / (end_time_min * 60), #v in m/s
         surface = ifelse(condition < 10, "sand", "concrete"))

#Separate particles on sand vs concrete
sand_shape <- shape_data %>% filter(condition < 10) 
conc_shape <- shape_data %>% filter(condition > 9) 

#Find mean particle velocities as a function of size
sand_shape_summary <- sand_shape %>%
  group_by(size_class) %>%
  summarize(
    mean_velocity  = mean(velocity)
  )
conc_shape_summary <- conc_shape %>%
  group_by(size_class) %>%
  summarize(
    mean_velocity  = mean(velocity)
  )

#Calculate normalized particle velocity: v_particle/v_mean (mean for a given size class + surface)
sand_shape <- sand_shape %>%
  left_join(sand_shape_summary, by = "size_class") %>%
  mutate(velocity_norm = velocity/mean_velocity)
conc_shape <- conc_shape %>%
  left_join(conc_shape_summary, by = "size_class") %>%
  mutate(velocity_norm = velocity/mean_velocity)

#Correlation tests between shape descriptors and normalized velocity
cor.test(sand_shape$elongation, sand_shape$velocity_norm, method = "spearman") #rho = 0.1156033
cor.test(sand_shape$flatness, sand_shape$velocity_norm, method = "spearman") #rho = 0.007444385
cor.test(sand_shape$Sphericity, sand_shape$velocity_norm, method = "spearman") #rho = 0.1136321
cor.test(sand_shape$Convexity, sand_shape$velocity_norm, method = "spearman") #rho = 0.06800014

cor.test(conc_shape$elongation, conc_shape$velocity_norm, method = "spearman") #rho = 0.09362008
cor.test(conc_shape$flatness, conc_shape$velocity_norm, method = "spearman") #rho = 0.1742604
cor.test(conc_shape$Sphericity, conc_shape$velocity_norm, method = "spearman") #rho = 0.1211792
cor.test(conc_shape$Convexity, conc_shape$velocity_norm, method = "spearman") #rho = 0.1635184


#Regression model
sand_shape_lm <- lm(velocity_norm ~ elongation + flatness + Sphericity + Convexity, data = sand_shape)
conc_shape_lm <- lm(velocity_norm ~ elongation + flatness + Sphericity + Convexity, data = conc_shape)

#R2 values
summary(sand_shape_lm)$r.squared #0.02311312
summary(conc_shape_lm)$r.squared #0.009892324

#Join shape data
shape_data <- rbind(sand_shape, conc_shape)

##Summary Stats ----
###V50 and mobilizable amount vs experimental conditions + particle size ----
#Data prep
df <- percent_mobilized %>%
  select(condition, size_class, perc_count_mobilized, perc_count_mobilized_sd, perc_volume_mobilized, perc_volume_mobilized_sd) %>%
  left_join(PSA_qtile_speeds %>% select(condition, size_class, Q50_count_v, Q50_vol_v, Q50_count_v_sd, Q50_vol_v_sd), by = c("condition", "size_class")) %>%
  left_join(metadata %>% select(-c(run, replicate)) %>% distinct(), by = "condition") %>%
  mutate(Q50_count_v = Q50_count_v * 100, #Convert m/s to cm/s
         Q50_vol_v = Q50_vol_v * 100,
         Q50_count_v_sd = Q50_count_v_sd * 100,
         Q50_vol_v_sd = Q50_vol_v_sd * 100,
         perc_count_mobilized = perc_count_mobilized/100, #Convert from percent to fraction
         perc_count_mobilized_sd = perc_count_mobilized_sd/100,
         perc_volume_mobilized = perc_volume_mobilized/100,
         perc_volume_mobilized_sd = perc_volume_mobilized_sd/100)

vars  <- c("Q50_count_v", "Q50_vol_v", "Q50_count_v_sd", "Q50_vol_v_sd", 
           "perc_count_mobilized","perc_volume_mobilized", "perc_count_mobilized_sd", "perc_volume_mobilized_sd")

vars_sd   <- grep("_sd$", vars, value = TRUE)
vars_mean <- setdiff(vars, vars_sd)
#Single-factor groupings you want summarized separately
group_vars <- c("surface","rainfall","slope", "size_class")

#Metrics for which you want min/max groups called out
minmax_targets <- c("Q50_count_v", "Q50_vol_v", "perc_count_mobilized","perc_volume_mobilized")

#Averages by each individual factor
averages_by_factor <- purrr::map(group_vars, function(gv) {
  df %>%
    dplyr::group_by(.data[[gv]]) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(vars),
        ~ mean(.x, na.rm = TRUE),
        .names = "avg_{.col}"
      ),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data[[gv]])
}) %>% purrr::set_names(group_vars)


#View each table
averages_by_factor$surface
averages_by_factor$rainfall
averages_by_factor$slope
averages_by_factor$size_class

#Min/max for each variable
pull_level_chr <- function(df, col) {
  as.character(df %>% dplyr::pull(dplyr::all_of(col)))
}

minmax_by_factor <- purrr::map_dfr(group_vars, function(gv) {
  avg_tbl <- averages_by_factor[[gv]]
  val_for      <- function(metric) paste0("avg_", metric)           #Mean column
  val_for_sd   <- function(metric) paste0("avg_", metric, "_sd")    #Paired SD col (if it exists)
  
  purrr::map_dfr(minmax_targets, function(metric) {
    valcol   <- val_for(metric)
    sdcol    <- val_for_sd(metric)
    has_sd   <- sdcol %in% names(avg_tbl)
    
    min_row <- avg_tbl %>% dplyr::slice_min(.data[[valcol]], n = 1, with_ties = FALSE)
    max_row <- avg_tbl %>% dplyr::slice_max(.data[[valcol]], n = 1, with_ties = FALSE)
    
    tibble::tibble(
      group       = gv,
      metric      = metric,
      min_level   = pull_level_chr(min_row, gv),
      min_value   = min_row %>% dplyr::pull(dplyr::all_of(valcol)),
      min_value_sd = if (has_sd) min_row %>% dplyr::pull(dplyr::all_of(sdcol)) else NA_real_,
      max_level   = pull_level_chr(max_row, gv),
      max_value   = max_row %>% dplyr::pull(dplyr::all_of(valcol)),
      max_value_sd = if (has_sd) max_row %>% dplyr::pull(dplyr::all_of(sdcol)) else NA_real_
    )
  })
})

#View
minmax_by_factor



##Significance testing ----
##Particle size significance testing ----
#Data prep
PSA_summary_data <- percent_mobilized %>%
  select(condition, size_class, perc_count_mobilized, perc_volume_mobilized) %>%
  left_join(PSA_qtile_speeds %>% select(condition, size_class, Q50_count_v, Q50_vol_v), by = c("condition", "size_class")) %>%
  left_join(metadata %>% select(-c(run, replicate)) %>% distinct(), by = "condition") %>%
  mutate(perc_count_mobilized = perc_count_mobilized/100, #Convert from percent to fraction
         perc_volume_mobilized = perc_volume_mobilized/100)

#Encode factors
PSA_summary_data <- PSA_summary_data %>%
  mutate(
    slope    = factor(slope),
    rainfall = factor(rainfall),
    surface  = factor(surface),
    size_class = factor(size_class)
  )

#ANOVA tests
aov_q50_v_count <- aov(Q50_count_v ~ surface + rainfall + slope + size_class +
                         surface:rainfall + surface:slope + rainfall:slope +
                         surface:size_class + rainfall:size_class + slope:size_class +
                         surface:rainfall:slope + surface:rainfall:size_class +
                         surface:slope:size_class + rainfall:slope:size_class,
                       data = PSA_summary_data)
aov_q50_v_vol <- aov(Q50_vol_v ~ surface + rainfall + slope + size_class +
                       surface:rainfall + surface:slope + rainfall:slope +
                       surface:size_class + rainfall:size_class + slope:size_class +
                       surface:rainfall:slope + surface:rainfall:size_class +
                       surface:slope:size_class + rainfall:slope:size_class,
                     data = PSA_summary_data)
aov_mobilized_count <- aov(perc_count_mobilized ~ surface + rainfall + slope + size_class +
                             surface:rainfall + surface:slope + rainfall:slope +
                             surface:size_class + rainfall:size_class + slope:size_class +
                             surface:rainfall:slope + surface:rainfall:size_class +
                             surface:slope:size_class + rainfall:slope:size_class,
                           data = PSA_summary_data)
aov_mobilized_vol <- aov(perc_volume_mobilized ~ surface + rainfall + slope + size_class +
                           surface:rainfall + surface:slope + rainfall:slope +
                           surface:size_class + rainfall:size_class + slope:size_class +
                           surface:rainfall:slope + surface:rainfall:size_class +
                           surface:slope:size_class + rainfall:slope:size_class,
                         data = PSA_summary_data)

summary(aov_q50_v_count)
summary(aov_q50_v_vol)
summary(aov_mobilized_count)
summary(aov_mobilized_vol)

##Data export ----
write.csv(PSA_qtile_speeds, "data/output_data/PSA_velocity.csv" , row.names = F)
write.csv(percent_mobilized, "data/output_data/PSA_percent_mobilized.csv" , row.names = F)
write.csv(shape_data, "data/output_data/shape_data.csv", row.names = F)
