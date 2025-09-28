#Load libraries ----
library(tidyverse)
library(broom)
library(minpack.lm)

#Load data ----
runoff_data <- read.csv("data/runoff_data.csv")
mass_data <- read.csv("data/mass_data.csv")
metadata <- read.csv("data/metadata.csv")
PSA_data <- read.csv("data/PSA_data.csv")
rain_data <- read.csv("data/rain_data.csv")

#Convert runoff mass to volume -----
# Jones (1992) function: density of air-saturated water
jones_density <- function(temp_C) {
  density_kg_m3 <- 999.84847 +
    6.337563e-2 * temp_C -
    8.523829e-4 * temp_C^2 +
    6.943248e-6 * temp_C^3 -
    3.821216e-8 * temp_C^4
  
  # Convert to g/mL
  density_g_mL <- density_kg_m3 / 1000
  return(density_g_mL)
}

# Calculate density
runoff_data$density_g_mL <- jones_density(runoff_data$water_temp_C)

# Replace missing density rows with mean of available values
mean_density <- mean(runoff_data$density_g_mL, na.rm = TRUE)
runoff_data$density_g_mL[is.na(runoff_data$density_g_mL)] <- mean_density

# Calculate volume (mL)
runoff_data$volume_mL <- runoff_data$water_mass_g / runoff_data$density_g_mL





#Blank calculations ----
#Filter data
blank_mass <- mass_data %>% filter(end_time_min == "Blank")
sample_mass <- mass_data %>% filter(!start_time_min == "N/A") %>%
  mutate(raw_mass_g = jar_and_sample_mass_g - jar_mas_g) #Correction for % calculation because sample_mass_g is already blank subtracted
blank_PSA <- PSA_data %>% filter(interval == "0")

#Mean and SD blank mass, % total sample mass
print(mean(blank_mass$sample_mass_g)) #mean = 0.005259625 g
print(sd(blank_mass$sample_mass_g)) #SD = 0.006978204 g
print(mean(blank_mass$sample_mass_g)/mean(sample_mass$raw_mass_g)) #blank % of avg sample mass = 4.06%
print(median(blank_PSA$FLength)) #D50 = 113.229
print(sd(blank_PSA$FLength)) #SD = 85.56953

#Average mobilized particle size
mobilized_PSA <- PSA_data %>% filter(interval < 7, interval > 0)
print(median(mobilized_PSA$FLength)) #D50 = 161.014
print(sd(mobilized_PSA$FLength)) #SD = 177.1187




#Variability between replicates ----
#Mass-based variability between individual samples
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

##KS tests for mass flux distributions ----
#Mass based KS test function
ks_mass <- function(df) {
  # Get all pairwise combinations of replicates
  replicates <- unique(df$replicate)
  combs <- combn(replicates, 2, simplify = FALSE)
  
  # Run ks.test for each replicate pair
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
  )

print(average_rain) #Average rain gauge data

#Runoff calculated water flux
water_flux <- runoff_data %>%
  filter(end_time_min > 5) %>%
  group_by(run) %>%
  summarize(mean_mL_min = (mean(volume_mL))) %>%
  left_join(metadata %>% select(run, slope, rainfall), by = "run") %>%
  mutate(mean_mL_min = ifelse(rainfall == "high", mean_mL_min * 2, mean_mL_min)) %>% #Correcting for high rainfall runoff collections being every 0.5 mins
  mutate(slope_rad = slope * (pi / 180)) %>%
  mutate(projected_area = (7200)*(cos(slope_rad))) %>% #Choi (2016) equation for projected (rain) area, 120 x 60 cm surface = 7200 cm2
  mutate(avg_mm_min = (mean_mL_min/projected_area)*10)

average_water_flux <- water_flux %>%
  group_by(rainfall) %>%
  summarize(
    n = n(),
    mean_mm_min = mean(avg_mm_min, na.rm = TRUE),
    sd   = sd(avg_mm_min, na.rm = TRUE)
  )

print(average_water_flux) #Average rain data from water flux


#Estimating flow depth ----
#Calculate Re to determine laminar vs turbulent flow
v = 1.0*10^(-6) #Kinematic viscosity (1.0×10^−6 m2/s)
water_flux <- water_flux %>% 
  mutate(q = ((mean_mL_min * 10^(-6))/60)/0.6) %>% #Unit discharge: convert mL/min -> m3/s, divide by flow width (0.6m)
  mutate(Re = q/v) #Re = unit discharge (m2/s) over kinematic viscosity (1.0×10^−6 m2/s)

print(min(water_flux$Re)) #Min Re
print(max(water_flux$Re)) #Max Re

#Calculate flow depth with laminar overland flow equation
g = 9.81 #Acceleration due to gravity (9.81 m/s2)
water_flux <- water_flux %>%
  mutate(S = tan(slope_rad)) %>% #Dimensionless slope
  mutate(y = ((3*v*q)/(g*S))^(1/3)) %>% #Flow depth in m
  mutate(depth_mm = y*1000, #Flow depth converted to mm
         V = q/y) #Flow velocity in m/s
  


#Fitting wash-off models ----
##Data cleaning ----
#Calculate cumulative sample mass
sample_mass_cumulative <- sample_mass %>%
  group_by(run) %>%
  mutate(start_time_min = as.numeric (start_time_min),
         end_time_min = as.numeric(end_time_min)) %>%
  arrange(end_time_min, .by_group = TRUE) %>%  # make sure time is sorted
  mutate(cumulative_mass_g = cumsum(sample_mass_g)) %>%
  ungroup()

w0 = 1 #Initial load, 1 g
wash_off_mass <- sample_mass_cumulative %>% 
  select(-c(jar_mas_g, jar_and_sample_mass_g, raw_mass_g)) %>%
  left_join(metadata %>% select(run, condition), by = "run") %>% #Add condition id
  left_join(water_flux %>% select(run, avg_mm_min), by = "run") %>% 
  rename(t = end_time_min, #t = time of wash-off measurement (min)
         i = avg_mm_min, #i = rainfall intensity (mm/min)
         wt = cumulative_mass_g) %>% #Wt = total wash-off mass (g) at time t
  mutate(frac = wt/w0) #Wash off fraction here equals Wt because W0 = 1 g, but here for future calculations

#Integrate across replicates
wash_off_avg <- wash_off_mass %>%
  group_by(condition, t) %>%
  summarise(
    frac_mean = mean(frac),
    i_mean = mean(i),
    .groups = "drop"
  )

#Uncoupled 2-factor model (capacity factor and k) (Egodawatta et al, 2007) ----
##Check model for single condition ----
model_check <- nls( #Nonlinear least squares
  frac_mean ~ cf * (1 - exp(-k * i_mean * t)),
  data = wash_off_avg %>% filter(condition == 1),
  start = list(cf = 0.5, k = 0.01)  # starting guesses
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
param_table <- map_dfr(fits, function(x) {
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


#Coupled model (f(k) and kprime) (Muthumasy et al, 2018) ----
#Build functions for iterative solving of c* and k'
rss_for_condition <- function(df_cond, cval, kprime) {
  i_val <- unique(df_cond$i_mean)[1]
  fcap  <- cval * kprime
  if (is.na(i_val) || cval <= 0 || kprime <= 0 || fcap > 1) return(Inf)
  pred <- fcap * (1 - exp(-kprime * i_val * df_cond$t))
  sum((df_cond$frac_mean - pred)^2, na.rm = TRUE)
}

fit_kprime_for_c <- function(df_cond, cval) {
  upper <- (1 / cval) - 1e-8
  if (upper <= 0) return(NULL)
  # 1-D optimize over k' in (0, 1/c)
  opt <- optimize(function(kp) rss_for_condition(df_cond, cval, kp),
                  interval = c(1e-8, upper))
  list(kprime = opt$minimum,
       rss    = opt$objective,
       fcap   = cval * opt$minimum)
}

#Solve for c and find min value (c*)
c_grid <- seq(1, 20, by = 0.2)
fit_for_c <- function(cval) {
  fits <- wash_off_avg %>%
    group_by(condition) %>%
    group_map(~{
      f <- fit_kprime_for_c(.x, cval)
      if (is.null(f)) return(NULL)
      tibble(condition = .y$condition, kprime = f$kprime, fcap = f$fcap, rss = f$rss)
    }) %>% bind_rows()
  tibble(c = cval,
         total_rss = sum(fits$rss),
         n_cond = nrow(fits))
}

c_sweep <- map_dfr(c_grid, fit_for_c)
c_star <- c_sweep$c[which.min(c_sweep$total_rss)]
print(c_star) #c* = 4.2 mm

#Fit all conditions with c*
refit_at_cstar <- function(df_cond, cval) {
  f <- fit_kprime_for_c(df_cond, cval)
  if (is.null(f)) return(tibble())
  rss <- f$rss
  n   <- nrow(df_cond)
  p   <- 1L
  rse <- sqrt(rss / pmax(1, n - p))
  
  tibble(
    c_mm      = cval,
    k_prime   = f$kprime,
    f_k       = f$fcap,
    rss       = rss,
    n_points  = n,
    rse       = rse,
    max_frac  = max(df_cond$frac_mean, na.rm = TRUE)
  )
}
#Apply function and extract parameter table
param_table_coupled <- wash_off_avg %>%
  dplyr::group_by(condition) %>%
  dplyr::group_modify(~ refit_at_cstar(.x, c_star)) %>%
  dplyr::ungroup()

#Per-condition intensity & time horizon
cond_meta <- wash_off_avg %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    i_val = unique(i_mean)[1],
    t_max = max(t, na.rm = TRUE),
    .groups = "drop"
  )

#Extract model predictions
predictions_coupled <- param_table_coupled %>%
  dplyr::inner_join(cond_meta, by = "condition") %>%
  purrr::pmap_dfr(function(condition, c_mm, k_prime, f_k, rss, n_points, rse, max_frac, i_val, t_max) {
    t_seq <- seq(0, t_max, length.out = 100)
    tibble::tibble(
      condition = condition,
      t = t_seq,
      pred_frac = f_k * (1 - exp(-k_prime * i_val * t_seq))
    )
  })

#Plot data points vs model predictions for visual confirmation
ggplot() +
  geom_point(data = wash_off_avg, aes(x = t, y = frac_mean), alpha = 0.7, color = "blue") +
  geom_line(data = predictions_coupled, aes(x = t, y = pred_frac), color = "red") +
  facet_wrap(~ condition, scales = "free_x") +
  labs(x = "Time (min)", y = "Wash-off fraction",
       title = "Coupled model (global c*): observed vs fitted by condition") +
  theme_minimal()







