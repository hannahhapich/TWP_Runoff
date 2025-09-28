#Load libraries ----
library(tidyverse)

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
  








