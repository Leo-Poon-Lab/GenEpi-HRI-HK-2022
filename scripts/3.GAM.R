# --- 0. Load Necessary Packages ---
# Ensure these are installed: install.packages(c("mgcv", "MASS", "tidyverse", "lubridate", "forecast", "car"))
library(mgcv)      # For GAMs
library(MASS)      # For glm.nb if needed as alternative/check
library(tidyverse) # For data manipulation (dplyr, tidyr)
library(lubridate) # For handling dates easily
library(forecast)  # For acf, pacf plots
library(car)       # For VIF calculation


# --- 1. Load and Prepare Your Specific Data ---

# --- 1a. Load Markov Jump Data (Community to Hospital) ---
data_markov_jump <- read_csv("results/ruixuan/Hospital infection analysis data_2/ba22_mj_timing.csv")
# Filter transitions and normalize by 501 trees
community_to_hospital <- data_markov_jump %>% 
  filter(startLocation == "community", endLocation == "hospital") %>% 
  mutate(transition = "community_to_hospital")
hospital_to_community <- data_markov_jump %>% 
  filter(startLocation == "hospital", endLocation == "community") %>% 
  mutate(transition = "hospital_to_community")

# Combine and normalize
combined <- bind_rows(community_to_hospital, hospital_to_community) %>% 
  mutate(rate = 1/501)  # Each jump contributes 1/501
combined$date <- ggtree::decimal2Date(combined$time)

# Create daily bins for 2022
all_days <- seq(ymd("2022-05-28"), ymd("2022-08-18"), by = "day")

# Create daily rates
daily_rates <- combined %>% 
  mutate(day = floor_date(date, "day")) %>% 
  group_by(day, transition) %>% 
  summarize(daily_rate = sum(rate), .groups = "drop") %>% 
  complete(day = all_days, transition, fill = list(daily_rate = 0)) %>% 
  group_by(transition) %>% 
  ungroup()

daily_rates_community_to_hospital <- daily_rates %>% 
  filter(transition == "community_to_hospital") %>% 
  select(day, daily_rate) %>% 
  rename(Date=day, MJ_value = daily_rate)

# --- 1b. Load Predictor Data ---
data_local_transport_plot <- read_csv("data/mobility_data/parsed/local_transport_long.csv")
data_cross_border_2022 <- read_csv("data/mobility_data/parsed/cross_border_2022.csv")
data_google_long <- read_csv("data/mobility_data/parsed/google_mobility_long.csv")

data_local_transport_plot <- data_local_transport_plot %>%
  pivot_wider(
    names_from = Transport,
    values_from = `No. of passengers`
  ) %>% 
  rename_with(~paste0("local_transport_", .), -Date) %>%
  filter(Date >= "2022-05-28" & Date <= "2022-08-18") 

data_local_transport_plot <- data_local_transport_plot %>% 
  mutate(data_local_transport_total = rowSums(dplyr::select(., starts_with("local_transport_")), na.rm = TRUE)) %>% 
  dplyr::select(Date, data_local_transport_total) # Keep only the total transport data, due to collinearity

data_cross_border_2022 <- data_cross_border_2022 %>%
  mutate(border_direction = paste0(Border, "_", `Arrival / Departure`)) %>% 
  dplyr::select(Date, border_direction, Total) %>% 
  pivot_wider(
    names_from = border_direction,
    values_from = Total
  ) %>% 
  rename_with(~paste0("cross_border_", .), -Date) %>% 
  filter(Date >= "2022-05-28" & Date <= "2022-08-18")
data_cross_border_2022 <- data_cross_border_2022 %>% 
  dplyr::select(-starts_with("cross_border_Harbour")) # Remove Harbour columns since they are mostly zeros

data_google_long <- data_google_long %>% 
  mutate(type_location = paste0(Type, "_", Location)) %>% 
  dplyr::select(Date=date, type_location, Value=Mobility) %>% 
  pivot_wider(names_from = type_location, values_from = Value) %>%
  rename_with(~paste0("google_index_", .), -Date) %>% 
  filter(Date >= "2022-05-28" & Date <= "2022-08-18")

data_predictors <- data_local_transport_plot %>% 
  left_join(data_cross_border_2022, by = "Date") %>% 
  left_join(data_google_long, by = "Date") %>%
  rename_with(~gsub(" ", "_", .))

# GAM
# --- 1c. Combine Response and Predictors ---
# Ensure Date columns are of the same type if necessary (they should be Date objects)
analysis_data <- daily_rates_community_to_hospital %>%
  left_join(data_predictors, by = "Date") %>%
  # Make sure data is ordered by date for time index creation
  arrange(Date) %>%
  filter(Date >= "2022-05-28" & Date <= "2022-08-18") %>% 
  # Create a numeric time index (1, 2, 3...)
  mutate(time_index = row_number()) %>%
  # Optional: Handle potential missing values in predictors if any (e.g., imputation or filtering)
  # analysis_data <- analysis_data %>% na.omit() # Simplest approach: remove rows with any NAs
  drop_na() # Alternative using drop_na()

# Check the structure of the final analysis dataset
glimpse(analysis_data)
summary(analysis_data$MJ_value) # Look at the distribution of your response

# --- 2. Assess Collinearity Among ALL Potential Predictors ---

# Identify all predictor column names in 'data_predictors' (excluding 'Date')
all_predictor_colnames <- setdiff(colnames(data_predictors), "Date")

# Check if all expected predictor columns exist in the joined analysis_data
missing_cols <- setdiff(all_predictor_colnames, colnames(analysis_data))
if (length(missing_cols) > 0) {
  warning("Some predictor columns were not found in the joined data: ", paste(missing_cols, collapse=", "))
  # Remove missing columns from the list to proceed
  all_predictor_colnames <- intersect(all_predictor_colnames, colnames(analysis_data))
}

# Select only the predictor columns from the prepared analysis_data
# Ensure no missing values ONLY within these columns for VIF/cor calculation
predictor_subset_for_diag <- analysis_data %>%
  dplyr::select(all_of(all_predictor_colnames)) %>%
  na.omit() # Use only rows complete for all predictors for diagnostics

# Calculate and print correlation matrix
if(nrow(predictor_subset_for_diag) > 1) {
    cor_matrix <- cor(predictor_subset_for_diag, use = "pairwise.complete.obs") # Use pairwise to handle potential remaining NAs robustly
    print("Correlation Matrix (rounded):")
    print(round(cor_matrix, 2))

    # Calculate VIFs using a temporary linear model
    # The outcome variable in the formula doesn't really matter for VIF calculation
    # Make sure the data used here has no NAs in the predictors
    # Need to ensure MJ_value exists for the formula structure, join it temporarily
    vif_data_temp <- analysis_data %>%
      dplyr::select(MJ_value, all_of(all_predictor_colnames)) %>%
      na.omit() # VIF requires complete cases for all variables in the formula

    if(nrow(vif_data_temp) > length(all_predictor_colnames) + 1) { # Need enough observations
        vif_model <- lm(MJ_value ~ ., data = vif_data_temp)
        vif_values <- tryCatch({
            vif(vif_model)
        }, error = function(e) {
            warning("VIF calculation failed. Check for perfect collinearity or other issues. Error: ", e$message)
            NULL
        })

        if (!is.null(vif_values)) {
            print("Variance Inflation Factors (VIFs):")
            print(vif_values)
        }
    } else {
      warning("Not enough complete observations to calculate VIFs after removing NAs.")
    }

} else {
  warning("Not enough data after removing NAs in predictors to calculate correlations or VIFs.")
}


# --- 3. Select Predictors Based on Collinearity ---

# USER ACTION: Examine the correlation matrix and VIF values carefully.
# Decide which predictors to include in the final model based on:
# - High VIFs (e.g., > 5 or 10 indicate problems)
# - High pairwise correlations (e.g., |r| > 0.7 or 0.8)
# - Theoretical relevance (which variables make most sense?)
# Create a vector 'selected_predictors' containing only the names of the columns you want to keep.

# Example (replace with your actual selection based on diagnostics):
selected_predictors <- c("data_local_transport_total", "cross_border_Airport_Arrival", "cross_border_Airport_Departure","cross_border_Mainland_border_Arrival", "cross_border_Mainland_border_Departure", "google_index_Public_Parks", "google_index_Public_Grocery_and_pharmacy", "google_index_Private_Workplaces")

lm(MJ_value ~ ., data = vif_data_temp %>% dplyr::select(MJ_value, all_of(selected_predictors))) %>% vif() # All VIFs are < 5 now.

# Make sure the selected predictors actually exist in the data
selected_predictors <- intersect(selected_predictors, colnames(analysis_data))
if(length(selected_predictors) == 0) {
    stop("No valid predictors selected or found in the data. Please check the 'selected_predictors' vector.")
}
print("Predictors selected for GAM model:")
print(selected_predictors)

# --- 4. Fit the GAM Model ---

# Define the model formula using selected predictors and the time smooth
# Choose k: Start based on data length (~sqrt(N) or ~10-20), check later
k_time_smooth <- max(10, min(30, floor(nrow(analysis_data) / 4))) # Heuristic starting k
gam_formula <- reformulate(
  termlabels = c(paste0("s(time_index, bs = 'cr', k = ", k_time_smooth, ")"),
                 selected_predictors),
  response = "MJ_value"
)

print("Fitting GAM with formula:")
print(gam_formula)

# NOTE on family: Your MJ_value is a rate (sum of 1/N_trees).
# Negative Binomial (nb) often works well for overdispersed count-like data.
# However, since it's not strictly counts, diagnostics are crucial.
# Alternatives if NB diagnostics look poor: gaussian(), Gamma(), quasipoisson().
gam_model <- gam(gam_formula,
                 data = analysis_data,
                 family = nb(),     # Using Negative Binomial - CHECK DIAGNOSTICS CAREFULLY
                 method = "REML")  # Use REML for better smoothness estimation


# --- 5. Check the Model ---

# Get a summary of the fitted model
summary(gam_model)
# Look for: Significance of predictors, edf for s(time_index).

# Perform diagnostic checks
par(mfrow = c(2, 2)) # Arrange plots in a 2x2 grid
gam.check(gam_model)
par(mfrow = c(1, 1)) # Reset plot layout
# Check: k-index p-value, residual plots, QQ plot. Is k adequate? Does NB seem okay?

# Check for remaining autocorrelation in residuals
# Using response residuals (difference between observed and fitted rates)
model_residuals <- residuals(gam_model, type = "response")
# Check if there are enough non-NA residuals to plot ACF/PACF
if(sum(!is.na(model_residuals)) > 20) { # Need a reasonable number of points
    acf(model_residuals, main = "ACF of Residuals", na.action = na.pass)
    pacf(model_residuals, main = "PACF of Residuals", na.action = na.pass)
} else {
    warning("Not enough non-NA residuals to plot ACF/PACF reliably.")
}
# Significant spikes (esp. at low lags) indicate remaining autocorrelation.


# --- 6. Interpret the Results ---

# Extract coefficients for parametric terms (your predictors)
model_summary <- summary(gam_model)
para_coeffs <- model_summary$p.coeff
para_se <- model_summary$se[names(para_coeffs)] # Match SEs using names
para_p_values <- model_summary$p.pv[names(para_coeffs)]

# Ensure we only select results for the predictors (exclude Intercept if needed)
predictor_indices <- names(para_coeffs) %in% selected_predictors
# Add Intercept back if you want to see it in the table
# predictor_indices <- names(para_coeffs) %in% c("(Intercept)", selected_predictors)


# Calculate IRRs (Incidence Rate Ratios for NB family) and 95% CIs
results_table <- tibble(
  predictor = names(para_coeffs[predictor_indices]),
  estimate_log = para_coeffs[predictor_indices],
  se = para_se[predictor_indices],
  z_value = estimate_log / se,
  p_value = para_p_values[predictor_indices],
  # IRR assumes a log link, which is default for nb()
  irr = exp(estimate_log),
  irr_ci_lower = exp(estimate_log - 1.96 * se),
  irr_ci_upper = exp(estimate_log + 1.96 * se)
)

print("Model Results (Selected Parametric Predictors):")
print(results_table, n=Inf) # Show all rows
results_table_round <- results_table %>% mutate_at(vars(estimate_log, se, z_value, p_value, irr, irr_ci_lower, irr_ci_upper), round, 3)
write_csv(results_table_round, "results/gam_model_results.csv")

# Interpretation Example (for a significant predictor 'X'):
# "A one-unit increase in 'X' is associated with an estimated Y% [increase/decrease]
# (where Y = (IRR - 1) * 100) in the daily rate of community-to-hospital Markov jumps,
# holding other included predictors constant and adjusting for underlying temporal
# trends captured by the smooth term (p = Z)."

# Plot the estimated smooth function of time
plot(gam_model, select = 1, shade = TRUE, residuals = TRUE, pch = 1, cex = 0.8,
  main = "Estimated Smooth Function of Time Index", ylab = "Effect on log(Expected Jumps)")
# Save the smooth function plot to a PDF file
pdf("results/gam_smooth_time_function.pdf", width = 8, height = 6)
plot(gam_model, select = 1, shade = TRUE, residuals = TRUE, pch = 1, cex = 0.8,
  main = "Estimated Smooth Function of Time Index", ylab = "Effect on log(Expected Jumps)")
dev.off()

# Save diagnostic plots to a PDF file
pdf("results/gam_diagnostics.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
gam.check(gam_model)
dev.off()

# Save ACF and PACF plots to a PDF file
pdf("results/gam_acf_pacf.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
acf(model_residuals, main = "ACF of Residuals", na.action = na.pass)
pacf(model_residuals, main = "PACF of Residuals", na.action = na.pass)
dev.off()