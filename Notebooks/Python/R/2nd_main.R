

# ==============================================================
# BDHS 2022 HYPERTENSION ANALYSIS: COMPLETE JOURNAL PIPELINE
# ==============================================================
# 1. LOAD LIBRARIES
# ==============================================================
library(survey)
library(car)
library(dplyr)
library(broom)
library(ggplot2)
library(forcats)
library(stringr)
library(tidyr)
library(scales)

# ==============================================================
# 2. DATA PREPARATION & SURVEY DESIGN
# ==============================================================
options(survey.lonely.psu = "adjust")

# Load your dataset
df <- read.csv("E:/download/BDHS_2022_HYPERTENSION_SELECTED_VARIABLES.csv")

target <- "hypertension"
categorical_vars <- c(
    "age_group", "sex", "education_level", "marital_status",
    "residence_type", "division", "wealth_index", "BMI_category",
    "diabetes", "household_size_cat", "crowding_category",
    "electricity", "improved_water", "improved_toilet", "mobile_phone"
)

# Convert variables to proper formats
df[[target]] <- as.numeric(df[[target]])
df[categorical_vars] <- lapply(df[categorical_vars], as.factor)

# Create the complex survey design object
bdhs_design <- svydesign(
    ids = ~cluster_id,
    strata = ~strata_id,
    weights = ~sample_weight,
    data = df,
    nest = TRUE
)

# ==============================================================
# 3. SECTIONS 4.1, 4.2 & 4.3: INTEGRATED TABLE 1 (BIVARIATE)
# ==============================================================
# Includes: N (Weighted %), Normal % (CI), Hypertension % (CI), Chi-sq, and P-value

master_bivariate_list <- list()
# ==============================================================
# 3. UPDATED BIVARIATE LOOP: INCLUDES FREQUENCY (WEIGHTED %)
# ==============================================================

for(var in categorical_vars) {
    
    # A. Survey Chi-Square Test (Rao-Scott Correction)
    chi_test <- svychisq(as.formula(paste("~", var, "+", target)), design = bdhs_design)
    chi_sq_val <- round(as.numeric(chi_test$statistic), 2)
    p_val <- chi_test$p.value
    p_label <- if(p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
    
    # B. Overall Weighted Distribution
    w_table <- as.data.frame(svytable(as.formula(paste0("~", var)), design = bdhs_design))
    colnames(w_table)[1] <- var 
    w_table$w_pct <- (w_table$Freq / sum(w_table$Freq)) * 100
    
    # C. Unweighted N
    n_table <- df %>% group_by(!!sym(var)) %>% summarise(n = n(), .groups = 'drop')
    
    # D. Weighted Prevalence and Counts by Hypertension Status
    # We use svyby with 'unwtd.count' and 'svymean' to get both N and %
    prev_stat <- svyby(as.formula(paste0("~", target)), 
                       as.formula(paste0("~", var)), 
                       design = bdhs_design, svymean, na.rm = TRUE)
    
    # Get Weighted Counts for Normal vs Hypertension
    # svytable gives the weighted frequencies for the cross-tab
    counts_weighted <- as.data.frame(svytable(as.formula(paste0("~", var, "+", target)), design = bdhs_design))
    
    # Pivot counts to wide format to merge easily
    counts_wide <- counts_weighted %>%
        mutate(status = ifelse(!!sym(target) == "1", "HTN", "Norm")) %>%
        select(-!!sym(target)) %>%
        pivot_wider(names_from = status, values_from = Freq)
    
    # E. Combine and Format with "Freq (Pct%)"
    var_summary <- n_table %>%
        left_join(w_table, by = var) %>%
        left_join(counts_wide, by = var) %>%
        mutate(
            Variable = var,
            Category = as.character(!!sym(var)),
            # Total Sample column
            Sample_N_Pct = paste0(n, " (", round(w_pct, 1), "%)"),
            # Normal Column: Weighted_N (Weighted_%) and CI
            Normal_Status = paste0(
                round(Norm, 0), " (", 
                round((1 - prev_stat[[target]]) * 100, 1), "%) [",
                round(((1 - prev_stat[[target]]) - (1.96 * prev_stat$se)) * 100, 1), "-",
                round(((1 - prev_stat[[target]]) + (1.96 * prev_stat$se)) * 100, 1), "%]"
            ),
            # Hypertension Column: Weighted_N (Weighted_%) and CI
            Hypertension_Status = paste0(
                round(HTN, 0), " (", 
                round(prev_stat[[target]] * 100, 1), "%) [",
                round((prev_stat[[target]] - (1.96 * prev_stat$se)) * 100, 1), "-",
                round((prev_stat[[target]] + (1.96 * prev_stat$se)) * 100, 1), "%]"
            ),
            Chi_Square = chi_sq_val,
            P_Value = p_label
        ) %>%
        # Standardize Category Labels
        mutate(Category = case_when(
            Variable == "wealth_index" & Category == "1" ~ "Poorest",
            Variable == "wealth_index" & Category == "2" ~ "Poorer",
            Variable == "wealth_index" & Category == "3" ~ "Middle",
            Variable == "wealth_index" & Category == "4" ~ "Richer",
            Variable == "wealth_index" & Category == "5" ~ "Richest",
            Variable == "division" & Category == "1" ~ "Barishal",
            Variable == "division" & Category == "2" ~ "Chattogram",
            Variable == "division" & Category == "3" ~ "Dhaka",
            Variable == "division" & Category == "4" ~ "Khulna",
            Variable == "division" & Category == "5" ~ "Mymensingh",
            Variable == "division" & Category == "6" ~ "Rajshahi",
            Variable == "division" & Category == "7" ~ "Rangpur",
            Variable == "division" & Category == "8" ~ "Sylhet",
            TRUE ~ Category
        )) %>%
        select(Variable, Category, Sample_N_Pct, Normal_Status, Hypertension_Status, Chi_Square, P_Value)
    
    # Clean up repeated stats for sub-rows
    if(nrow(var_summary) > 1) {
        var_summary$Chi_Square[2:nrow(var_summary)] <- NA
        var_summary$P_Value[2:nrow(var_summary)] <- ""
    }
    
    master_bivariate_list[[var]] <- var_summary
}

final_table_1 <- bind_rows(master_bivariate_list)
write.csv(final_table_1, "BDHS_Table1_Bivariate_Integrated.csv", row.names = FALSE, na = "")





# ==============================================================
# FIGURE 1: PREVALENCE BY AGE GROUP AND SEX
# ==============================================================

# 1. Calculate Weighted Prevalence for the Interaction (Age x Sex)
age_sex_prev <- svyby(as.formula(paste0("~", target)), 
                      ~age_group + sex, 
                      design = bdhs_design, 
                      svymean, 
                      na.rm = TRUE)

# 2. Prepare Data for Plotting
plot_data <- as.data.frame(age_sex_prev) %>%
    mutate(
        Prevalence = !!sym(target) * 100,
        Lower = (!!sym(target) - (1.96 * se)) * 100,
        Upper = (!!sym(target) + (1.96 * se)) * 100,
        # Ensure Sex labels are clean for the legend
        sex = case_when(
            sex == "1" ~ "Male",
            sex == "2" ~ "Female",
            TRUE ~ as.character(sex)
        )
    )

# 3. Create the Bar Chart
fig1 <- ggplot(plot_data, aes(x = age_group, y = Prevalence, fill = sex)) +
    # Grouped bars
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    # Error bars (95% CI)
    geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                  position = position_dodge(width = 0.8), 
                  width = 0.25, color = "black") +
    # Formatting
    scale_fill_manual(values = c("Male" = "#1f78b4", "Female" = "#e31a1c")) +
    scale_y_continuous(labels = label_number(suffix = "%"), limits = c(0, 70)) +
    theme_minimal() +
    labs(
        title = "Prevalence of Hypertension by Age Group and Sex",
        x = "Age Group (Years)",
        y = "Weighted Prevalence (%)",
        fill = "Gender",
        caption = "Error bars represent 95% Confidence Intervals (CIs)"
    ) +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.minor = element_blank()
    )

# 4. Save High-Resolution Figure
ggsave("Figure1_Age_Sex_Prevalence.png", fig1, width = 10, height = 6, dpi = 600)

# Display the plot
print(fig1)