# ==============================================================
# BDHS 2022 HYPERTENSION ANALYSIS
# SURVEY-WEIGHTED JOURNAL-READY PIPELINE
# ==============================================================
# INCLUDED:
# 1. Survey Design
# 2. Weighted Prevalence
# 3. Bivariate Analysis (Chi-square + p-value)
# 4. Prevalence Visualization
# 5. Survey-Weighted Logistic Regression
# 6. Adjusted Odds Ratio (AOR) + 95% CI
# 7. Multicollinearity (VIF)
# 8. Publication-Quality Forest Plot
# ==============================================================

# ==============================================================
# 0. INSTALL PACKAGES (RUN ONLY ONCE)
# ==============================================================

install.packages(c(
   "survey",
   "car",
  "dplyr",
   "broom",
   "ggplot2",
   "forcats",
   "stringr",
   "readr",
   "scales"
 ))

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
library(readr)
library(scales)

# ==============================================================
# 2. LOAD DATA
# ==============================================================

df <- read.csv(
    "E:/download/BDHS_2022_HYPERTENSION_SELECTED_VARIABLES.csv"
)

# ==============================================================
# 3. DEFINE VARIABLES
# ==============================================================

target <- "hypertension"

categorical_vars <- c(
    "age_group",
    "sex",
    "education_level",
    "marital_status",
    "residence_type",
    "division",
    "wealth_index",
    "BMI_category",
    "diabetes",
    "household_size_cat",
    "crowding_category",
    "electricity",
    "improved_water",
    "improved_toilet",
    "mobile_phone"
)

# ==============================================================
# 4. DATA PREPARATION
# ==============================================================

df[[target]] <- as.numeric(df[[target]])

for(col in categorical_vars){
    
    df[[col]] <- as.factor(df[[col]])
}

# 1. Ensure the package is loaded in this session first
library(survey) 

# 2. Set options
options(survey.lonely.psu = "adjust")

# 3. Create the survey design object
bdhs_design <- svydesign(
    ids = ~cluster_id,
    strata = ~strata_id,
    weights = ~sample_weight,
    data = df,
    nest = TRUE
)
# ==============================================================
# 6. OVERALL WEIGHTED PREVALENCE
# ==============================================================

overall_prev <- svymean(
    ~hypertension,
    bdhs_design,
    na.rm = TRUE
)

cat("\n=================================================\n")
cat("OVERALL HYPERTENSION PREVALENCE\n")
cat("=================================================\n")

print(overall_prev)

# ==============================================================
# 7. BIVARIATE ANALYSIS
# ==============================================================

# 1. Ensure the required packages are loaded in this session
library(dplyr)
library(survey)

cat("\n=================================================\n")
cat("BIVARIATE ANALYSIS\n")
cat("=================================================\n")

bivariate_results <- data.frame()

for(var in categorical_vars){
    
    # ============================================================
    # SURVEY CHI-SQUARE TEST
    # ============================================================
    
    chi_test <- svychisq(
        as.formula(
            paste("~", var, "+", target)
        ),
        design = bdhs_design
    )
    
    p_value <- chi_test$p.value
    
    # ============================================================
    # WEIGHTED PREVALENCE
    # ============================================================
    
    prev_table <- svyby(
        as.formula(
            paste0("~", target)
        ),
        as.formula(
            paste0("~", var)
        ),
        bdhs_design,
        svymean,
        na.rm = TRUE
    )
    
    # ============================================================
    # CLEAN TABLE
    # ============================================================
    
    prev_table <- prev_table %>%
        mutate(
            Variable = var,
            Category = .data[[var]],
            Prevalence = hypertension * 100,
            P_value = round(p_value, 4)
        ) %>%
        select(
            Variable,
            Category,
            Prevalence,
            P_value
        )
    
    # ============================================================
    # APPEND
    # ============================================================
    
    bivariate_results <- bind_rows(
        bivariate_results,
        prev_table
    )
}

# ==============================================================
# ROUND VALUES & SAVE
# ==============================================================

bivariate_results$Prevalence <- round(bivariate_results$Prevalence, 2)

write.csv(
    bivariate_results,
    "Bivariate_Analysis_Table.csv",
    row.names = FALSE
)

print(bivariate_results)




# ==============================================================
# 7. ENHANCED BIVARIATE ANALYSIS (PUBLICATION STYLE)
# ==============================================================

library(survey)
library(dplyr)
library(tidyr)

# Create an empty list to store results for each variable
all_vars_list <- list()

for(var in categorical_vars) {
    
    # 1. Get Weighted Chi-Square P-value
    chi_test <- svychisq(as.formula(paste("~", var, "+", target)), design = bdhs_design)
    p_val <- chi_test$p.value
    p_label <- if(p_val < 0.001) "p < 0.001" else if(p_val < 0.01) "p < 0.01" else sprintf("%.3f", p_val)
    
    # 2. Get Unweighted Frequencies and Percentages (Total Column)
    freq_total <- df %>%
        group_by(!!sym(var)) %>%
        summarise(n_total = n()) %>%
        mutate(pct_total = (n_total / sum(n_total)) * 100)
    
    # 3. Get Weighted Prevalence by Hypertension Status
    # This calculates the row-wise percentage based on the survey weights
    weighted_prev <- svyby(as.formula(paste0("~", target)), 
                           as.formula(paste0("~", var)), 
                           design = bdhs_design, 
                           svymean, na.rm = TRUE)
    
    # 4. Combine everything into a clean table for this variable
    var_table <- freq_total %>%
        mutate(
            # Hypertension weighted prevalence (from survey design)
            Prev_High = weighted_prev[[target]] * 100,
            Prev_Normal = (1 - weighted_prev[[target]]) * 100,
            # Add metadata
            Variable = var,
            P_Value = p_label
        ) %>%
        rename(Category = 1) %>%
        mutate(
            Frequency_N_Pct = paste0(n_total, " (", round(pct_total, 1), "%)"),
            Normal_n_Pct = paste0(round(Prev_Normal, 1), "%"),
            Hypertension_n_Pct = paste0(round(Prev_High, 1), "%")
        ) %>%
        select(Variable, Category, Frequency_N_Pct, Normal_n_Pct, Hypertension_n_Pct, P_Value)
    
    # Only keep the P-value for the first row of each variable to look like a journal table
    var_table$P_Value[2:nrow(var_table)] <- ""
    
    all_vars_list[[var]] <- var_table
}

# Combine all variables into one final dataframe
final_bivariate_table <- bind_rows(all_vars_list)

# ==============================================================
# VIEW & SAVE
# ==============================================================

# Print to console
print(final_bivariate_table, row.names = FALSE)

# Save to CSV
write.csv(final_bivariate_table, "Journal_Bivariate_Table.csv", row.names = FALSE)


# ==============================================================
# 8. PUBLICATION-READY PREVALENCE VISUALIZATION (FULL CODE)
# ==============================================================

library(survey)
library(dplyr)
library(ggplot2)
library(scales)

# 1. Define variables to visualize
plot_vars <- c("age_group", "BMI_category", "wealth_index", "division")

# 2. Extract Weighted Prevalence with Standard Errors from the Survey Design
plot_list <- list()

for(v in plot_vars) {
    # Calculate mean and SE using the survey design object (bdhs_design)
    stat <- svyby(as.formula(paste0("~", target)), 
                  as.formula(paste0("~", v)), 
                  bdhs_design, svymean, na.rm = TRUE)
    
    # Calculate 95% CI: Mean +/- (1.96 * SE)
    df_stat <- data.frame(
        Variable = v,
        Category = as.character(stat[[1]]),
        Prevalence = stat[[target]] * 100,
        SE = stat$se * 100
    ) %>%
        mutate(
            Lower = Prevalence - (1.96 * SE),
            Upper = Prevalence + (1.96 * SE)
        )
    
    plot_list[[v]] <- df_stat
}

# Combine all variables into one dataframe
final_plot_data <- bind_rows(plot_list)

# 3. Map Numeric Codes to Descriptive Labels (BDHS 2022 Standards)
final_plot_data <- final_plot_data %>%
    mutate(Category = case_when(
        # Wealth Index Mapping
        Variable == "wealth_index" & Category == "1" ~ "Poorest",
        Variable == "wealth_index" & Category == "2" ~ "Poorer",
        Variable == "wealth_index" & Category == "3" ~ "Middle",
        Variable == "wealth_index" & Category == "4" ~ "Richer",
        Variable == "wealth_index" & Category == "5" ~ "Richest",
        
        # Division Mapping
        Variable == "division" & Category == "1" ~ "Barishal",
        Variable == "division" & Category == "2" ~ "Chattogram",
        Variable == "division" & Category == "3" ~ "Dhaka",
        Variable == "division" & Category == "4" ~ "Khulna",
        Variable == "division" & Category == "5" ~ "Mymensingh",
        Variable == "division" & Category == "6" ~ "Rajshahi",
        Variable == "division" & Category == "7" ~ "Rangpur",
        Variable == "division" & Category == "8" ~ "Sylhet",
        
        # Keep other variable categories (Age, BMI) as they are
        TRUE ~ Category
    ))

# 4. Set Factor Levels to maintain logical order (Poorest -> Richest)
# This prevents R from sorting categories alphabetically
final_plot_data$Category <- factor(final_plot_data$Category, levels = c(
    "Poorest", "Poorer", "Middle", "Richer", "Richest",
    "Barishal", "Chattogram", "Dhaka", "Khulna", "Mymensingh", "Rajshahi", "Rangpur", "Sylhet",
    unique(final_plot_data$Category[!final_plot_data$Variable %in% c("wealth_index", "division")])
))

# 5. Create the Final Plot
prevalence_plot <- ggplot(final_plot_data, aes(x = Category, y = Prevalence, fill = Variable)) +
    # Bar layer
    geom_col(alpha = 0.8, color = "black", width = 0.7) +
    # Error bar layer (95% CI)
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "black") +
    # Faceting to separate groups
    facet_wrap(~Variable, scales = "free_x", ncol = 2) +
    # Labels and scales
    labs(
        title = "Weighted Prevalence of Hypertension (BDHS 2022)",
        subtitle = "Adjusted for Survey Weights with 95% Confidence Intervals",
        x = "Socio-demographic Factors",
        y = "Prevalence (%)",
        caption = "Data Source: BDHS 2022 | Error bars represent 95% CI"
    ) +
    scale_y_continuous(
        labels = label_number(suffix = "%"), 
        expand = expansion(mult = c(0, 0.1))
    ) +
    scale_fill_brewer(palette = "Set2") + 
    # High-quality theme for publication
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 16),
        strip.background = element_rect(fill = "grey95", color = "black"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey95")
    )

# 6. Show and Save the Plot
print(prevalence_plot)

ggsave(
    "Hypertension_Prevalence_Final.png",
    prevalence_plot,
    width = 12,
    height = 9,
    dpi = 600
)

cat("\nPlot saved successfully as 'Hypertension_Prevalence_Final.png'\n")

# ==============================================================
# 9. MODEL FORMULA
# ==============================================================

formula_text <- paste(
    target,
    "~",
    paste(categorical_vars, collapse = " + ")
)

model_formula <- as.formula(
    formula_text
)

# ==============================================================
# 10. SURVEY-WEIGHTED LOGISTIC REGRESSION
# ==============================================================

svy_model <- svyglm(
    model_formula,
    design = bdhs_design,
    family = quasibinomial()
)

cat("\n=================================================\n")
cat("SURVEY-WEIGHTED LOGISTIC REGRESSION\n")
cat("=================================================\n")

summary(svy_model)

# ==============================================================
# 11. MULTICOLLINEARITY CHECK (VIF)
# ==============================================================
# ==============================================================
# 11. MULTICOLLINEARITY CHECK (VIF) - FIXED
# ==============================================================

# 1. Load the 'car' library
library(car)

# 2. Run the GLM (unweighted is standard for VIF)
glm_model <- glm(
    model_formula,
    data = df,
    family = binomial()
)

# 3. Calculate VIF
vif_values <- vif(glm_model)

cat("\n=================================================\n")
cat("MULTICOLLINEARITY ASSESSMENT (VIF)\n")
cat("=================================================\n")

print(round(vif_values, 2))
# ==============================================================
# SAVE VIF TABLE
# ==============================================================

vif_table <- as.data.frame(vif_values)

vif_table$Variable <- rownames(
    vif_table
)

rownames(vif_table) <- NULL

vif_table <- vif_table %>%
    select(
        Variable,
        everything()
    )

write.csv(
    vif_table,
    "VIF_Table.csv",
    row.names = FALSE
)

# ==============================================================
# 12-16. REGRESSION RESULTS & CLEANING
# ==============================================================
library(broom)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)

# Extract tidy results with Survey-Adjusted Standard Errors
results <- tidy(svy_model) %>%
    filter(term != "(Intercept)") %>%
    mutate(
        OR = exp(estimate),
        Lower_CI = exp(estimate - 1.96 * std.error),
        Upper_CI = exp(estimate + 1.96 * std.error),
        P_value = p.value
    )

# Feature Label Mapping
feature_labels <- c(
    "age_group" = "Age Group: ",
    "BMI_category" = "BMI Status: ",
    "wealth_index" = "Wealth Quintile: ",
    "division" = "Division: ",
    "sex" = "Gender: ",
    "education_level" = "Education: ",
    "marital_status" = "Marital Status: ",
    "residence_type" = "Residence: ",
    "diabetes" = "Diabetes: ",
    "improved_water" = "Water Source: ",
    "improved_toilet" = "Sanitation: "
)

# Clean Variable Names (e.g., 'age_group35-39' -> 'Age Group: 35-39')
results$Variable <- results$term
for(pattern in names(feature_labels)){
    results$Variable <- str_replace(
        results$Variable, 
        pattern = paste0("^", pattern), 
        replacement = feature_labels[pattern]
    )
}

# Create Forest Plot Display Labels
results <- results %>%
    mutate(
        display_label = paste0(
            Variable, " | AOR: ", 
            sprintf("%.2f", OR), " (", 
            sprintf("%.2f", Lower_CI), "-", 
            sprintf("%.2f", Upper_CI), ")"
        )
    )

# Save Table
write.csv(results %>% select(Variable, OR, Lower_CI, Upper_CI, P_value), 
          "Survey_Weighted_Logistic_Regression_Results.csv", row.names = FALSE)

# ==============================================================
# 17. UPDATED PUBLICATION-READY FOREST PLOT
# ==============================================================

forest_plot <- ggplot(results, aes(x = OR, y = fct_reorder(display_label, OR))) +
    # Reference line at OR = 1
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    # Confidence Intervals
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.3, linewidth = 0.7, color = "midnightblue") +
    # Odds Ratio points
    geom_point(size = 3.5, color = "midnightblue") +
    # Formatting
    labs(
        title = "Factors Associated with Hypertension (BDHS 2022)",
        subtitle = "Adjusted Odds Ratios (AOR) with 95% Confidence Intervals",
        x = "Adjusted Odds Ratio (Log Scale)",
        y = NULL,
        caption = "Note: AOR > 1 indicates increased odds of hypertension."
    ) +
    # Log scale is statistically more accurate for Odds Ratios
    scale_x_log10(breaks = c(0.5, 1, 1.5, 2, 3, 5)) +
    theme_classic(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey30"),
        axis.text.y = element_text(size = 9, family = "mono"), # Monospace keeps the AOR numbers aligned
        panel.grid.major.x = element_line(color = "grey90", linetype = "dotted")
    )

# Print and Save
print(forest_plot)

ggsave(
    "Forest_Plot_AOR_Publication.png",
    forest_plot,
    width = 12,
    height = 10,
    dpi = 600
)












# ==============================================================
# 1. CLEANING THE REGRESSION RESULTS
# ==============================================================
library(broom)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(survey)

# Extract tidy results
results <- tidy(svy_model) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR = exp(estimate),
    Lower_CI = exp(estimate - 1.96 * std.error),
    Upper_CI = exp(estimate + 1.96 * std.error)
  )

# 2. DEFINING THE FEATURE MAPPING (Exact as per your list)
feature_map <- c(
  "age_group" = "Age Group",
  "sex" = "Sex",
  "education_level" = "Education",
  "marital_status" = "Marital Status",
  "residence_type" = "Residence",
  "division" = "Division",
  "wealth_index" = "Wealth",
  "BMI_category" = "BMI",
  "diabetes" = "Diabetes",
  "household_size_cat" = "HH Size",
  "crowding_category" = "Crowding",
  "electricity" = "Electricity",
  "improved_water" = "Water",
  "improved_toilet" = "Toilet",
  "mobile_phone" = "Mobile"
)

# 3. ADVANCED LABEL CLEANING
# This loop identifies the variable, removes it from the term, and creates "Variable: Category"
results <- results %>%
  mutate(Clean_Label = term)

for(var_name in names(feature_map)) {
  results <- results %>%
    mutate(Clean_Label = ifelse(
      str_detect(Clean_Label, paste0("^", var_name)),
      paste0(feature_map[var_name], ": ", str_remove(Clean_Label, paste0("^", var_name))),
      Clean_Label
    ))
}

# 4. FORMATTING THE STATISTICAL TEXT (BOLD AOR & CI)
# Format: **1.02 (0.83–0.90)**
results <- results %>%
  mutate(
    display_stat = paste0(
      sprintf("%.2f", OR), " (", 
      sprintf("%.2f", Lower_CI), "–", 
      sprintf("%.2f", Upper_CI), ")"
    )
  )

# ==============================================================
# 2. FINAL PUBLICATION-READY FOREST PLOT
# ==============================================================

forest_plot <- ggplot(results, aes(x = OR, y = fct_reorder(Clean_Label, OR))) +
  # 1. Reference line at AOR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.6) +
  
  # 2. Confidence Interval bars
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), 
                 height = 0.2, linewidth = 0.8, color = "#2c3e50") +
  
  # 3. Odds Ratio points (Diamonds are preferred in journals)
  geom_point(size = 3.5, shape = 18, color = "#2c3e50") +
  
  # 4. BOLD Statistical Text (placed to the right of the points)
  geom_text(aes(label = display_stat), 
            hjust = -0.15, vjust = -0.6, 
            size = 3.5, fontface = "bold", color = "black") +
  
  # 5. Scaling and Theme
  scale_x_log10(breaks = c(0.5, 1, 2, 3, 5), 
                limits = c(0.4, 8), # Adjusted limits to fit text
                expand = expansion(mult = c(0.05, 0.15))) +
  theme_classic(base_size = 12) +
  labs(
    title = "Factors Associated with Hypertension: BDHS 2022",
    subtitle = "Adjusted Odds Ratios (AOR) with 95% Confidence Intervals",
    x = "Adjusted Odds Ratio (95% CI) - Log Scale",
    y = NULL,
    caption = "Note: Estimates adjusted for all variables in the plot. Ref. categories = baseline (OR=1.00)."
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    panel.grid.major.x = element_line(color = "grey95"),
    plot.margin = margin(10, 20, 10, 10)
  )

# ==============================================================
# 3. SAVE FOR MANUSCRIPT
# ==============================================================
print(forest_plot)

ggsave(
  "Final_Forest_Plot_BDHS2022.png",
  forest_plot,
  width = 11,
  height = 10,
  dpi = 600 # High resolution for printing
)








# ==============================================================
# 1. SETUP & THE PRECISE LOOKUP TABLE
# ==============================================================
library(broom)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(ggtext) # Critical for bold headers

# Define the exact mapping of your categories
label_lookup <- c(
  "age_group18–29" = "  18–29 (Ref)", "age_group30–39" = "  30–39",
  "age_group40–49" = "  40–49", "age_group50–59" = "  50–59",
  "age_group60–69" = "  60–69", "age_group≥70" = "  ≥70",
  "sex1" = "  Male (Ref)", "sex2" = "  Female",
  "education_level1" = "  No education (Ref)", "education_level2" = "  Primary",
  "education_level3" = "  Secondary", "education_level4" = "  Higher",
  "marital_status1" = "  Never married (Ref)", "marital_status2" = "  Married",
  "marital_status3" = "  Widowed/Divorced", "marital_status4" = "  Separated/Other",
  "residence_type1" = "  Urban (Ref)", "residence_type2" = "  Rural",
  "division1" = "  Barishal (Ref)", "division2" = "  Chattogram",
  "division3" = "  Dhaka", "division4" = "  Khulna",
  "division5" = "  Mymensingh", "division6" = "  Rajshahi",
  "division7" = "  Rangpur", "division8" = "  Sylhet",
  "wealth_index1" = "  Poorest (Ref)", "wealth_index2" = "  Poorer",
  "wealth_index3" = "  Middle", "wealth_index4" = "  Richer",
  "wealth_index5" = "  Richest",
  "BMI_category2" = "  Normal (Ref)", "BMI_category1" = "  Underweight",
  "BMI_category3" = "  Overweight", "BMI_category4" = "  Obese",
  "diabetes0" = "  No (Ref)", "diabetes1" = "  Yes",
  "household_size_cat1" = "  Large (Ref)", "household_size_cat2" = "  Medium", "household_size_cat3" = "  Small",
  "crowding_category1" = "  High (Ref)", "crowding_category2" = "  Moderate", "crowding_category3" = "  Low",
  "electricity0" = "  No (Ref)", "electricity1" = "  Yes",
  "improved_water0" = "  No (Ref)", "improved_water1" = "  Yes",
  "improved_toilet0" = "  No (Ref)", "improved_toilet1" = "  Yes",
  "mobile_phone0" = "  No (Ref)", "mobile_phone1" = "  Yes"
)

header_groups <- c(
  "Age Group" = "age_group", "Sex" = "sex", "Education Level" = "education_level",
  "Marital Status" = "marital_status", "Residence Type" = "residence_type",
  "Administrative Division" = "division", "Wealth Quintile" = "wealth_index",
  "BMI Category" = "BMI_category", "Diabetes Status" = "diabetes",
  "Household Size" = "household_size_cat", "Crowding Category" = "crowding_category",
  "Electricity Access" = "electricity", "Improved Water Source" = "improved_water",
  "Improved Sanitation Facility" = "improved_toilet", "Mobile Phone Ownership" = "mobile_phone"
)

# ==============================================================
# 2. DATA PREPARATION (Injecting References and Headers)
# ==============================================================

# Extract results
results_raw <- tidy(svy_model) %>%
  filter(term != "(Intercept)") %>%
  mutate(OR = exp(estimate), Lower = exp(estimate - 1.96 * std.error), Upper = exp(estimate + 1.96 * std.error))

# Create Reference rows
ref_terms <- names(label_lookup)[grepl("Ref", label_lookup)]
ref_data <- data.frame(term = ref_terms, OR = 1.00, Lower = 1.00, Upper = 1.00)

plot_data <- bind_rows(results_raw, ref_data) %>%
  mutate(
    Clean_Label = ifelse(term %in% names(label_lookup), label_lookup[term], term),
    Variable_Group = case_when(
      str_detect(term, "^age_group") ~ "Age Group", str_detect(term, "^sex") ~ "Sex",
      str_detect(term, "^education_level") ~ "Education Level", str_detect(term, "^marital_status") ~ "Marital Status",
      str_detect(term, "^residence_type") ~ "Residence Type", str_detect(term, "^division") ~ "Administrative Division",
      str_detect(term, "^wealth_index") ~ "Wealth Quintile", str_detect(term, "^BMI_category") ~ "BMI Category",
      str_detect(term, "^diabetes") ~ "Diabetes Status", str_detect(term, "^household_size_cat") ~ "Household Size",
      str_detect(term, "^crowding_category") ~ "Crowding Category", str_detect(term, "^electricity") ~ "Electricity Access",
      str_detect(term, "^improved_water") ~ "Improved Water Source", str_detect(term, "^improved_toilet") ~ "Improved Sanitation Facility",
      str_detect(term, "^mobile_phone") ~ "Mobile Phone Ownership", TRUE ~ "Other"
    ),
    bold_stat = ifelse(OR == 1.00, "1.00 (Ref)",
                       paste0(sprintf("%.2f", OR), " (", sprintf("%.2f", Lower), "–", sprintf("%.2f", Upper), ")"))
  )

# Add Headers to data
headers_df <- data.frame(
  Clean_Label = names(header_groups), Variable_Group = names(header_groups),
  OR = NA, Lower = NA, Upper = NA, bold_stat = "AOR (95% CI)", is_header = TRUE
)

plot_data <- plot_data %>% mutate(is_header = FALSE) %>% bind_rows(headers_df)
plot_data$Variable_Group <- factor(plot_data$Variable_Group, levels = names(header_groups))

# THE KEY: Create a wide label combining Category and Stats
plot_data <- plot_data %>%
  arrange(Variable_Group, desc(is_header), term) %>%
  mutate(
    sort_order = row_number(),
    # We use HTML/Markdown to create a "Table" look on the Y-axis
    # <span style='display:inline-block; width:200px'> mimics a table cell
    final_y_label = ifelse(is_header, 
                           paste0("**", Clean_Label, "**"), 
                           paste0(Clean_Label, "   |   *", bold_stat, "*"))
  )

# ==============================================================
# 3. FINAL PLOT
# ==============================================================

forest_plot <- ggplot(plot_data, aes(x = OR, y = fct_reorder(final_y_label, -sort_order))) +
  # Reference line
  geom_vline(xintercept = 1, linetype = "solid", color = "black", linewidth = 0.6) +
  
  # Error bars and Points
  geom_errorbarh(data = filter(plot_data, !is_header),
                 aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = 0.8, color = "#1f78b4") +
  geom_point(data = filter(plot_data, !is_header),
             size = 3.5, shape = 15, color = "#e31a1c") +
  
  # Scaling
  scale_x_log10(position = "top", breaks = c(0.1, 0.5, 1, 2, 5, 10), limits = c(0.1, 15)) +
  
  theme_classic(base_size = 12) +
  labs(x = "Adjusted Odds Ratio (95% CI)", y = NULL) +
  theme(
    # Render labels as markdown for the bold headers and column look
    axis.text.y = element_markdown(color = "black", size = 10, hjust = 0),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.grid.major.x = element_line(color = "grey95")
  )

# ==============================================================
# 4. SAVE (MUST USE TALL HEIGHT)
# ==============================================================
ggsave("Exact_Forest_Table_BDHS.png", forest_plot, width = 12, height = 18, dpi = 600)




































# ==============================================================
# 1. SETUP & LIBRARIES
# ==============================================================
library(broom)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(patchwork) # Essential for side-by-side alignment

# ==============================================================
# 2. DATA PREPARATION
# ==============================================================

# Note: This assumes results_raw is already generated from your svy_model.
# If you are using a standard dataframe, ensure it has: term, OR, Lower, Upper.

# Create Reference rows based on your lookup
ref_terms <- names(label_lookup)[grepl("Ref", label_lookup)]
ref_data <- data.frame(
  term = ref_terms, 
  OR = 1.00, 
  Lower = 1.00, 
  Upper = 1.00,
  estimate = 0,
  std.error = 0
)

plot_data <- bind_rows(results_raw, ref_data) %>%
  mutate(
    Clean_Label = ifelse(term %in% names(label_lookup), label_lookup[term], term),
    Variable_Group = case_when(
      str_detect(term, "^age_group") ~ "Age Group", 
      str_detect(term, "^sex") ~ "Sex",
      str_detect(term, "^education_level") ~ "Education Level", 
      str_detect(term, "^marital_status") ~ "Marital Status",
      str_detect(term, "^residence_type") ~ "Residence Type", 
      str_detect(term, "^division") ~ "Administrative Division",
      str_detect(term, "^wealth_index") ~ "Wealth Quintile", 
      str_detect(term, "^BMI_category") ~ "BMI Category",
      str_detect(term, "^diabetes") ~ "Diabetes Status", 
      str_detect(term, "^household_size_cat") ~ "Household Size",
      str_detect(term, "^crowding_category") ~ "Crowding Category", 
      str_detect(term, "^electricity") ~ "Electricity Access",
      str_detect(term, "^improved_water") ~ "Improved Water Source", 
      str_detect(term, "^improved_toilet") ~ "Improved Sanitation Facility",
      str_detect(term, "^mobile_phone") ~ "Mobile Phone Ownership", 
      TRUE ~ "Other"
    ),
    # Formatting the text column
    bold_stat = ifelse(OR == 1.00, "1.00 (Ref)",
                       paste0(sprintf("%.2f", OR), " (", 
                              sprintf("%.2f", Lower), "–", 
                              sprintf("%.2f", Upper), ")"))
  )

# Add Group Headers
headers_df <- data.frame(
  Clean_Label = names(header_groups), 
  Variable_Group = names(header_groups),
  OR = NA, Lower = NA, Upper = NA, 
  bold_stat = "OR (95% CI)", 
  is_header = TRUE
)

plot_data <- plot_data %>% 
  mutate(is_header = FALSE) %>% 
  bind_rows(headers_df) %>%
  mutate(Variable_Group = factor(Variable_Group, levels = names(header_groups))) %>%
  arrange(Variable_Group, desc(is_header), term) %>%
  mutate(sort_order = row_number())

# Clean labels for display: Headers stay left, sub-items get indented
plot_data <- plot_data %>%
  mutate(
    display_label = ifelse(is_header, Clean_Label, paste0("   ", Clean_Label))
  )

# ==============================================================
# 3. COMPONENT 1: THE TEXT PANEL (Variable Names & Values)
# ==============================================================
# We create a plot that contains only text to act as our table columns
table_panel <- ggplot(plot_data, aes(y = fct_reorder(display_label, -sort_order))) +
  # Variable Column
  geom_text(aes(x = 0, label = display_label, 
                fontface = ifelse(is_header, "bold", "plain")), 
            hjust = 0, size = 3.5) +
  # Statistics Column
  geom_text(aes(x = 2.5, label = bold_stat, 
                fontface = ifelse(is_header, "bold", "plain")), 
            hjust = 0, size = 3.5) +
  theme_void() +
  coord_cartesian(xlim = c(0, 5)) +
  labs(title = "Variable") +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.05),
    plot.margin = margin(5, 0, 5, 5)
  )

# ==============================================================
# 4. COMPONENT 2: THE FOREST PLOT PANEL
# ==============================================================
plot_panel <- ggplot(plot_data, aes(x = OR, y = fct_reorder(display_label, -sort_order))) +
  # Reference line at 1 (for Odds Ratios)
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  # Confidence Interval Bars
  geom_errorbarh(data = filter(plot_data, !is_header),
                 aes(xmin = Lower, xmax = Upper), 
                 height = 0.2, linewidth = 0.6, color = "#4682B4") +
  # Point Estimates
  geom_point(data = filter(plot_data, !is_header),
             size = 2.5, shape = 15, color = "#CD5C5C") +
  # X-Axis Styling
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
  labs(x = "Adjusted Odds Ratio (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 10),
    plot.margin = margin(5, 5, 5, 0)
  )

# ==============================================================
# 5. MERGE AND SAVE
# ==============================================================
# Combine the text panel and plot panel using patchwork
# widths = c(2, 1) means the text takes twice the space of the plot
final_plot <- table_panel + plot_panel + plot_layout(widths = c(1.8, 1))

# Save with enough height to prevent text overlapping
ggsave("Forest_Plot_Final.png", final_plot, width = 14, height = 18, dpi = 300, bg = "white")























# ==============================================================
# 1. SETUP & LIBRARIES
# ==============================================================
library(broom)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(patchwork) # Essential for side-by-side alignment

# ==============================================================
# 2. DATA PREPARATION
# ==============================================================

# Note: This assumes results_raw is already generated from your svy_model.
# If you are using a standard dataframe, ensure it has: term, OR, Lower, Upper.

# Create Reference rows based on your lookup
ref_terms <- names(label_lookup)[grepl("Ref", label_lookup)]
ref_data <- data.frame(
  term = ref_terms, 
  OR = 1.00, 
  Lower = 1.00, 
  Upper = 1.00,
  estimate = 0,
  std.error = 0
)

plot_data <- bind_rows(results_raw, ref_data) %>%
  mutate(
    Clean_Label = ifelse(term %in% names(label_lookup), label_lookup[term], term),
    Variable_Group = case_when(
      str_detect(term, "^age_group") ~ "Age Group", 
      str_detect(term, "^sex") ~ "Sex",
      str_detect(term, "^education_level") ~ "Education Level", 
      str_detect(term, "^marital_status") ~ "Marital Status",
      str_detect(term, "^residence_type") ~ "Residence Type", 
      str_detect(term, "^division") ~ "Administrative Division",
      str_detect(term, "^wealth_index") ~ "Wealth Quintile", 
      str_detect(term, "^BMI_category") ~ "BMI Category",
      str_detect(term, "^diabetes") ~ "Diabetes Status", 
      str_detect(term, "^household_size_cat") ~ "Household Size",
      str_detect(term, "^crowding_category") ~ "Crowding Category", 
      str_detect(term, "^electricity") ~ "Electricity Access",
      str_detect(term, "^improved_water") ~ "Improved Water Source", 
      str_detect(term, "^improved_toilet") ~ "Improved Sanitation Facility",
      str_detect(term, "^mobile_phone") ~ "Mobile Phone Ownership", 
      TRUE ~ "Other"
    ),
    # Formatting the text column
    bold_stat = ifelse(OR == 1.00, "1.00 (Ref)",
                       paste0(sprintf("%.2f", OR), " (", 
                              sprintf("%.2f", Lower), "–", 
                              sprintf("%.2f", Upper), ")"))
  )

# Add Group Headers
headers_df <- data.frame(
  Clean_Label = names(header_groups), 
  Variable_Group = names(header_groups),
  OR = NA, Lower = NA, Upper = NA, 
  bold_stat = "OR (95% CI)", 
  is_header = TRUE
)

plot_data <- plot_data %>% 
  mutate(is_header = FALSE) %>% 
  bind_rows(headers_df) %>%
  mutate(Variable_Group = factor(Variable_Group, levels = names(header_groups))) %>%
  arrange(Variable_Group, desc(is_header), term) %>%
  mutate(sort_order = row_number())

# Clean labels for display: Headers stay left, sub-items get indented
plot_data <- plot_data %>%
  mutate(
    display_label = ifelse(is_header, Clean_Label, paste0("   ", Clean_Label))
  )

# ==============================================================
# 3. COMPONENT 1: THE TEXT PANEL (Variable Names & Values)
# ==============================================================
# We create a plot that contains only text to act as our table columns
table_panel <- ggplot(plot_data, aes(y = fct_reorder(display_label, -sort_order))) +
  # Variable Column
  geom_text(aes(x = 0, label = display_label, 
                fontface = ifelse(is_header, "bold", "plain")), 
            hjust = 0, size = 3.5) +
  # Statistics Column
  geom_text(aes(x = 2.5, label = bold_stat, 
                fontface = ifelse(is_header, "bold", "plain")), 
            hjust = 0, size = 3.5) +
  theme_void() +
  coord_cartesian(xlim = c(0, 5)) +
  labs(title = "Variable") +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.05),
    plot.margin = margin(5, 0, 5, 5)
  )

# ==============================================================
# 4. COMPONENT 2: THE FOREST PLOT PANEL
# ==============================================================
plot_panel <- ggplot(plot_data, aes(x = OR, y = fct_reorder(display_label, -sort_order))) +
  # Reference line at 1 (for Odds Ratios)
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  # Confidence Interval Bars
  geom_errorbarh(data = filter(plot_data, !is_header),
                 aes(xmin = Lower, xmax = Upper), 
                 height = 0.2, linewidth = 0.6, color = "#4682B4") +
  # Point Estimates
  geom_point(data = filter(plot_data, !is_header),
             size = 2.5, shape = 15, color = "#CD5C5C") +
  # X-Axis Styling
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
  labs(x = "Adjusted Odds Ratio (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 10),
    plot.margin = margin(5, 5, 5, 0)
  )

# ==============================================================
# 5. MERGE AND SAVE
# ==============================================================
# Combine the text panel and plot panel using patchwork
# widths = c(2, 1) means the text takes twice the space of the plot
final_plot <- table_panel + plot_panel + plot_layout(widths = c(1.8, 1))

# Save with enough height to prevent text overlapping
ggsave("Forest_Plot_Final.png", final_plot, width = 14, height = 18, dpi = 300, bg = "white")




















# ==============================================================
# 1. SETUP & LIBRARIES
# ==============================================================
library(broom)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(patchwork)

# ==============================================================
# 2. DATA PREPARATION & CLEANING
# ==============================================================

# Extract results from your model
results_raw <- tidy(svy_model, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR = exp(estimate),
    Lower = exp(conf.low),
    Upper = exp(conf.high)
  )

# Add Reference rows
ref_terms <- names(label_lookup)[grepl("Ref", label_lookup)]
ref_data <- data.frame(
  term = ref_terms, OR = 1.00, Lower = 1.00, Upper = 1.00
)

# Calculate "Overall" Row (Weighted Mean OR or from model summary)
# Usually "Overall" is a summary of the whole population (OR = 1 by definition or a calculated rate)
overall_row <- data.frame(
  term = "overall_total",
  Clean_Label = "Overall",
  Variable_Group = "Overall",
  OR = exp(coef(svy_model)[1]), # Adjust this based on your specific 'Overall' requirement
  Lower = exp(confint(svy_model)[1,1]),
  Upper = exp(confint(svy_model)[1,2]),
  is_header = FALSE
)

plot_data <- bind_rows(results_raw, ref_data) %>%
  mutate(
    # FIX: Ensure Clean_Label is pulled from your lookup table
    Clean_Label = ifelse(term %in% names(label_lookup), label_lookup[term], term),
    # Strip "(Ref)" for logic but keep it for display if needed
    Clean_Label = str_replace(Clean_Label, " \\(Ref\\)", ""),
    
    Variable_Group = case_when(
      str_detect(term, "age_group") ~ "Age Group",
      str_detect(term, "sex") ~ "Sex",
      str_detect(term, "education_level") ~ "Education Level",
      str_detect(term, "marital_status") ~ "Marital Status",
      str_detect(term, "residence_type") ~ "Residence Type",
      str_detect(term, "division") ~ "Administrative Division",
      str_detect(term, "wealth_index") ~ "Wealth Quintile",
      str_detect(term, "BMI_category") ~ "BMI Category",
      str_detect(term, "diabetes") ~ "Diabetes Status",
      str_detect(term, "household_size") ~ "Household Size",
      str_detect(term, "crowding") ~ "Crowding Category",
      str_detect(term, "electricity") ~ "Electricity Access",
      str_detect(term, "improved_water") ~ "Improved Water Source",
      str_detect(term, "improved_toilet") ~ "Improved Sanitation Facility",
      str_detect(term, "mobile_phone") ~ "Mobile Phone Ownership",
      TRUE ~ "Other"
    )
  )

# Add Group Headers
headers_df <- data.frame(
  Clean_Label = names(header_groups),
  Variable_Group = names(header_groups),
  OR = NA, Lower = NA, Upper = NA,
  is_header = TRUE
)

# Combine everything
plot_data <- plot_data %>%
  mutate(is_header = FALSE) %>%
  bind_rows(headers_df) %>%
  # Append the Overall row at the very end
  bind_rows(overall_row %>% mutate(Variable_Group = "Overall", is_header = FALSE)) %>%
  mutate(Variable_Group = factor(Variable_Group, levels = c(names(header_groups), "Overall"))) %>%
  arrange(Variable_Group, desc(is_header), term) %>%
  mutate(
    sort_order = row_number(),
    # Formatting the stats column
    bold_stat = case_when(
      is_header ~ "OR (95% CI)",
      OR == 1.00 ~ "1.00 (Ref)",
      is.na(OR) ~ "",
      TRUE ~ paste0(sprintf("%.2f", OR), " (", sprintf("%.2f", Lower), "–", sprintf("%.2f", Upper), ")")
    ),
    display_label = ifelse(is_header, Clean_Label, paste0("   ", Clean_Label))
  )

# ==============================================================
# 3. COMPONENT 1: THE TEXT TABLE
# ==============================================================
table_panel <- ggplot(plot_data, aes(y = fct_reorder(display_label, -sort_order))) +
  geom_text(aes(x = 0, label = display_label, 
                fontface = ifelse(is_header | Clean_Label == "Overall", "bold", "plain")), 
            hjust = 0, size = 3.2) +
  geom_text(aes(x = 3, label = bold_stat, 
                fontface = ifelse(is_header | Clean_Label == "Overall", "bold", "plain")), 
            hjust = 0, size = 3.2) +
  theme_void() +
  coord_cartesian(xlim = c(0, 6)) +
  labs(title = "Variable") +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.02))

# ==============================================================
# 4. COMPONENT 2: THE FOREST PLOT
# ==============================================================
plot_panel <- ggplot(plot_data, aes(x = OR, y = fct_reorder(display_label, -sort_order))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_errorbarh(data = filter(plot_data, !is_header),
                 aes(xmin = Lower, xmax = Upper), height = 0.2, color = "#4682B4") +
  geom_point(data = filter(plot_data, !is_header),
             size = 2, shape = 15, color = "#CD5C5C") +
  # Use log scale if OR values vary widely, or continuous for "Rate" style
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20)) +
  labs(x = "Adjusted Odds Ratio (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 10)
  )

# ==============================================================
# 5. FINAL MERGE
# ==============================================================
final_plot <- table_panel + plot_panel + plot_layout(widths = c(1.5, 1))

ggsave("Forest_Plot_Improved.png", final_plot, width = 12, height = 18, dpi = 300, bg = "white")




































# ==============================================================
# 1. SETUP & LIBRARIES
# ==============================================================
library(broom)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(patchwork)

# ==============================================================
# 2. DATA PREPARATION
# ==============================================================

# Extract results with confidence intervals
results_raw <- tidy(svy_model, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR = exp(estimate),
    Lower = exp(conf.low),
    Upper = exp(conf.high)
  )

# Add Reference rows
ref_terms <- names(label_lookup)[grepl("Ref", label_lookup)]
ref_data <- data.frame(term = ref_terms, OR = 1.00, Lower = 1.00, Upper = 1.00)

# Create "Overall" Data (Manually set this to your population rate/mean OR)
# Note: If you want the population average, calculate it from your survey design object
overall_val <- 1.00 # Set this to your desired 'Overall' point
overall_row <- data.frame(
  term = "Overall", Clean_Label = "Overall", Variable_Group = "Overall",
  OR = overall_val, Lower = overall_val, Upper = overall_val, is_header = FALSE
)

plot_data <- bind_rows(results_raw, ref_data) %>%
  mutate(
    # --- IMPROVED CLEANING LOGIC ---
    # This logic matches your technical names to your pretty labels
    Clean_Label = case_when(
      term %in% names(label_lookup) ~ label_lookup[term],
      # Fuzzy match: if "improved_toilet" is in the term, use the toilet label
      str_detect(term, "improved_toilet") ~ str_replace(term, "improved_toilet", ""),
      str_detect(term, "mobile_phone") ~ str_replace(term, "mobile_phone", ""),
      TRUE ~ term
    ),
    # Remove "(Ref)" for clean grouping logic, but we'll add it back to text
    Clean_Label = str_replace(Clean_Label, " \\(Ref\\)", ""),
    
    Variable_Group = case_when(
      str_detect(term, "age_group") ~ "Age Group",
      str_detect(term, "sex") ~ "Sex",
      str_detect(term, "education") ~ "Education Level",
      str_detect(term, "marital") ~ "Marital Status",
      str_detect(term, "residence") ~ "Residence Type",
      str_detect(term, "division") ~ "Administrative Division",
      str_detect(term, "wealth") ~ "Wealth Quintile",
      str_detect(term, "BMI") ~ "BMI Category",
      str_detect(term, "diabetes") ~ "Diabetes Status",
      str_detect(term, "household_size") ~ "Household Size",
      str_detect(term, "crowding") ~ "Crowding Category",
      str_detect(term, "electricity") ~ "Electricity Access",
      str_detect(term, "improved_water") ~ "Improved Water Source",
      str_detect(term, "improved_toilet") ~ "Improved Sanitation Facility",
      str_detect(term, "mobile_phone") ~ "Mobile Phone Ownership",
      TRUE ~ "Other"
    )
  )

# Add Group Headers
headers_df <- data.frame(
  Clean_Label = names(header_groups),
  Variable_Group = names(header_groups),
  OR = NA, Lower = NA, Upper = NA, is_header = TRUE
)

# Combine and Sort
plot_data <- plot_data %>%
  mutate(is_header = FALSE) %>%
  bind_rows(headers_df) %>%
  bind_rows(overall_row) %>%
  mutate(Variable_Group = factor(Variable_Group, levels = c(names(header_groups), "Overall"))) %>%
  arrange(Variable_Group, desc(is_header), term) %>%
  mutate(
    sort_order = row_number(),
    # Formatting text column
    bold_stat = case_when(
      is_header ~ "OR (95% CI)",
      term %in% ref_terms ~ "1.00 (Ref)",
      is.na(OR) ~ "",
      TRUE ~ paste0(sprintf("%.2f", OR), " (", sprintf("%.2f", Lower), "–", sprintf("%.2f", Upper), ")")
    ),
    display_label = ifelse(is_header, Clean_Label, paste0("   ", Clean_Label))
  )

# ==============================================================
# 3. PLOTTING (The "Table" Look)
# ==============================================================

# Text Panel
table_panel <- ggplot(plot_data, aes(y = fct_reorder(display_label, -sort_order))) +
  geom_text(aes(x = 0, label = display_label, 
                fontface = ifelse(is_header | Clean_Label == "Overall", "bold", "plain")), 
            hjust = 0, size = 3.5) +
  geom_text(aes(x = 3.5, label = bold_stat, 
                fontface = ifelse(is_header | Clean_Label == "Overall", "bold", "plain")), 
            hjust = 0, size = 3.5) +
  theme_void() +
  coord_cartesian(xlim = c(0, 7)) +
  labs(title = "Variable") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.02))

# Forest Panel
plot_panel <- ggplot(plot_data, aes(x = OR, y = fct_reorder(display_label, -sort_order))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_errorbarh(data = filter(plot_data, !is_header),
                 aes(xmin = Lower, xmax = Upper), height = 0.2, color = "#4682B4") +
  geom_point(data = filter(plot_data, !is_header), size = 2.5, shape = 15, color = "#CD5C5C") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20)) +
  labs(x = "Adjusted Odds Ratio (95% CI)", y = NULL) +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

# Merge
final_plot <- table_panel + plot_panel + plot_layout(widths = c(1.8, 1))

ggsave("Forest_Plot_Fixed_Final.png", final_plot, width = 14, height = 18, dpi = 300, bg = "white")






























# ==============================================================
# 1. SETUP & LIBRARIES
# ==============================================================
library(dplyr)
library(ggplot2)
library(patchwork)
library(forcats)

# ==============================================================
# 2. COMPLETE DATASET (EXACTLY FROM YOUR TABLE)
# ==============================================================
plot_data_raw <- tribble(
  ~Variable, ~Category, ~AOR, ~Lower, ~Upper,
  "Age Group", "18–29", 1.00, 1.00, 1.00,
  "Age Group", "30–39", 3.23, 2.59, 4.03,
  "Age Group", "40–49", 6.04, 4.86, 7.51,
  "Age Group", "50–59", 10.99, 8.72, 13.84,
  "Age Group", "60–69", 15.81, 12.33, 20.29,
  "Age Group", "≥70", 23.76, 17.86, 31.61,
  
  "Sex", "Male", 1.00, 1.00, 1.00,
  "Sex", "Female", 1.57, 1.39, 1.77,
  
  "Education Level", "No education", 1.00, 1.00, 1.00,
  "Education Level", "Primary", 0.99, 0.85, 1.16,
  "Education Level", "Secondary", 1.03, 0.88, 1.20,
  "Education Level", "Higher", 1.23, 0.99, 1.52,
  
  "Marital Status", "Never married", 1.00, 1.00, 1.00,
  "Marital Status", "Married", 0.76, 0.57, 1.03,
  "Marital Status", "Widowed/Divorced", 1.02, 0.72, 1.45,
  "Marital Status", "Separated/Other", 1.00, 0.61, 1.62,
  
  "Residence Type", "Urban", 1.00, 1.00, 1.00,
  "Residence Type", "Rural", 0.91, 0.79, 1.06,
  
  "Administrative Division", "Barishal", 1.00, 1.00, 1.00,
  "Administrative Division", "Chattogram", 1.09, 0.87, 1.37,
  "Administrative Division", "Dhaka", 1.00, 0.80, 1.26,
  "Administrative Division", "Khulna", 1.10, 0.88, 1.38,
  "Administrative Division", "Mymensingh", 1.02, 0.80, 1.29,
  "Administrative Division", "Rajshahi", 1.34, 1.05, 1.70,
  "Administrative Division", "Rangpur", 1.26, 1.00, 1.58,
  "Administrative Division", "Sylhet", 1.37, 1.07, 1.77,
  
  "Wealth Quintile", "Poorest", 1.00, 1.00, 1.00,
  "Wealth Quintile", "Poorer", 1.18, 0.97, 1.42,
  "Wealth Quintile", "Middle", 1.07, 0.88, 1.30,
  "Wealth Quintile", "Richer", 1.25, 1.03, 1.53,
  "Wealth Quintile", "Richest", 1.37, 1.09, 1.73,
  
  "BMI Category", "Normal", 1.00, 1.00, 1.00,
  "BMI Category", "Underweight", 0.71, 0.60, 0.84,
  "BMI Category", "Overweight", 2.03, 1.79, 2.31,
  "BMI Category", "Obese", 3.39, 2.88, 3.98,
  
  "Diabetes Status", "No", 1.00, 1.00, 1.00,
  "Diabetes Status", "Yes", 1.38, 1.17, 1.63,
  
  "Household Size", "Large", 1.00, 1.00, 1.00,
  "Household Size", "Medium", 1.08, 0.91, 1.29,
  "Household Size", "Small", 1.23, 1.00, 1.52,
  
  "Crowding Category", "High", 1.00, 1.00, 1.00,
  "Crowding Category", "Moderate", 1.22, 0.99, 1.50,
  "Crowding Category", "Low", 1.35, 1.10, 1.65,
  
  "Electricity Access", "No", 1.00, 1.00, 1.00,
  "Electricity Access", "Yes", 0.75, 0.47, 1.20,
  
  "Improved Water Source", "No", 1.00, 1.00, 1.00,
  "Improved Water Source", "Yes", 0.89, 0.59, 1.33,
  
  "Improved Sanitation Facility", "No", 1.00, 1.00, 1.00,
  "Improved Sanitation Facility", "Yes", 0.84, 0.57, 1.26,
  
  "Mobile Phone Ownership", "No", 1.00, 1.00, 1.00,
  "Mobile Phone Ownership", "Yes", 0.81, 0.56, 1.17
)

# ==============================================================
# 3. DATA PROCESSING
# ==============================================================
# 1. Create header rows for each variable
header_rows <- plot_data_raw %>%
  distinct(Variable) %>%
  mutate(Category = Variable, is_header = TRUE, AOR = NA, Lower = NA, Upper = NA)

# 2. Create the "Overall" summary row
overall_row <- data.frame(
  Variable = "Overall", Category = "Overall", 
  AOR = 1.00, Lower = 1.00, Upper = 1.00, is_header = FALSE
)

# 3. Combine everything
plot_data <- bind_rows(plot_data_raw %>% mutate(is_header = FALSE), header_rows) %>%
  bind_rows(overall_row) %>%
  mutate(Variable = factor(Variable, levels = c(unique(plot_data_raw$Variable), "Overall"))) %>%
  arrange(Variable, desc(is_header)) %>%
  # Filter out duplicate rows where Category equals Variable name unless it's the Header
  filter(!(Category == Variable & !is_header & Variable != "Overall")) %>%
  mutate(
    sort_order = row_number(),
    stat_text = case_when(
      is_header ~ "AOR (95% CI)",
      AOR == 1.00 & Variable != "Overall" ~ "1.00 (Ref)",
      Variable == "Overall" ~ "1.00 (Ref)",
      TRUE ~ paste0(sprintf("%.2f", AOR), " (", sprintf("%.2f", Lower), "–", sprintf("%.2f", Upper), ")")
    ),
    display_label = ifelse(is_header, Category, paste0("   ", Category))
  )

# ==============================================================
# 4. COMPONENT 1: THE TEXT TABLE
# ==============================================================
table_panel <- ggplot(plot_data, aes(y = fct_reorder(display_label, -sort_order))) +
  geom_text(aes(x = 0, label = display_label, 
                fontface = ifelse(is_header | Variable == "Overall", "bold", "plain")), 
            hjust = 0, size = 3.2) +
  geom_text(aes(x = 3.8, label = stat_text, 
                fontface = ifelse(is_header | Variable == "Overall", "bold", "plain")), 
            hjust = 0, size = 3.2) +
  theme_void() +
  coord_cartesian(xlim = c(0, 8)) +
  labs(title = "Variable") +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.05))

# ==============================================================
# 5. COMPONENT 2: THE FOREST PLOT
# ==============================================================
plot_panel <- ggplot(plot_data, aes(x = AOR, y = fct_reorder(display_label, -sort_order))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_errorbarh(data = filter(plot_data, !is_header),
                 aes(xmin = Lower, xmax = Upper), height = 0.2, color = "#4682B4") +
  geom_point(data = filter(plot_data, !is_header),
             size = 2, shape = 15, color = "#CD5C5C") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 35)) +
  labs(x = "Adjusted Odds Ratio (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 10)
  )

# ==============================================================
# 6. COMBINE AND SAVE
# ==============================================================
final_figure <- table_panel + plot_panel + plot_layout(widths = c(1.8, 1))

ggsave("Forest_Plot_Fixed_Complete.png", final_figure, width = 14, height = 20, dpi = 300, bg = "white")