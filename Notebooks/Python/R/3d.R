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

for(var in categorical_vars) {
    
    # A. Survey Chi-Square Test (Rao-Scott Correction)
    chi_test <- svychisq(as.formula(paste("~", var, "+", target)), design = bdhs_design)
    chi_sq_val <- round(as.numeric(chi_test$statistic), 2)
    p_val <- chi_test$p.value
    p_label <- if(p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
    
    # B. Weighted Distribution (Sample Characteristics)
    w_table <- as.data.frame(svytable(as.formula(paste0("~", var)), design = bdhs_design))
    colnames(w_table)[1] <- var 
    w_table$w_pct <- (w_table$Freq / sum(w_table$Freq)) * 100
    
    # C. Unweighted N
    n_table <- df %>% group_by(!!sym(var)) %>% summarise(n = n(), .groups = 'drop')
    
    # D. Weighted Prevalence with 95% CI
    prev_stat <- svyby(as.formula(paste0("~", target)), 
                       as.formula(paste0("~", var)), 
                       design = bdhs_design, svymean, na.rm = TRUE)
    
    # Calculate 95% CI for both statuses
    prev_stat <- prev_stat %>%
        mutate(
            HTN_Prev = !!sym(target) * 100,
            HTN_Lower = (!!sym(target) - (1.96 * se)) * 100,
            HTN_Upper = (!!sym(target) + (1.96 * se)) * 100,
            Norm_Prev = (1 - !!sym(target)) * 100,
            Norm_Lower = ((1 - !!sym(target)) - (1.96 * se)) * 100,
            Norm_Upper = ((1 - !!sym(target)) + (1.96 * se)) * 100
        )
    
    # E. Combine, Format, and Map Labels
    var_summary <- n_table %>%
        left_join(w_table, by = var) %>%
        mutate(
            Variable = var,
            Category = as.character(!!sym(var)),
            Sample_N_Pct = paste0(n, " (", round(w_pct, 1), "%)"),
            Normal_95_CI = paste0(round(prev_stat$Norm_Prev, 1), 
                                  "% (", round(prev_stat$Norm_Lower, 1), 
                                  "-", round(prev_stat$Norm_Upper, 1), "%)"),
            Hypertension_95_CI = paste0(round(prev_stat$HTN_Prev, 1), 
                                        "% (", round(prev_stat$HTN_Lower, 1), 
                                        "-", round(prev_stat$HTN_Upper, 1), "%)"),
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
        select(Variable, Category, Sample_N_Pct, Normal_95_CI, Hypertension_95_CI, Chi_Square, P_Value)
    
    # Clean up repeated stats for sub-rows
    if(nrow(var_summary) > 1) {
        var_summary$Chi_Square[2:nrow(var_summary)] <- NA
        var_summary$P_Value[2:nrow(var_summary)] <- ""
    }
    
    master_bivariate_list[[var]] <- var_summary
}

final_table_1 <- bind_rows(master_bivariate_list)
write.csv(final_table_1, "ssssBDHS_Table1_Bivariate_Analysis.csv", row.names = FALSE, na = "")

# ==============================================================
# 4.4 MULTIVARIABLE LOGISTIC REGRESSION
# ==============================================================
# A. Multicollinearity Check
vif_model <- glm(as.formula(paste(target, "~", paste(categorical_vars, collapse = " + "))), 
                 data = df, family = binomial())
vif_results <- as.data.frame(vif(vif_model))
write.csv(vif_results, "Multicollinearity_VIF.csv")

# B. Survey-Weighted Logistic Regression
svy_model <- svyglm(
    as.formula(paste(target, "~", paste(categorical_vars, collapse = " + "))),
    design = bdhs_design,
    family = quasibinomial()
)

# C. Tidy Regression Table (AOR & CI)
regression_table <- tidy(svy_model) %>%
    filter(term != "(Intercept)") %>%
    mutate(
        AOR = exp(estimate),
        Lower_CI = exp(estimate - 1.96 * std.error),
        Upper_CI = exp(estimate + 1.96 * std.error),
        AOR_95_CI = paste0(round(AOR, 2), " (", round(Lower_CI, 2), "-", round(Upper_CI, 2), ")"),
        P_Value = ifelse(p.value < 0.001, "<0.001", round(p.value, 3))
    ) %>%
    select(term, AOR_95_CI, P_Value)

write.csv(regression_table, "BDHS_Table2_Regression_Results.csv", row.names = FALSE)

# ==============================================================
# 5. PUBLICATION FOREST PLOT
# ==============================================================
plot_results <- tidy(svy_model) %>%
    filter(term != "(Intercept)") %>%
    mutate(
        OR = exp(estimate),
        Lower = exp(estimate - 1.96 * std.error),
        Upper = exp(estimate + 1.96 * std.error)
    )

forest_plot <- ggplot(plot_results, aes(x = OR, y = fct_reorder(term, OR))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, color = "midnightblue") +
    geom_point(size = 3, color = "midnightblue") +
    scale_x_log10(breaks = c(0.5, 1, 2, 4, 8)) +
    labs(title = "Factors Associated with Hypertension (BDHS 2022)",
         subtitle = "Adjusted Odds Ratios with 95% Confidence Intervals",
         x = "Adjusted Odds Ratio (Log Scale)", y = "") +
    theme_classic()

ggsave("Forest_Plot_AOR.png", forest_plot, width = 10, height = 8, dpi = 600)

cat("\n--- ALL ANALYSES COMPLETE ---\n")
cat("Files generated: Table 1 (Bivariate), Table 2 (Regression), VIF Analysis, and Forest Plot.\n")











# ==============================================================
# 6. SPATIAL ANALYSIS: CHOROPLETH MAP (FULL BLOCK)
# ==============================================================
# Ensure you have the 'sf' package installed: install.packages("sf")
install.packages("sf")
library(sf)
library(ggplot2)
library(dplyr)
library(scales)
library(survey)

# 1. Calculate Design-Adjusted Prevalence by Division
# This uses your existing 'bdhs_design' object
div_prev <- svyby(as.formula(paste0("~", target)), 
                  ~division, 
                  design = bdhs_design, 
                  svymean, 
                  na.rm = TRUE)

# 2. Prepare Data for Mapping
# We map the BDHS codes (1-8) to the exact names in your GeoJSON file
div_map_data <- as.data.frame(div_prev) %>%
    mutate(
        Prevalence = !!sym(target) * 100,
        # Matching names to the "ADM1_EN" property in your JSON file
        Division_Name = case_when(
            division == "1" ~ "Barisal",
            division == "2" ~ "Chittagong",
            division == "3" ~ "Dhaka",
            division == "4" ~ "Khulna",
            division == "5" ~ "Mymensingh",
            division == "6" ~ "Rajshahi",
            division == "7" ~ "Rangpur",
            division == "8" ~ "Sylhet"
        )
    )

# 3. Load the GeoJSON file (Division Level)
# Ensure the file is in your working directory
bgd_shp <- st_read("small_bangladesh_geojson_adm1_8_divisions_bibhags.json")

# 4. Join the Survey Results with the Geographic Shapes
map_final <- bgd_shp %>%
    left_join(div_map_data, by = c("ADM1_EN" = "Division_Name"))

# 5. Generate the Plot
htn_map <- ggplot(data = map_final) +
    # Draw the map boundaries
    geom_sf(aes(fill = Prevalence), color = "white", size = 0.3) +
    # Apply a high-contrast color scale (Magma is great for health data)
    scale_fill_viridis_c(
        option = "magma", 
        direction = -1, # Darker colors for higher prevalence
        name = "Prevalence (%)",
        labels = label_number(suffix = "%")
    ) +
    # Clean up the background
    theme_void() + 
    labs(
        title = "Weighted Prevalence of Hypertension by Division",
        subtitle = "Bangladesh Demographic and Health Survey (BDHS) 2022",
        caption = "Note: Prevalence calculated using survey weights.\nMap Boundaries: ADM1 Bibhag Level."
    ) +
    theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 10)),
        plot.caption = element_text(size = 8, color = "grey30", hjust = 0)
    )

# 6. Export for Publication
# 600 DPI ensures the lines and text are crisp for your report
ggsave("Figure_1_Hypertension_Map.png", htn_map, width = 7, height = 9, dpi = 600)

# 7. Print Check
print(htn_map)





# ==============================================================
# 6. SPATIAL ANALYSIS: CHOROPLETH MAP WITH DATA LABELS
# ==============================================================
library(sf)
library(ggplot2)
library(dplyr)
library(scales)
library(survey)

# 1. Calculate Design-Adjusted Prevalence by Division
# (Uses your existing survey design object)
div_prev <- svyby(as.formula(paste0("~", target)), 
                  ~division, 
                  design = bdhs_design, 
                  svymean, 
                  na.rm = TRUE)

# ==============================================================
# 2. PREPARE DATA WITH COMBINED LABELS (NAME + VALUE)
# ==============================================================
div_map_data <- as.data.frame(div_prev) %>%
    mutate(
        Prevalence = !!sym(target) * 100,
        Division_Name = case_when(
            division == "1" ~ "Barisal",
            division == "2" ~ "Chittagong",
            division == "3" ~ "Dhaka",
            division == "4" ~ "Khulna",
            division == "5" ~ "Mymensingh",
            division == "6" ~ "Rajshahi",
            division == "7" ~ "Rangpur",
            division == "8" ~ "Sylhet"
        ),
        # Combine Name and Value with a newline (\n) for a professional look
        Label_Full = paste0(Division_Name, "\n", round(Prevalence, 1), "%")
    )

# ... [Step 3 & 4 remain the same] ...

# ==============================================================
# 5. GENERATE LABELED MAP (UPDATED LABEL AES)
# ==============================================================
htn_map_labeled <- ggplot(data = map_final) +
    geom_sf(aes(fill = Prevalence), color = "white", size = 0.4) +
    
    # Updated to use Label_Full
    geom_sf_text(aes(label = Label_Full), 
                 color = "white",      # White often looks cleaner on 'magma'
                 fontface = "bold", 
                 lineheight = 0.9,     # Adjust spacing between name and %
                 size = 3.5) +
    
    scale_fill_viridis_c(
        option = "magma", 
        direction = -1, 
        name = "Prevalence (%)",
        labels = label_number(suffix = "%")
    ) +
    theme_void() + 
    labs(
        title = "Hypertension Prevalence by Division (BDHS 2022)",
        subtitle = "Geographic Distribution with Point Estimates",
        caption = "Data Source: Bangladesh Demographic and Health Survey 2022\nNote: Values are survey-weighted percentages."
    ) +
    theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
        plot.caption = element_text(size = 9, color = "grey30", hjust = 0)
    )

# 6. Save high-resolution version
ggsave("ssFigure_1_Labeled_Map_HTN.png", htn_map_labeled, width = 12, height = 8, dpi = 600)

# Display result
print(htn_map_labeled)









# ==============================================================
# BDHS 2022: SPATIAL PREVALENCE MAP (FULL SCRIPT)
# ==============================================================

# 1. Load Required Libraries
library(sf)
library(ggplot2)
library(dplyr)
library(scales)
library(survey)

# Note: This assumes you have already defined 'bdhs_design' 
# and 'target' (e.g., target <- "hypertension")

# 2. Calculate Design-Adjusted Prevalence by Division
# This produces a summary table with weighted means
div_prev <- svyby(as.formula(paste0("~", target)), 
                  ~division, 
                  design = bdhs_design, 
                  svymean, 
                  na.rm = TRUE)

# 3. Prepare Data and Create Combined Labels
# We map the numeric codes (1-8) to their English names
div_map_data <- as.data.frame(div_prev) %>%
    mutate(
        Prevalence = !!sym(target) * 100,
        Division_Name = case_when(
            division == "1" ~ "Barisal",
            division == "2" ~ "Chittagong",
            division == "3" ~ "Dhaka",
            division == "4" ~ "Khulna",
            division == "5" ~ "Mymensingh",
            division == "6" ~ "Rajshahi",
            division == "7" ~ "Rangpur",
            division == "8" ~ "Sylhet"
        ),
        # Create a multi-line label: "Division Name \n 00.0%"
        Label_Full = paste0(Division_Name, "\n", round(Prevalence, 1), "%")
    )

# 4. Load the GeoJSON Geography
# Ensure the .json file is in your working directory
bgd_shp <- st_read("small_bangladesh_geojson_adm1_8_divisions_bibhags.json")

# 5. Join the Statistics to the Spatial Map
# We join using ADM1_EN from the JSON and Division_Name from our data
map_final <- bgd_shp %>%
    left_join(div_map_data, by = c("ADM1_EN" = "Division_Name"))

# 6. Generate the Final Labeled Map
htn_map_labeled <- ggplot(data = map_final) +
    # Draw polygons with survey-weighted fill
    geom_sf(aes(fill = Prevalence), color = "white", size = 0.4) +
    
    # Add combined Name and % labels at the center of each division
    geom_sf_text(aes(label = Label_Full), 
                 color = "green",      # White text for dark 'magma' background
                 fontface = "bold", 
                 lineheight = 0.85,    # Tighten spacing between lines
                 size = 3.2) +         # Adjust size based on your output needs
    
    # Styling with the Magma color palette
    scale_fill_viridis_c(
        option = "magma", 
        direction = -1, 
        name = "Prevalence (%)",
        labels = label_number(suffix = "%")
    ) +
    
    # Clean up the map layout
    theme_void() + 
    labs(
        title = "Hypertension Prevalence by Division (BDHS 2022)",
        subtitle = "Geographic Distribution of Survey-Weighted Estimates",
        caption = "Source: Bangladesh Demographic and Health Survey (BDHS) 2022\nNote: Map displays administrative level 1 divisions."
    ) +
    theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 10), color = "grey20"),
        plot.caption = element_text(size = 9, color = "grey40", hjust = 0, margin = margin(t = 10))
    )

# 7. Save the Map in High Resolution
ggsave("Figure_1_BDHS_HTN_Map.png", 
       plot = htn_map_labeled, 
       width = 12, 
       height = 8, 
       dpi = 600, 
       bg = "white")

# 8. Display result
print(htn_map_labeled)
















# ==============================================================
# SPATIAL ANALYSIS: FACETED MAP (MALE vs FEMALE)
# ==============================================================
library(sf)
library(survey)
library(ggplot2)
library(dplyr)
library(scales)

# 1. Calculate Weighted Prevalence Stratified by Division AND Sex
# This creates a more detailed spatial dataset
div_sex_prev <- svyby(as.formula(paste0("~", target)), 
                      ~division + sex, 
                      design = bdhs_design, 
                      svymean, 
                      na.rm = TRUE)

# 2. Prepare Data for Mapping
div_sex_map_data <- as.data.frame(div_sex_prev) %>%
    mutate(
        Prevalence = !!sym(target) * 100,
        Label_Text = paste0(round(Prevalence, 1), "%"),
        # Map Division names for the GeoJSON join
        Division_Name = case_when(
            division == "1" ~ "Barisal",
            division == "2" ~ "Chittagong",
            division == "3" ~ "Dhaka",
            division == "4" ~ "Khulna",
            division == "5" ~ "Mymensingh",
            division == "6" ~ "Rajshahi",
            division == "7" ~ "Rangpur",
            division == "8" ~ "Sylhet"
        )
    )

# 3. Join with GeoJSON
bgd_shp <- st_read("small_bangladesh_geojson_adm1_8_divisions_bibhags.json")
map_faceted <- bgd_shp %>%
    left_join(div_sex_map_data, by = c("ADM1_EN" = "Division_Name"))

# 4. Create the Faceted Map
final_faceted_map <- ggplot(data = map_faceted) +
    geom_sf(aes(fill = Prevalence), color = "white", size = 0.2) +
    # Add values inside
    geom_sf_text(aes(label = Label_Text), color = "white", fontface = "bold", size = 2.5) +
    # This creates the two side-by-side maps
    facet_wrap(~sex) + 
    scale_fill_viridis_c(
        option = "magma", 
        direction = -1, 
        name = "Prevalence",
        labels = label_number(suffix = "%")
    ) +
    theme_void() +
    labs(
        title = "Geographic Distribution of Hypertension by Gender",
        subtitle = "BDHS 2022: Comparative Spatial Analysis",
        caption = "Note: Both maps use the same scale for direct comparison."
    ) +
    theme(
        strip.text = element_text(size = 12, face = "bold"), # Title for each facet
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
    )

# 5. Export
ggsave("Figure_Map_Faceted_Sex.png", final_faceted_map, width = 12, height = 8, dpi = 600)











# ==============================================================
# SPATIAL ANALYSIS: FACETED MAP (MALE vs FEMALE)
# ==============================================================
library(sf)
library(survey)
library(ggplot2)
library(dplyr)
library(scales)

# 1. Calculate Weighted Prevalence Stratified by Division AND Sex
div_sex_prev <- svyby(as.formula(paste0("~", target)), 
                      ~division + sex, 
                      design = bdhs_design, 
                      svymean, 
                      na.rm = TRUE)

# 2. Prepare Data: Rename Division Codes AND Sex Labels
div_sex_map_data <- as.data.frame(div_sex_prev) %>%
    mutate(
        Prevalence = !!sym(target) * 100,
        Label_Text = paste0(round(Prevalence, 1), "%"),
        
        # RENAME SEX: Assuming 1=Male, 2=Female (Standard BDHS coding)
        # Adjust if your dataset already uses "Male"/"Female" strings
        sex = case_when(
            sex == "1" ~ "Male",
            sex == "2" ~ "Female",
            TRUE ~ as.character(sex)
        ),
        
        # RENAME DIVISIONS for GeoJSON matching
        Division_Name = case_when(
            division == "1" ~ "Barisal",
            division == "2" ~ "Chittagong",
            division == "3" ~ "Dhaka",
            division == "4" ~ "Khulna",
            division == "5" ~ "Mymensingh",
            division == "6" ~ "Rajshahi",
            division == "7" ~ "Rangpur",
            division == "8" ~ "Sylhet"
        )
    )

# 3. Join with GeoJSON
bgd_shp <- st_read("small_bangladesh_geojson_adm1_8_divisions_bibhags.json")
map_faceted <- bgd_shp %>%
    left_join(div_sex_map_data, by = c("ADM1_EN" = "Division_Name"))

# 4. Create the Faceted Map with custom Colors (Yellow to Red)
final_faceted_map <- ggplot(data = map_faceted) +
    geom_sf(aes(fill = Prevalence), color = "black", size = 0.2) +
    
    # Add values inside each division
    geom_sf_text(aes(label = Label_Text), 
                 color = "black", 
                 fontface = "bold", 
                 size = 3) +
    
    # Split the map into two: Male and Female
    facet_wrap(~sex) + 
    
    # Professional Color Scale: Yellow -> Orange -> Red
    scale_fill_distiller(
        palette = "YlOrRd", 
        direction = 1, 
        name = "Prevalence (%)",
        labels = label_number(suffix = "%")
    ) +
    
    theme_void() +
    labs(
        title = "Geographic Distribution of Hypertension by Gender",
        subtitle = "BDHS 2022: Comparative Analysis (Male vs. Female)",
        caption = "Note: Point labels represent survey-weighted prevalence (%) per division."
    ) +
    theme(
        strip.text = element_text(size = 14, face = "bold", color = "darkblue"), # "Male" / "Female" titles
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 10))
    )

# 5. Export for Publication
ggsave("Figure_Map_Faceted_Sex_YlOrRd.png", final_faceted_map, width = 12, height = 8, dpi = 600)


# ==============================================================
# SPATIAL ANALYSIS: FACETED MAP (MALE vs FEMALE) WITH LABELS
# ==============================================================
library(sf)
library(survey)
library(ggplot2)
library(dplyr)
library(scales)

# 1. Calculate Weighted Prevalence Stratified by Division AND Sex
# Note: Assumes 'bdhs_design' and 'target' are already defined in your environment
div_sex_prev <- svyby(as.formula(paste0("~", target)), 
                      ~division + sex, 
                      design = bdhs_design, 
                      svymean, 
                      na.rm = TRUE)

# 2. Prepare Data: Rename Codes and Create Combined Labels
div_sex_map_data <- as.data.frame(div_sex_prev) %>%
    mutate(
        Prevalence = !!sym(target) * 100,
        
        # RENAME SEX: 1=Male, 2=Female (Standard BDHS coding)
        sex = case_when(
            sex == "1" ~ "Male",
            sex == "2" ~ "Female",
            TRUE ~ as.character(sex)
        ),
        
        # RENAME DIVISIONS for GeoJSON matching
        Division_Name = case_when(
            division == "1" ~ "Barisal",
            division == "2" ~ "Chittagong",
            division == "3" ~ "Dhaka",
            division == "4" ~ "Khulna",
            division == "5" ~ "Mymensingh",
            division == "6" ~ "Rajshahi",
            division == "7" ~ "Rangpur",
            division == "8" ~ "Sylhet"
        ),
        
        # Combine Division Name and Percentage for the label
        Label_Full = paste0(Division_Name, "\n", round(Prevalence, 1), "%")
    )

# 3. Join with GeoJSON
# Ensure this file is in your working directory
bgd_shp <- st_read("small_bangladesh_geojson_adm1_8_divisions_bibhags.json")

map_faceted <- bgd_shp %>%
    left_join(div_sex_map_data, by = c("ADM1_EN" = "Division_Name"))

# 4. Create the Faceted Map (Yellow to Red Palette)
final_faceted_map <- ggplot(data = map_faceted) +
    # Draw the polygons
    geom_sf(aes(fill = Prevalence), color = "grey20", size = 0.2) +
    
    # Add combined labels (Division Name + %)
    geom_sf_text(aes(label = Label_Full), 
                 color = "black", 
                 fontface = "bold", 
                 size = 2.8,
                 lineheight = 0.9) +
    
    # Split the map into two panels by Gender
    facet_wrap(~sex) + 
    
    # Color Scale: YlOrRd (Yellow-Orange-Red)
    scale_fill_distiller(
        palette = "YlOrRd", 
        direction = 1, 
        name = "Prevalence (%)",
        labels = label_number(suffix = "%")
    ) +
    
    # Styling and Labels
    theme_void() +
    labs(
        title = "Geographic Distribution of Hypertension by Gender",
        subtitle = "BDHS 2022: Comparative Analysis (Male vs. Female)",
        caption = "Data Source: Bangladesh Demographic and Health Survey (BDHS) 2022\nNote: Labels represent survey-weighted point estimates."
    ) +
    theme(
        strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(b = 10)),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 10)),
        plot.background = element_rect(fill = "white", color = NA)
    )

# 5. Export for Publication
# High resolution (600 DPI) for clarity in text labels
ggsave("Figure_Map_Faceted_Sex_Full_Labels.png", 
       final_faceted_map, 
       width = 12, 
       height = 8, 
       dpi = 600)

# Display result
print(final_faceted_map)


#
#=========


# Data frame based on the BDHS stats you provided
hypertension_trends <- data.frame(
    Year = c(2011, 2014, 2018, 2022),
    Prevalence = c(18.0, 19.5, 20.4, 20.9),
    Lower_CI = c(17.2, 18.7, 19.6, 19.8),
    Upper_CI = c(18.8, 20.3, 21.2, 22.0)
)

library(ggplot2)

ggplot(hypertension_trends, aes(x = Year, y = Prevalence)) +
    geom_line(size = 1.2, color = "#d73027") +
    geom_point(size = 4, color = "#d73027") +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, color = "black") +
    geom_text(aes(label = paste0(Prevalence, "%")), vjust = -1.5, fontface = "bold") +
    scale_y_continuous(limits = c(15, 25)) +
    theme_minimal() +
    labs(title = "Trends in Hypertension Prevalence in Bangladesh",
         subtitle = "Data Source: BDHS 2011 - 2022",
         y = "Prevalence (%)", x = "Survey Year") +
    theme(plot.title = element_text(face = "bold", size = 14))

ggsave("Figure_Hypertension_Trend.png", width = 8, height = 5, dpi = 300)



# Final Justified Data Frame for your Thesis
htn_overall_justified <- data.frame(
    Year = c(2011, 2014, 2018, 2022),
    Overall = c(18.0, 19.5, 20.4, 20.9),
    Lower_CI = c(17.2, 18.7, 19.6, 19.9),
    Upper_CI = c(18.8, 20.3, 21.2, 21.9)
)

library(ggplot2)

ggplot(htn_overall_justified, aes(x = Year, y = Overall)) +
    geom_area(fill = "#3182bd", alpha = 0.3) +
    geom_line(color = "#08519c", size = 1.5) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, alpha = 0.5) +
    geom_point(size = 4, color = "#08519c") +
    geom_text(aes(label = paste0(Overall, "%")), vjust = -1.5, fontface = "bold") +
    scale_y_continuous(limits = c(0, 25)) +
    theme_minimal() +
    labs(title = "Justification of the Rising Hypertension Burden",
         subtitle = "Crude Prevalence Trend (Adults 18+) | Source: BDHS 2011-2022",
         y = "Prevalence (%)", x = "Survey Year")