# Tutorial: IAS Hotspot Model - Modelling Trial 2025

This tutorial guides you through the process of running the Invasive Alien Species (IAS) Hotspot Model for the "Modelling Trial 2025". This trial focuses on re-running the model with new data for 80 non-indigenous species (NIS) to map their potential suitable habitats in European coastal waters.

**Note on Workflow:** This tutorial is based on the workflow described in `Modelling_trial_2025/README.md`, which involves custom scripts for data processing and a modeling approach combining MCMC-based variable assessment with Random Forest models. Some scripts mentioned in the README for modeling (e.g., `modelisation.teet.GA.R`) currently contain a `biomod2`-based workflow, which appears to be an alternative or developmental version. This tutorial will focus on the custom workflow utilizing functions from `scripts/SEanalytics.functions2025.r`.

## Prerequisites

1.  **R Environment:** Install R and RStudio.
2.  **R Packages:** Ensure the following R packages are installed. You can install them using `install.packages("package_name")`.
    *   `rgbif`
    *   `CoordinateCleaner`
    *   `speciesgeocodeR`
    *   `sf`
    *   `raster`
    *   `terra`
    *   `dplyr`
    *   `readr`
    *   `ggplot2`
    *   `maps`
    *   `glue`
    *   `randomForest`
    *   `fBasics`
    *   `devtools` (for `biooracler` if not on CRAN)
    *   `biooracler` (install using `devtools::install_github("bio-oracle/biooracler")`)
    *   `rnaturalearth`
    *   `rnaturalearthdata`
    *   `ncdf4`
    *   `RNetCDF`
    *   `readxl`
    *   `parallel` (mentioned in README for MCMC)
3.  **Project Structure:** Clone or download this repository and maintain its directory structure.
4.  **GBIF Credentials:** You will need a GBIF account (username, password) and an associated email address.

## Step 1: Configuration

1.  **Set Working Directory:** Open RStudio and set your working directory to the root of the `Modelling_trial_2025` folder or ensure all script paths are correctly pointing to your local repository structure.
2.  **Edit Script Paths (Crucial):** Many scripts contain hardcoded paths like `path <- "~/Dokument/Projekt/HAV2025/"`. You **must** review and update these paths in all relevant scripts to point to the correct locations in your cloned repository.
    *   Key scripts to check for paths:
        *   `scripts/gbif.download.2025.r`
        *   `scripts/gbif.download.check.2025.r`
        *   `scripts/access biooracle.r`
        *   `scripts/arrange.datalayers.2025.r`
        *   `scripts/focal_interpolation.2025.R`
        *   (And any master scripts you create or use for modeling/projection steps)
3.  **GBIF Credentials:**
    *   Open `scripts/gbif.download.2025.r`.
    *   Update the `user`, `pwd`, and `email` variables with your GBIF credentials.
    ```R
    # In scripts/gbif.download.2025.r
    user <- "your_gbif_username"
    pwd <- "your_gbif_password"
    email <- "your_gbif_email"
    ```

## Step 2: Acquire Species Occurrence Data

1.  **Input Species List:** The list of species to model is defined in `input data/NIS_list_combined_Mar2025_v2.csv`. This file is used by `scripts/gbif.download.2025.r` to fetch Taxon Keys.
2.  **Run GBIF Download Script:**
    *   Execute `scripts/gbif.download.2025.r`.
    ```R
    # In R console, assuming working directory is .../Modelling_trial_2025/
    source("scripts/gbif.download.2025.r")
    ```
    *   This script initiates a download request to GBIF. GBIF will process this and send an email with a download link for a ZIP file.
3.  **Retrieve and Place Data:**
    *   Download the ZIP file from the link in your email.
    *   Extract the CSV file (usually named something like `occurrence.txt` or a numeric code `.csv`) from the ZIP.
    *   Place this CSV file into the `Modelling_trial_2025/input data/species.data/` directory (create `species.data` if it doesn't exist, though the scripts seem to use `data/species.data/` relative to the `path` variable).
    *   **Important:** Note the exact filename of this downloaded CSV.
4.  **Update Filename in Check Script:**
    *   Open `scripts/gbif.download.check.2025.r`.
    *   Find the line: `file <- "data/species.data/0018230-250310093411724.csv"`
    *   Change the filename to match the CSV file you just downloaded and placed. For example:
        ```R
        # In scripts/gbif.download.check.2025.r
        file <- "data/species.data/your_downloaded_gbif_data.csv"
        ```

## Step 3: Acquire and Prepare Environmental Data Layers

This step involves downloading data from Bio-Oracle, performing interpolation to fill coastal gaps, and creating raster stacks.

1.  **Define Variables and Scenarios:** Environmental layers are specified in `data/dataset och variabler 2025.xlsx` (relative to the main project path variable). This file is used by `scripts/access biooracle.r`.
2.  **Download Bio-Oracle Layers:**
    *   Execute `scripts/access biooracle.r`.
    ```R
    source("scripts/access biooracle.r")
    ```
    *   This will download the specified layers as `.nc` files into subdirectories like `data/Biooracle.download/datalayer.nc/SCENARIO_NAME/`.
3.  **Perform Focal Interpolation (Baseline Scenario First):**
    *   The script `scripts/focal_interpolation.2025.R` masks land and interpolates coastal data gaps. It's crucial for generating the `filled_layers_new.tif` used in species data cleaning.
    *   Currently, this script is set to process only the first scenario (usually 'baseline'). Ensure 'baseline' is the first in `dataset_scenarios` within the script or modify the script to explicitly process 'baseline'.
    ```R
    # In R console
    source("scripts/focal_interpolation.2025.R")
    ```
    *   **Output:** This creates `data/Biooracle.download/rasterstacks/baseline/filled_layers_new.tif`, among other files.
4.  **Interpolate Other Scenarios (for projections):**
    *   To generate interpolated environmental layers for future scenarios (e.g., SSP119, SSP245), you'll need to modify `scripts/focal_interpolation.2025.R` to loop through all desired `dataset_scenarios` or run it multiple times, changing `sel.sen` each time.
    *   The output for each scenario will be placed in `data/Biooracle.download/rasterstacks/SCENARIO_NAME/`.

## Step 4: Clean Species Data and Generate Pseudoabsences

Now that you have the raw species occurrences and the baseline interpolated environmental data, you can clean the species data.

1.  **Inputs for Cleaning Script:**
    *   Your downloaded and renamed GBIF data CSV (path updated in the script).
    *   `data/Biooracle.download/rasterstacks/baseline/filled_layers_new.tif` (for marine filtering).
    *   `input data/NIS_list_combined_Mar2025_v2.csv` (for species categorization).
2.  **Run Cleaning and Preparation Script:**
    ```R
    source("scripts/gbif.download.check.2025.r")
    ```
3.  **Outputs:**
    *   Individual cleaned species CSV files in `data/Indata2025/speciesIndata2025/`.
    *   `data/species.data/data.table.mars2025.csv`: A crucial manifest file listing species, their data files, and pseudoabsence files. This will be the main input for the modeling loop.
    *   Pseudoabsence CSV files (e.g., `pseudoabsences.marine.excludeboxCATEGORY.csv`) in `data/Indata2025/speciesIndata2025/`.
    *   Various `.rda` files and diagnostic plots.

## Step 5: Species Distribution Modeling (Custom Workflow)

This step involves iterating through each species listed in `data.table.mars2025.csv` and performing the modeling using functions from `SEanalytics.functions2025.r`.

**Note:** The `Modelling_trial_2025/README.md` implies scripts like `modelisation.teet.GA.R` drive this, but their content is `biomod2`-based. You will likely need to create or use a **master R script** that implements the loop and calls the functions from `SEanalytics.functions2025.r` as described below and in the README.

**Conceptual Master Script (`run_all_species_modeling.R` - You may need to create this):**

```R
# run_all_species_modeling.R (Conceptual - adapt paths and create)
library(raster)
library(randomForest)
# ... other necessary libraries

source("scripts/SEanalytics.functions2025.r") # Load helper functions

# Define base paths (ensure these are correct and folders exist)
base_path <- "~/Dokument/Projekt/HAV2025/" # Or your actual base path
indata_path_for_extracted_env_data <- paste0(base_path, "data/Indata2025/extracted_env_data/") 
iterations_path <- paste0(base_path, "data/Indata2025/iterations_cv/")
model_output_path <- paste0(base_path, "results/models_custom/")
roc_plot_path <- paste0(base_path, "results/Plots_ROC_curves_custom/")
prediction_map_path <- paste0(base_path, "results/Rastermaps_custom/current/")
prediction_plot_path <- paste0(base_path, "results/Plots_distribution_custom/current/")
mcmc_plot_path <- paste0(base_path, "results/Plots_variables_custom/")

# Create output directories if they don't exist
dir.create(indata_path_for_extracted_env_data, recursive = TRUE, showWarnings = FALSE)
dir.create(iterations_path, recursive = TRUE, showWarnings = FALSE)
dir.create(model_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(roc_plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(prediction_map_path, recursive = TRUE, showWarnings = FALSE)
dir.create(prediction_plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(mcmc_plot_path, recursive = TRUE, showWarnings = FALSE)


# Load the main species data table
species_table_file <- paste0(base_path, "data/species.data/data.table.mars2025.csv")
data_table <- read.csv2(species_table_file, stringsAsFactors = FALSE) # Use read.csv2 if it's semicolon separated

# Load the baseline environmental raster stack (interpolated)
env_stack_file <- paste0(base_path, "data/Biooracle.download/rasterstacks/baseline/Biooracle.filled.layers.Europe2025.tif") # Or global
env_stack <- stack(env_stack_file)

# Define species data path
species_data_path <- paste0(base_path, "data/Indata2025/speciesIndata2025/")


# --- Loop through species ---
for (i in 1:nrow(data_table)) {
    species_name <- data_table$species[i]
    print(paste("Processing species:", species_name))

    # 1. Read and Extract Environmental Data
    # Note: read.and.extract in SEanalytics.functions2025.r uses slightly different pathing logic
    # It expects speciespath, stackpath, plotpath, outpath.
    # Here, we're simplifying to assume it returns the data and we save it.
    extracted_data_list <- read.and.extract(
        data.table = data_table,
        species = species_name,
        stack = env_stack, # Pass the loaded stack object
        speciespath = species_data_path,
        stackpath = dirname(env_stack_file), # Path to the stack directory
        plotpath = prediction_plot_path, # General plot path
        outpath = indata_path_for_extracted_env_data # Where to save outputs if function does so
    )
    
    complete_points <- extracted_data_list$complete.points
    species_stats <- extracted_data_list$stats
    
    # Save the extracted data for the species
    write.csv(complete_points, paste0(indata_path_for_extracted_env_data, species_name, "_env_data.csv"), row.names = FALSE)
    print(paste("Saved extracted env data for", species_name))

    # Filter out species with too few occurrences (as per README)
    if (species_stats[5] < 5) { # npos.unique.complete < 5
        print(paste("Skipping", species_name, "due to insufficient (<5) unique complete positive findings."))
        next
    }

    # 2. Prepare Cross-Validation Folds
    # This function saves .rda files to iterations.path
    cv_site_counts <- split.data(
        species = species_name,
        indata.path = indata_path_for_extracted_env_data, # Path where species_env_data.csv is
        iterations.path = iterations_path
    )
    print(paste("CV folds prepared for", species_name, "- Positive sites:", cv_site_counts[1], "Negative sites:", cv_site_counts[2]))
    
    # Filter out species where CV splitting might fail (as per README criteria: <5 distinct areas)
    if (cv_site_counts[1] < 5) { # Assuming pos.sites < 5
         print(paste("Skipping", species_name, "due to insufficient (<5) distinct positive areas for CV."))
         next
    }

    # 3. MCMC Feature Selection (for variable assessment)
    # The README mentions 25 iterations. MCMC.process is complex and might need parallel setup.
    # This is a simplified placeholder for the concept.
    # Actual implementation would require careful setup of MCMC.process calls.
    print(paste("Performing MCMC variable assessment for", species_name, "(conceptual step)"))
    # RI_results <- list()
    # for (m_iter in 1:25) {
    #    # Subset data or use CV folds appropriately for MCMC.process
    #    # mcmc_output <- MCMC.process(...) 
    #    # RI_results[[m_iter]] <- mcmc_output$RI
    # }
    # Aggregate RI_results and plot variable importance (e.g., boxplots of RI scores)
    # Save plot to mcmc_plot_path
    
    # 4. Train Random Forest Models (CV and Final)
    # run.random.forests saves the final model to an .rda file.
    # It expects indata.path to point to where species_env_data.csv is.
    rf_model_list <- run.random.forests(
        species = species_name,
        selvar = "all", # Use all variables as per README
        indata.path = indata_path_for_extracted_env_data,
        iterations.path = iterations_path
    )
    
    # Save the full model explicitly if not done by run.random.forests in a known way
    final_model <- rf_model_list$RF.results.alldata$RF.selected
    save(final_model, file = paste0(model_output_path, species_name, "_final_model.rda"))
    print(paste("Trained and saved final RF model for", species_name))

    # 5. Calculate and Plot ROC Curves
    # Prepare true classes for ROC calculation from the CV test sets
    load(paste0(iterations_path, "/occurance.iters_", species_name, ".rda")) # loads all.occurance.iters
    true_classes_all_reps <- list()
    
    for (rep_idx in 1:length(all.occurance.iters)) {
        true_classes_one_rep <- c()
        for (fold_idx in 1:length(all.occurance.iters[[rep_idx]])) {
            test_ids <- all.occurance.iters[[rep_idx]][[fold_idx]]
            test_data_subset <- complete_points[complete_points$ID %in% test_ids, ]
            true_classes_one_rep <- c(true_classes_one_rep, setNames(as.character(test_data_subset$occurrenceStatus), test_data_subset$ID))
        }
        true_classes_all_reps[[rep_idx]] <- true_classes_one_rep
    }

    all_roc_stats <- list()
    mean_auc_values <- numeric(length(rf_model_list$RF.results.CV))
    
    for (k_rep in 1:length(rf_model_list$RF.results.CV)) { # For each repetition of 5-fold CV
        roc_output <- calc.ROC(
            .RF.result = rf_model_list$RF.results.CV[[k_rep]], # Pass one full CV set (5 folds)
            .true.class = true_classes_all_reps[[k_rep]] # True classes for this repetition
        )
        all_roc_stats[[k_rep]] <- roc_output
        mean_auc_values[k_rep] <- roc_output$AUC
    }
    
    # Create a mean ROC from all_roc_stats (average sensitivities at specificities)
    # This is a simplified placeholder; actual averaging needs care.
    mean_roc_for_plotting <- all_roc_stats[[1]] # Placeholder
    
    plot.ROC(
        .ROC.path = roc_plot_path,
        .my.species = species_name,
        .mean.ROC = mean_roc_for_plotting, # This should be properly averaged
        .all.ROC = all_roc_stats, # list of ROC outputs from each rep
        method = "allvars" # as per SEanalytics default
    )
    print(paste("Generated ROC plots for", species_name, "Mean AUCs:", paste(round(mean_auc_values,3), collapse=", ")))

    # 6. Spatial Prediction (Current Conditions)
    # predict.maps needs the path to the saved model .rda file
    # The function SEanalytics.functions2025.r::predict.maps expects a model *path*, not object.
    # It loads a variable 'rf.output' which should contain $RF.selected.
    # We need to ensure the saved model RDA is in the format predict.maps expects or adapt.
    # For simplicity, let's assume we pass the loaded model object and env_stack directly,
    # or we modify predict.maps. The original predict.maps loads 'bar' (env stack) itself.
    
    # To use predict.maps as is, we'd save rf.output correctly:
    rf.output <- rf_model_list$RF.results.alldata 
    save(rf.output, file=paste0(model_output_path,"/RF.model.and.predictions.eur.wt.",species_name,".rda"))

    current_prediction_map <- predict.maps(
        species = species_name,
        modelpath = model_output_path # Path to where RF.model.and.predictions... is saved
        # `bar` (env_stack) is loaded inside predict.maps from a hardcoded path, which needs to be fixed.
        # For now, let's assume `bar` should be `env_stack` and passed or loaded flexibly.
    )
    # We need to assign `bar` globally or modify predict.maps
    assign("bar", env_stack, envir = .GlobalEnv) # Hacky, better to pass to function

    current_prediction_map_filename <- paste0(prediction_map_path, "linear.prob.", species_name, ".tif")
    writeRaster(current_prediction_map, current_prediction_map_filename, overwrite = TRUE)
    print(paste("Generated current prediction map for", species_name))

    # 7. Plot Current Prediction Maps
    # plot.maps uses the generated raster map file.
    # It also loads species data from indata.path (where species_env_data.csv is)
    # and a shapefile (hardcoded path in function, needs fixing).
    # Define colors and breaks for plotting
    brk <- seq(0, 1, by = 0.01)
    colors <- rev(rainbow(length(brk) - 1, start = 0, end = 0.7))
    world_shape <- st_read(paste0(base_path, "data/ref-countries-2020-01m.shp/CNTR_RG_01M_2020_4326.shp")) # Load shapefile
    
    # Assign shape globally or modify plot.maps
    assign("shape", world_shape, envir = .GlobalEnv) # Hacky

    plot.maps(
        .Species = species_name,
        indata.path = indata_path_for_extracted_env_data, # Path to species_env_data.csv
        .rastermap = current_prediction_map_filename,
        .plotpath = prediction_plot_path,
        .colors = colors,
        .brk = brk,
        .shape = world_shape, # Pass shapefile object
        xlim = c(10, 32), # Sweden X extent for Swe plot
        ylim = c(54, 69)  # Sweden Y extent for Swe plot
    )
    print(paste("Plotted current distribution for", species_name))
    
    # Clean up global assignments if any
    if(exists("bar", envir = .GlobalEnv)) remove("bar", envir = .GlobalEnv)
    if(exists("shape", envir = .GlobalEnv)) remove("shape", envir = .GlobalEnv)

    print(paste("Finished processing for species:", species_name))
}

print("All species processed.")
```

**Explanation of Conceptual Master Script Sections:**
*   **Setup:** Loads `SEanalytics.functions2025.r`, defines paths, loads main species table and baseline environmental stack.
*   **Loop:** Iterates through each `species_name` from `data_table`.
*   **`read.and.extract`:** Prepares `complete_points` by merging species occurrences with environmental data. Skips species with <5 positive unique complete findings.
*   **`split.data`:** Prepares CV folds based on spatial blocking. Skips species if <5 distinct positive areas for robust CV.
*   **MCMC (Conceptual):** Placeholder for the MCMC variable assessment loop described in the README. This would involve running `MCMC.process` multiple times and aggregating RI scores for plotting.
*   **`run.random.forests`:** Trains RF models (for each CV fold and one final model on all data). Saves the final model.
*   **`calc.ROC` & `plot.ROC`:** Calculates ROC metrics from CV results and plots them.
*   **`predict.maps`:** Generates the spatial prediction (probability map) for current conditions using the final model and baseline environmental stack.
*   **`plot.maps`:** Creates PNG plots of the current prediction map with occurrence data overlaid.

## Step 6: Future Projections (Custom Workflow)

Similar to Step 5, you'll need a master script to loop through species and future scenarios, using the trained models and functions from `SEanalytics.functions2025.r`.

**Conceptual Master Script for Projections (`run_all_species_projections.R`):**
```R
# run_all_species_projections.R (Conceptual - adapt paths and create)
library(raster)
# ... other necessary libraries

source("scripts/SEanalytics.functions2025.r")

# Define base paths
base_path <- "~/Dokument/Projekt/HAV2025/" # Or your actual base path
model_input_path <- paste0(base_path, "results/models_custom/") # Where final_model.rda are
future_prediction_map_path_base <- paste0(base_path, "results/Rastermaps_custom/future/")
future_prediction_plot_path_base <- paste0(base_path, "results/Plots_species_custom/future/") # For 4-panel plots

dir.create(future_prediction_map_path_base, recursive = TRUE, showWarnings = FALSE)
dir.create(future_prediction_plot_path_base, recursive = TRUE, showWarnings = FALSE)

# Load the main species data table
species_table_file <- paste0(base_path, "data/species.data/data.table.mars2025.csv")
data_table <- read.csv2(species_table_file, stringsAsFactors = FALSE)

# Define future scenarios and their corresponding environmental data paths
# Assumes focal_interpolation.2025.R has been run for these scenarios
future_scenarios <- list(
    SSP119_2050 = paste0(base_path, "data/Biooracle.download/rasterstacks/ssp119/Biooracle.filled.layers.Europe2025_2050.tif"), # Example path, adapt
    SSP119_2100 = paste0(base_path, "data/Biooracle.download/rasterstacks/ssp119/Biooracle.filled.layers.Europe2025_2100.tif")  # Example path, adapt
    # Add other SSPs and time periods as needed
)
# Ensure these future environmental stack files exist and names match your processing.
# The `focal_interpolation.2025.R` currently creates files like `Biooracle.filled.layers.Europe2025.tif`
# within scenario folders. You might need to rename or point to these, and ensure they represent 2050/2100.

# --- Loop through species ---
for (i in 1:nrow(data_table)) {
    species_name <- data_table$species[i]
    print(paste("Processing future projections for species:", species_name))

    # Load the final model for the species
    model_file <- paste0(model_input_path, species_name, "_final_model.rda")
    if (!file.exists(model_file)) {
        print(paste("Model file not found for", species_name, "- skipping projections."))
        next
    }
    load(model_file) # Loads `final_model` object (or `rf.output` if saved that way)
    
    # Ensure the model object is named `final_model` or adapt below.
    # If `rf.output` was saved, then `final_model <- rf.output$RF.selected`

    # Load current prediction map for comparison plots
    current_map_file <- paste0(base_path, "results/Rastermaps_custom/current/linear.prob.", species_name, ".tif")
    if (!file.exists(current_map_file)){
        print(paste("Current map not found for", species_name, "- 4-panel plots might be incomplete."))
        current_pred_map <- NULL
    } else {
        current_pred_map <- raster(current_map_file)
    }
    
    species_future_maps <- list(current = current_pred_map)

    # --- Loop through future scenarios ---
    for (scenario_name in names(future_scenarios)) {
        env_stack_future_file <- future_scenarios[[scenario_name]]
        if (!file.exists(env_stack_future_file)) {
            print(paste("Future env stack not found for scenario", scenario_name, "- skipping."))
            next
        }
        env_stack_future <- stack(env_stack_future_file)
        
        # Predict for the future scenario
        # Adapt predict.maps or ensure global `bar` is set to env_stack_future
        # For now, conceptual call:
        # future_prediction <- predict(env_stack_future, final_model, type="prob")
        # The actual predict.maps from SEanalytics is more complex.
        # It expects rf.output to be saved in a specific way.
        # Re-saving the model in the expected format for predict.maps:
        rf.output_temp <- list(RF.selected = final_model) # Assuming final_model is the randomForest object
        temp_model_path_for_predict <- paste0(model_input_path, "temp_model_for_predict.rda")
        save(rf.output_temp, file=temp_model_path_for_predict)
        rm(rf.output_temp) # clean up

        assign("bar", env_stack_future, envir = .GlobalEnv) # Hacky
        
        prediction_future_map <- predict.maps(
            species = species_name, # Species name is used for loading specific rf.output
            modelpath = model_input_path # Path where the temp_model_for_predict.rda could be placed, or actual species model
            # The species argument in predict.maps actually loads "RF.model.and.predictions.eur.wt.SPECIES.rda"
            # So, we need to use the original model path and ensure `bar` is the future stack.
        )
        
        future_map_filename <- paste0(future_prediction_map_path_base, scenario_name, "_linear.prob.", species_name, ".tif")
        writeRaster(prediction_future_map, future_map_filename, overwrite = TRUE)
        species_future_maps[[scenario_name]] <- prediction_future_map
        print(paste("Generated future prediction map for", species_name, "-", scenario_name))
        
        if(exists("bar", envir = .GlobalEnv)) remove("bar", envir = .GlobalEnv)
    }

    # Plot 4-panel map (Current, 2050, 2100, Difference) - Conceptual
    # This requires a custom plotting function not explicitly in SEanalytics.functions2025.r
    # but described in the README.
    if (!is.null(species_future_maps$current) && 
        !is.null(species_future_maps$SSP119_2050) && 
        !is.null(species_future_maps$SSP119_2100)) {
        
        png(paste0(future_prediction_plot_path_base, species_name, "_ssp119_comparison.png"), width = 1000, height = 1000)
        par(mfrow = c(2, 2))
        
        plot(species_future_maps$current, main = paste(species_name, "- Current"))
        # Add world boundaries, etc.
        
        plot(species_future_maps$SSP119_2050, main = "SSP119 - 2050")
        # Add world boundaries
        
        plot(species_future_maps$SSP119_2100, main = "SSP119 - 2100")
        # Add world boundaries
        
        # Difference map (e.g., 2100 - current)
        diff_map <- species_future_maps$SSP119_2100 - species_future_maps$current
        plot(diff_map, main = "Difference (2100 - Current)")
        # Add world boundaries
        
        dev.off()
        par(mfrow = c(1, 1))
        print(paste("Generated 4-panel comparison plot for", species_name))
    }
}
print("All species future projections processed.")
```

## Step 7: Aggregate Results and Final Outputs

After individual species modeling and projections:
1.  **Average Probability Maps:** Calculate mean probability across all species for current and future scenarios.
2.  **Species Richness Hotspots:** Sum binary presence/absence maps (thresholded probability maps) to identify areas with high numbers of potential invaders.
3.  **Combine with Traffic Layers:** Overlay hotspot maps with shipping traffic data (external data source) to identify high-risk introduction and establishment zones.

These aggregation steps will require further custom R scripting, potentially building upon the individual GeoTIFF outputs.

## Troubleshooting & Notes

*   **Paths:** Incorrect file paths are the most common source of errors. Double-check all path variables in scripts.
*   **Memory:** Processing large raster stacks and many species can be memory-intensive.
*   **Hardcoded Values:** Be aware of hardcoded filenames, species names, or parameters within the scripts and adjust as needed for your specific run.
*   **`SEanalytics.functions2025.r` Dependencies:** Some functions within `SEanalytics.functions2025.r` might have implicit dependencies on global variables (like `bar` for the raster stack or `shape` for plotting) or expect files to be in very specific locations or formats. These might need adjustments for robust execution. The conceptual master scripts above try to handle some of this by assigning them globally before function calls, which is not ideal but reflects potential script structure.
*   **`biomod2` Alternative:** If the custom workflow proves too difficult to reconstruct from the existing functions, the `biomod2`-based scripts (`modelisation.teet.GA.R`, `projection&forecasting.mod.R`) provide an alternative, well-documented framework for SDM. However, they would need to be adapted to loop through the project's species list and use the centrally prepared data.

This tutorial provides a detailed plan. Some parts, especially the master scripts for steps 5 and 6, are conceptual and would need to be fully implemented by the user, carefully calling the functions from `SEanalytics.functions2025.r` and managing data flow.

