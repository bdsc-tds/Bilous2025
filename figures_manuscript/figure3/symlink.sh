#!/bin/bash

# --- Configuration ---
# Get the absolute path to the directory where this script resides
SCRIPT_REAL_PATH=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_REAL_PATH")

# --- DEFINE THE BASE DIRECTORY FOR FIGURES ---
FIGURES_BASE_DIR="/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/figures/" # <--- EDIT THIS LINE

# --- DEFINE YOUR RELATIVE TARGET PATHS HERE ---
RELATIVE_PATHS=(
  # umaps lung 5um
  "embed_panel/10x_5um/NSCLC/lung/lognorm/umap_data_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "resolvi_embed_panel/10x_5um/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "resolvi_supervised_embed_panel/10x_5um/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "ovrlpy_correction_signal_integrity_threshold=0.5_embed_panel/10x_5um/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "ovrlpy_correction_signal_integrity_threshold=0.7_embed_panel/10x_5um/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "split_fully_purified_embed_panel/10x_5um/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  
  # umaps 5k 5um
  "embed_panel/10x_5um/NSCLC/5k/lognorm/umap_data_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "resolvi_embed_panel/10x_5um/NSCLC/5k/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "resolvi_supervised_embed_panel/10x_5um/NSCLC/5k/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "split_fully_purified_embed_panel/10x_5um/NSCLC/5k/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"

  # umaps lung proseg
  "embed_panel/proseg_expected/NSCLC/lung/lognorm/umap_data_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "resolvi_embed_panel/proseg_expected/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "resolvi_supervised_embed_panel/proseg_expected/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "ovrlpy_correction_signal_integrity_threshold=0.5_embed_panel/proseg_expected/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "ovrlpy_correction_signal_integrity_threshold=0.7_embed_panel/proseg_expected/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "split_fully_purified_embed_panel/proseg_expected/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  
  # umaps 5k proseg
  "embed_panel/proseg_expected/NSCLC/5k/lognorm/umap_data_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "resolvi_embed_panel/proseg_expected/NSCLC/5k/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "resolvi_supervised_embed_panel/proseg_expected/NSCLC/5k/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "split_fully_purified_embed_panel/proseg_expected/NSCLC/5k/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"


  # bio/batch conservation barplots
  "scib_metrics_panel_plot/NSCLC/5k/lognorm/scib_metrics_Leiden_ARI_data_matched_reference_combo_rctd_class_aware_Level2.1_n_comps=50_max_n_cells=100000.png"
  "scib_metrics_panel_plot/NSCLC/lung/lognorm/scib_metrics_Leiden_ARI_data_matched_reference_combo_rctd_class_aware_Level2.1_n_comps=50_max_n_cells=100000.png"

  "scib_metrics_panel_plot/NSCLC/5k/lognorm/scib_metrics_iLISI_data_matched_reference_combo_rctd_class_aware_Level2.1_n_comps=50_max_n_cells=100000.png"
  "scib_metrics_panel_plot/NSCLC/lung/lognorm/scib_metrics_iLISI_data_matched_reference_combo_rctd_class_aware_Level2.1_n_comps=50_max_n_cells=100000.png"

  # contamination boxplots
  "contamination_metrics_diffexpr_radius=10_n_permutations=30_n_repeats=5_top_n=20_f1_specificity_boxplot/NSCLC/lung/lognorm/data_matched_reference_combo_rctd_class_aware_Level2.1/lung_T cell_contaminated_by_malignant cell_-log10pvals_x_logfoldchanges_n_hits_top_n=20.png"
  "contamination_metrics_diffexpr_radius=10_n_permutations=30_n_repeats=5_top_n=20_f1_specificity_boxplot/NSCLC/5k/lognorm/data_matched_reference_combo_rctd_class_aware_Level2.1/5k_T cell_contaminated_by_malignant cell_-log10pvals_x_logfoldchanges_n_hits_top_n=20.png"
  
  # summary stats
  "contamination_metrics_diffexpr_radius=10_n_permutations=30_n_repeats=5_top_n=20_f1_sensitivity_boxplot/NSCLC/lung/lognorm/data_matched_reference_combo_rctd_class_aware_Level2.1/lung_median_n_genes.png"
  "contamination_metrics_diffexpr_radius=10_n_permutations=30_n_repeats=5_top_n=20_f1_sensitivity_boxplot/NSCLC/lung/lognorm/data_matched_reference_combo_rctd_class_aware_Level2.1/lung_n_cells.png"
  "contamination_metrics_diffexpr_radius=10_n_permutations=30_n_repeats=5_top_n=20_f1_sensitivity_boxplot/NSCLC/5k/lognorm/data_matched_reference_combo_rctd_class_aware_Level2.1/5k_median_n_genes.png"
  "contamination_metrics_diffexpr_radius=10_n_permutations=30_n_repeats=5_top_n=20_f1_sensitivity_boxplot/NSCLC/5k/lognorm/data_matched_reference_combo_rctd_class_aware_Level2.1/5k_n_cells.png"
  
)

# --- Define Max Filename Length (Common limit, adjust if needed for your filesystem) ---
MAX_FILENAME_LEN=250 # Use 255 for many Linux filesystems, 250 is safer

echo "Script directory (where links will be created): $SCRIPT_DIR"
echo "Figures base directory: $FIGURES_BASE_DIR"
echo "Max link name length: $MAX_FILENAME_LEN"
echo "---"

# --- Loop and Create Links ---
for RELATIVE_PATH in "${RELATIVE_PATHS[@]}"; do
  # Construct the full absolute path to the target file/directory
  FULL_TARGET_PATH="${FIGURES_BASE_DIR%/}/${RELATIVE_PATH#/}"

  echo "Processing relative target: '$RELATIVE_PATH'"

  # --- Initial Checks ---
   # Handle edge case: Empty relative path string
  if [ -z "$RELATIVE_PATH" ]; then
      echo "  Skipping: Relative path string is empty."
      echo "---"
      continue
  fi

   # Check if the *actual* target file exists using the FULL path
  if [ ! -e "$FULL_TARGET_PATH" ]; then
    echo "  Error: Target path '$FULL_TARGET_PATH' does not exist. Skipping link creation."
    echo "---"
    continue
  fi

  # --- Generate and Select Link Name ---
  LINK_NAME="" # Initialize link name

  # 1. Try the default method (replace '/' with '_')
  DEFAULT_LINK_NAME="${RELATIVE_PATH//\//_}"
  DEFAULT_LINK_NAME="${DEFAULT_LINK_NAME#_}" # Remove potential leading '_'

  # Check if the default name is valid and within length limits
  if [ -n "$DEFAULT_LINK_NAME" ] && [ ${#DEFAULT_LINK_NAME} -le $MAX_FILENAME_LEN ]; then
      LINK_NAME="$DEFAULT_LINK_NAME"
      # echo "  Using default link name: '$LINK_NAME'" # Optional: uncomment for debugging
  else
      # Default name is too long or invalid, try fallback (basename)
      if [ -n "$DEFAULT_LINK_NAME" ]; then # Only print length warning if it was generated
         echo "  Warning: Default name ('${DEFAULT_LINK_NAME:0:50}...') too long (${#DEFAULT_LINK_NAME} chars > $MAX_FILENAME_LEN)."
      fi

      FALLBACK_LINK_NAME=$(basename "$RELATIVE_PATH")

      # Check if basename is valid
      if [ -z "$FALLBACK_LINK_NAME" ] || [ "$FALLBACK_LINK_NAME" == "." ] || [ "$FALLBACK_LINK_NAME" == ".." ]; then
          echo "  Error: Could not generate a valid fallback link name (basename) for '$RELATIVE_PATH'. Skipping."
          echo "---"
          continue
      else
          # Check fallback name length (unlikely to be too long, but good practice)
          if [ ${#FALLBACK_LINK_NAME} -gt $MAX_FILENAME_LEN ]; then
              echo "  Error: Fallback link name (basename) '$FALLBACK_LINK_NAME' is *still* too long (${#FALLBACK_LINK_NAME} chars). Skipping."
              echo "---"
              continue
          fi
          LINK_NAME="$FALLBACK_LINK_NAME"
          echo "  Using fallback link name (basename): '$LINK_NAME'"
      fi
  fi

  # Final safety check in case LINK_NAME is somehow empty
  if [ -z "$LINK_NAME" ]; then
      echo "  Error: Failed to determine a valid link name for '$RELATIVE_PATH'. Skipping."
      echo "---"
      continue
  fi

  # Construct the full path for the new symbolic link within the script's directory
  LINK_PATH="$SCRIPT_DIR/$LINK_NAME"

  # --- Collision Check ---
  if [ -e "$LINK_PATH" ]; then
    # Check if it's already a symlink pointing to the correct FULL target path
    # This check works whether the default or fallback name was used.
    if [ -L "$LINK_PATH" ] && [ "$(readlink "$LINK_PATH")" == "$FULL_TARGET_PATH" ]; then
      echo "  Skipping: Link '$LINK_NAME' already exists and points correctly."
    else
      # This is a real collision (same link name points elsewhere or is not a link)
      echo "  Skipping: A file or different link named '$LINK_NAME' already exists."
      # Add extra context if the collision happened with a fallback name
      if [ "$LINK_NAME" != "$DEFAULT_LINK_NAME" ]; then # Check if fallback name was used
            echo "            (Note: This collision occurred using the fallback basename. Check for duplicate filenames in RELATIVE_PATHS)"
      fi
      echo "            Existing item path: $LINK_PATH"
      echo "            Current target was: $FULL_TARGET_PATH"
    fi
    echo "---"
    continue # Skip to the next path in the list
  fi

  # --- Create Symlink ---
  echo "  Creating link: '$LINK_PATH' -> '$FULL_TARGET_PATH'"
  ln -s "$FULL_TARGET_PATH" "$LINK_PATH"

  # --- Verification ---
  if [ $? -eq 0 ]; then
    echo "  Successfully created symbolic link."
  else
    echo "  Error: Failed to create symbolic link '$LINK_NAME'. Check permissions or filesystem limits."
    echo "         Target was: '$FULL_TARGET_PATH'"
  fi
  echo "---" # Separator for clarity
done

echo "Finished processing all paths."
exit 0