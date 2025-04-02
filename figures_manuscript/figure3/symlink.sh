#!/bin/bash

# --- Configuration ---
# Get the absolute path to the directory where this script resides
SCRIPT_REAL_PATH=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_REAL_PATH")

# --- DEFINE THE BASE DIRECTORY FOR FIGURES ---
# !!! IMPORTANT: Set this to the ABSOLUTE path of the directory *containing* the 'figures' folder !!!
# Example: If your figures are in /home/user/my_project/figures/..., set this to /home/user/my_project
FIGURES_BASE_DIR="/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/figures/" # <--- EDIT THIS LINE

# --- DEFINE YOUR RELATIVE TARGET PATHS HERE ---
# Paths should be relative to FIGURES_BASE_DIR defined above.
# They should typically start with 'figures/' if FIGURES_BASE_DIR is the parent.
RELATIVE_PATHS=(
  # umaps lung
  "resolvi_embed_panel/10x_5um/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "resolvi_supervised_embed_panel/10x_5um/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "ovrlpy_correction_signal_integrity_threshold=0.5_embed_panel/10x_5um/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "ovrlpy_correction_signal_integrity_threshold=0.7_embed_panel/10x_5um/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "split_fully_purified_embed_panel/10x_5um/NSCLC/lung/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  
  # umaps 5k
  "resolvi_embed_panel/10x_5um/NSCLC/5k/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "resolvi_supervised_embed_panel/10x_5um/NSCLC/5k/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"
  "split_fully_purified_embed_panel/10x_5um/NSCLC/5k/lognorm/umap_n_comps=50_n_neighbors=50_min_dist=0.5_euclidean_matched_reference_combo_rctd_class_aware_Level2.1.png"

  # bio/batch conservation barplots
  "scib_metrics_panel_plot/NSCLC/5k/lognorm/scib_metrics_Leiden_ARI_data_matched_reference_combo_rctd_class_aware_Level2.1_n_comps=50_max_n_cells=100000.png"
  "scib_metrics_panel_plot/NSCLC/lung/lognorm/scib_metrics_Leiden_ARI_data_matched_reference_combo_rctd_class_aware_Level2.1_n_comps=50_max_n_cells=100000.png"

  "scib_metrics_panel_plot/NSCLC/5k/lognorm/scib_metrics_iLISI_data_matched_reference_combo_rctd_class_aware_Level2.1_n_comps=50_max_n_cells=100000.png"
  "scib_metrics_panel_plot/NSCLC/lung/lognorm/scib_metrics_iLISI_data_matched_reference_combo_rctd_class_aware_Level2.1_n_comps=50_max_n_cells=100000.png"

  # contamination boxplots
  
  # summary stats

)

echo "Script directory (where links will be created): $SCRIPT_DIR"
echo "Figures base directory: $FIGURES_BASE_DIR"
echo "---"

# --- Loop and Create Links ---
for RELATIVE_PATH in "${RELATIVE_PATHS[@]}"; do
  # Construct the full absolute path to the target file/directory
  # Remove potential double slashes if FIGURES_BASE_DIR ends with / and RELATIVE_PATH starts with /
  FULL_TARGET_PATH="${FIGURES_BASE_DIR%/}/${RELATIVE_PATH#/}"

  echo "Processing relative target: '$RELATIVE_PATH'"
  # Optional: uncomment to always show the full target path
  # echo "  Computed full target path: '$FULL_TARGET_PATH'"

  # --- Generate Custom Link Name (Based on RELATIVE path) ---
  # 1. Replace all '/' characters with '_'
  SANITIZED_NAME="${RELATIVE_PATH//\//_}"
  # 2. Remove leading '_' (unlikely needed if paths start like 'figures/..', but safe)
  LINK_NAME="${SANITIZED_NAME#_}"
  # 3. Handle edge case: Empty relative path
  if [ -z "$LINK_NAME" ]; then
      echo "  Skipping: Cannot generate a valid link name for empty relative path."
      echo "---"
      continue # Skip to the next path
  fi

  # Construct the full path for the new symbolic link within the script's directory
  LINK_PATH="$SCRIPT_DIR/$LINK_NAME"

  # --- Pre-checks ---
  # Check if the *actual* target file exists using the FULL path
  if [ ! -e "$FULL_TARGET_PATH" ]; then
    # Changed to output an error and skip by default, as dangling links are often undesirable
    echo "  Error: Target path '$FULL_TARGET_PATH' does not exist. Skipping link creation."
    # If you WANT dangling links, comment out the echo/continue below and uncomment the 'ln -s' line
    # echo "  Warning: Target path '$FULL_TARGET_PATH' does not exist. Creating dangling symlink."
    echo "---"
    continue
  fi

  # Check if a file/directory/link with the generated name already exists
  if [ -e "$LINK_PATH" ]; then
    # Check if it's already a symlink pointing to the correct FULL target path
    if [ -L "$LINK_PATH" ] && [ "$(readlink "$LINK_PATH")" == "$FULL_TARGET_PATH" ]; then
      echo "  Skipping: Link '$LINK_NAME' already exists and points correctly."
    else
      echo "  Skipping: A file or different link named '$LINK_NAME' already exists in the script directory."
      echo "            Existing item path: $LINK_PATH"
    fi
    echo "---"
    continue # Skip to the next path in the list
  fi

  # --- Create Symlink ---
  # Point the link TO the FULL absolute target path
  echo "  Creating link: '$LINK_PATH' -> '$FULL_TARGET_PATH'"
  # Use quotes to handle paths/names with spaces or special characters
  ln -s "$FULL_TARGET_PATH" "$LINK_PATH"

  # --- Verification ---
  if [ $? -eq 0 ]; then
    echo "  Successfully created symbolic link."
  else
    # This error might occur if permissions are wrong, etc.
    echo "  Error: Failed to create symbolic link for '$FULL_TARGET_PATH'."
  fi
  echo "---" # Separator for clarity
done

echo "Finished processing all paths."
exit 0