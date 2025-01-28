import argparse
import os
import shutil
import pathlib
import tiledbsoma as soma
import tiledbsoma.io


# Function to run Python portion
def soma_to_h5ad(experiment_dir, h5ad=None):
    if h5ad is None:
        h5ad = pathlib.Path(experiment_dir).stem + ".h5ad"

    with soma.open(experiment_dir) as experiment:
        soma.io.to_h5ad(
            experiment=experiment,
            measurement_name="RNA",
            X_layer_name="counts",
            h5ad_path=h5ad,
        )
    print(f"Finished Python tasks: Saved H5AD at {h5ad}")


# Function to delete SOMA directory
def delete_soma(soma):
    if os.path.exists(soma):
        shutil.rmtree(soma)
        print(f"Deleted SOMA object at {soma}")
    else:
        print(f"Warning: SOMA directory {soma} not found.")


# Main function to handle argument parsing and call respective functions
def main():
    parser = argparse.ArgumentParser(description="Convert SOMA object to H5AD format")

    parser.add_argument("--soma", required=True, help="Path to SOMA input folder")
    parser.add_argument(
        "--h5ad", required=False, default=None, help="Output path for H5AD file"
    )
    parser.add_argument(
        "--delete-soma",
        action="store_true",
        help="Delete SOMA object after conversion to H5AD",
    )

    args = parser.parse_args()

    # Run Python script to convert SOMA to H5AD
    soma_to_h5ad(args.soma, args.h5ad)

    # If delete_soma is set, delete the SOMA object
    if args.delete_soma:
        delete_soma(args.soma)


if __name__ == "__main__":
    main()
