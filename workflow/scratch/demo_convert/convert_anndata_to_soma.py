import scanpy as sc
import tiledb
import tiledbsoma
import tiledbsoma.io
import argparse
import pathlib


def h5ad_to_soma(h5ad, soma=None, rm_uns=False):
    if soma is None:
        soma = pathlib.Path(h5ad).stem + ".soma"
    print("Reading anndata...")
    ad = sc.read(h5ad)
    ad.obs_names = ad.obs_names.astype(str)
    ad.obs.index.name = None

    if rm_uns:
        for k in list(ad.uns.keys()):
            del ad.uns[k]

    for c in ad.obs.columns:
        if ad.obs[c].dtype == "category":
            ad.obs[c] = ad.obs[c].astype("str")

    print(f"Saving SOMA at {soma}")
    tiledbsoma.logging.info()
    vfs = tiledb.VFS()
    tiledbsoma.io.from_anndata(experiment_uri=soma, anndata=ad, measurement_name="RNA")
    print(f"Finished Python tasks: Saved SOMA at {soma}")


def main():
    parser = argparse.ArgumentParser(description="Convert H5AD object to SOMA format")

    parser.add_argument("--h5ad", required=True, help="Path to SOMA input folder")
    parser.add_argument(
        "--soma", required=False, default=None, help="Output path for H5AD file"
    )
    parser.add_argument(
        "--rm-uns",
        action="store_true",
        help="Remove uns fields in H5AD object before saving to SOMA",
    )

    args = parser.parse_args()

    # Run Python script to convert SOMA to H5AD
    h5ad_to_soma(args.h5ad, args.soma, args.rm_uns)


if __name__ == "__main__":
    main()
