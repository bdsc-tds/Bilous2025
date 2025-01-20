library(Seurat)

library(yaml)
# install.packages("this.path")

script_dir = this.path::here()
config_path = paste0(script_dir,"/../../config/config.yml")

config <- function(path=config_path) {
    # Read the configuration file and return a dictionary of config values.

    # Parameters
    # ----------
    # path : str
    #     The path to the configuration file. Defaults to the value of
    #     `config_path` if not provided.

    # Returns
    # -------
    # cfg : dict
    #     A dictionary of configuration values. All values are strings and
    #     have been converted to absolute paths by prepending the value of
    #     `cfg["base_dir"]`.
    
  cfg <- read_yaml(path)
  for (k in names(cfg)) {
    if (k != "base_dir") {
      cfg[[k]] <- paste0(cfg$base_dir, cfg[[k]])
    }
  }
  return(cfg)
}

soma_to_seurat <- function( soma_input, 
                            counts = "counts", 
                            data = "data", 
                            obs_index="obs_id", 
                            var_index="var_id", 
                            measurement_name = "RNA", 
                            return_experiment = FALSE) {
  library(tiledbsoma)
  # Export a SOMA experiment to a Seurat object.
    
  # Open the SOMA experiment
  experiment <- SOMAExperimentOpen(soma_input)
  
  # Create a new SOMA axis query for the specified measurement
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = measurement_name
  )
  

  X_layers <- c()
  for (name in query$ms$X$names()) {
      X_layers[name] <- name
  }


  # Convert the query to a Seurat object
  seurat_object <- query$to_seurat(
    X_layers = X_layers,
    obs_index = obs_index,
    var_index = var_index
  )
  
  # Optionally return the experiment as well
  if (return_experiment) {
    return(list(seurat_object = seurat_object, experiment = experiment))
  }
  
  # Return only the Seurat object
  return(seurat_object)
}
