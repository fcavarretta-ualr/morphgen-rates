import pandas as pd

_file_list = {
  'aPC':{
    'SL':{
      'apical':{
        'fname_sholl':'sl_sholl_plot',
        'fname_bif':'sl_bifurcations'
        }
      },
    'PYR':{
      'apical':{
        'fname_sholl':'pyr_sholl_plot',
        'fname_bif':'pyr_bifurcations'
        }
      },
    }
  'NEOC':{
    'PYR':{
      'apical':{
        'fname_sholl':'mitral_sholl_plot',
        'fname_bif':'mitral_bifurcations'
        }
      },
    }
  'OB':{
    'MITRAL':{
      'lateral':{
        'fname_sholl':'mitral_sholl_plot',
        'fname_bif':'mitral_bifurcations'
        }
      },
    'TUFTED':{
      'lateral':{
        'fname_sholl':'tufted_sholl_plot',
        'fname_bif':'tufted_bifurcations'
        }
      },
    }
  }


def _local_data_path(filename, ext="csv"):
  """
  Build a path like: <this_file_dir>/data/<filename>.<ext>

  Parameters
  ----------
  filename : str
      Base filename (without extension)
  ext : str, default "csv"
      File extension (without the dot)

  Returns
  -------
  pathlib.Path
      Full path to the data file
  """
  work_dir = Path(__file__).resolve().parent
  return work_dir / "data" / f"{filename}.{ext}"


def _get_by_path(d, path, sep="/"):
  """
  Retrieve a value from a nested dictionary using a path-like key string.

  Example
  -------
  d = {"a": {"b": {"c": 123}}}
  get_by_path(d, "a/b/c") -> 123

  Parameters
  ----------
  d : dict
      Nested dictionary
  path : str
      Path of keys, e.g. "keys1/keys2/keys3"
  sep : str, default "/"
      Path separator

  Returns
  -------
  object
      The value stored at the given path

  Raises
  ------
  KeyError
      If any key along the path is missing
  TypeError
      If a non-dict is encountered before the final key
  """
  cur = d
  for k in path.split(sep):
    if not isinstance(cur, dict):
      raise TypeError(f"Expected dict at '{k}', got {type(cur).__name__}")
    cur = cur[k]
  return cur


def _get_data(fname_sholl, fname_bif):
  """
  Retrieve a dataset entry using a key-path of the form
  "<brain region>/<neuron class>/<subcellular section>".

  The argument `data_path` is interpreted as a slash-separated path of keys used
  to traverse a nested dataset dictionary. The selected dataset is expected to
  contain both Sholl-plot statistics and bifurcation statistics; when both are
  available, this function returns a standardized dictionary compatible with
  `compute_rates`.

  Parameters
  ----------
  data_path : str
      Dataset identifier expressed as a key path:

      "<brain region>/<neuron class>/<subcellular section>"

      Examples:
      - "CTX/pyr/apical"
      - "HPC/pyr/basal"

      Each component is used as a successive key lookup into the nested dataset
      container.

  Returns
  -------
  dict
      If both Sholl and bifurcation information are present for the selected dataset,
      returns:

      data = {
        "sholl": {
          "bin_size": float,
          "mean": numpy.ndarray,   # shape (K,)
          "var":  numpy.ndarray,   # shape (K,)
        },
        "bifurcations": {
          "mean": float,
          "var":  float,
        },
      }

      Where:
      - `data["sholl"]["bin_size"]` is the spatial bin size used to define Sholl shells
      - `data["sholl"]["mean"]` is the mean Sholl intersection count per radial bin
      - `data["sholl"]["var"]` is the variance of the Sholl intersection count per bin
      - `data["bifurcations"]["mean"]` is the mean bifurcation count
      - `data["bifurcations"]["var"]` is the variance of the bifurcation count

  Raises
  ------
  KeyError
      If any key along `data_path` is missing (brain region, neuron class, or section)
  ValueError
      If the selected dataset does not contain both Sholl and bifurcation data, or
      if the provided arrays have incompatible shapes

  Notes
  -----
  - `data_path` is a *key path*, not a filesystem path
  - The function assumes the dataset entry referenced by `data_path` includes:
      - Sholl bin size, mean array, variance array
      - Bifurcation mean and variance

  Examples
  --------
  >>> data = get("CTX/pyr/apical")
  >>> data["sholl"]["bin_size"]
  50.0
  >>> data["bifurcations"]["mean"]
  12.3
  """
  
  data = {}
  
  # Load Sholl plot summary statistics (bin counts + variance) from CSV
  if fname_sholl:
    df_sholl = pd.read_csv(_local_data_path(fname_sholl), index_col=0)
    # manipulate the data
    df_sholl = df_sholl.T.describe().T[['mean', 'std']]
    df_sholl = df_sholl[(df_sholl != 0).all(axis=1)]
    bin_size = df_sholl.index[1] - df_sholl.index[0]
    df_sholl = df_sholl.to_numpy()

    data['sholl'] = {
      'bin_size':bin_size, 
      'mean':df_sholl[:, 0],
      'var':df_sholl[:, 1] ** 2,
      }

  if fname_bif:
    # Load bifurcation summary statistics from CSV
    df_bif = pd.read_csv(_local_data_path(fname_bif), index_col=0).to_numpy()

    # Bundle inputs exactly as loaded (no preprocessing)
    data["bifurcations"] = {
      'mean':df_bif.mean(),
      'var':df_bif.var()
      }

  return data


def get(data_path):
  return _get_data(**_get_by_path(_file_list, data_path))
  
