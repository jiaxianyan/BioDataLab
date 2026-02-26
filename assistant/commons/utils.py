from joblib import Parallel, delayed, cpu_count
import yaml
from tqdm import tqdm
from typing import *

def pmap_multi(pickleable_fn, data, n_jobs=None, verbose=1, desc=None, **kwargs):
  if n_jobs is None:
    n_jobs = cpu_count() - 1

  with Parallel(n_jobs=n_jobs, verbose=verbose, backend='loky') as parallel:
      # The generator is now passed to the 'parallel' object inside the 'with' block
      results = parallel(
          delayed(pickleable_fn)(*d, **kwargs) for d in tqdm(data, desc=desc)
      )

  return results
