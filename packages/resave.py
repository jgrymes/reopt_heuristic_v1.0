import cloudpickle as pickle
import pandas as pd
import gzip
# with gzip.open(filename, 'rb') as f:
with gzip.open("data/08f599_1.pkl", "rb") as f:
    obj = pickle.load(f,fix_imports=False)

# Re-save using current pandas version
with gzip.open("data/08f599_2.pkl", "wb") as f:
    pickle.dump(obj, f)