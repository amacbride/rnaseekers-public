#!/usr/bin/env python

# run inference over the data

import sys
import anndata
import scanpy as sc
from joint import joint, utils

import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import load_model

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

import numpy as np
from collections import Counter
import os

import warnings
from sklearn.exceptions import DataConversionWarning

# Suppress the warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="tensorflow")
warnings.filterwarnings("ignore", category=FutureWarning, message=".*disable_resource_variables.*")

# Suppress specific warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Variable names are not unique")
warnings.filterwarnings("ignore", category=UserWarning, message="`Model.state_updates` will be removed")
warnings.filterwarnings("ignore", category=DataConversionWarning, message="A column-vector y was passed")

# 0 = all messages are logged (default behavior)
# 1 = INFO messages are not printed
# 2 = INFO and WARNING messages are not printed
# 3 = INFO, WARNING, and ERROR messages are not printed

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
tf.get_logger().setLevel('ERROR')

filename = sys.argv[1]
print(filename)
adata = anndata.read_h5ad(filename)
adata.var_names_make_unique()

# Normalize and log-transform each dataset
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

model = load_model("trained_model.h5")

# Extract features (X)
X = adata.X

# Run inference
y_pred = model.predict(X)

# Convert probabilities to class labels
y_pred_labels = (y_pred > 0.5).astype(int)

# Initialize the LabelEncoder
encoder = LabelEncoder()
encoder.classes_ = np.array(['normal', 'tumor']) # Assumes 'normal' is 0 and 'tumor' is 1

# Convert encoded labels back to original labels
y_pred_labels = encoder.inverse_transform(y_pred_labels)

# Print the classifications
#for i, label in enumerate(y_pred_labels):
#    print(f"Sample {i}: {label}")

label_counts = Counter(y_pred_labels)

#for label, count in label_counts.items():
# print(f"Label '{label}' appears {count} times.")

lc = len(y_pred_labels)
tc = label_counts['tumor']
nc = label_counts['normal']

npct = nc/lc*100
tpct = tc/lc*100

# calculate a rough "tumor score"
print(f"score: {npct:.2f} % normal, {tpct:.2f} % tumor")
