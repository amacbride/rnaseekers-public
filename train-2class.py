#!/usr/bin/env python

# run training over teh PDAC datasets, holding back a tumor and normal

import anndata
from joint import joint, utils
import scanpy as sc
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
import numpy as np

# load datasets from disk
normal_base = "AdjNorm_TISSUE_"
# skip sample 3
normals = ['1', '2']

tumor_base = "PDAC_TISSUE_"
# skip sample 2
tumors = ['1', '3', '4', '5', '6', '7', '8', '9', '10', '11A', '11B', '12', '13', '15']

adatas = []

for n in normals:
    filename = normal_base + n + "/" + normal_base + n + "_output.h5ad"
    print(filename)
    na = anndata.read_h5ad(filename)
    na.obs['tissue_type'] = 'normal'
    na.var_names_make_unique()
    adatas.append(na)

for t in tumors:
    filename = tumor_base + n + "/"  + tumor_base + n + "_output.h5ad"
    print(filename)
    ta = anndata.read_h5ad(filename)
    ta.obs['tissue_type'] = 'tumor'
    ta.var_names_make_unique()
    adatas.append(ta)

print("length = ")
print(len(adatas))

print ("Normalize and log-transform each dataset")
# Normalize and log-transform each dataset
for adata in adatas:
    print(adata)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

# Concatenate all datasets
combined = anndata.concat(adatas)

# Extract features (X) and target (y)
X = combined.X
y = combined.obs['tissue_type']

# Encode labels
encoder = LabelEncoder()
y_encoded = encoder.fit_transform(y)

# Split data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.2, random_state=42)

# Define neural network architecture
model = Sequential()
model.add(Dense(128, input_dim=X_train.shape[1], activation='relu')) 
model.add(Dense(64, activation='relu'))
model.add(Dense(32, activation='relu'))
model.add(Dense(1, activation='sigmoid')) # Use sigmoid for binary classification

# Compile the model
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

# Train the model
model.fit(X_train, y_train, epochs=50, batch_size=32, validation_data=(X_test, y_test))

# Evaluate the model
_, accuracy = model.evaluate(X_test, y_test)
print('Accuracy: %.2f' % (accuracy*100))

# write out the model to disk
model.save("trained_model.h5")
