#%%
import pickle
from typing import Iterable
from collections import Counter


import pandas as pd
from pandas import DataFrame, Series
import matplotlib.pyplot as plt

#%%
df = pd.read_excel("data/export2.xlsx")

# %%
df_gbm = df.loc[df["Histologic Type ICD-O-3"].isin([9440, 9441, 9442])] # df.isin() method used since "in" operator doesn't work in this context
df_gbm.to_excel("data/export2_gbm.xlsx")
print(df_gbm["Patient ID"].nunique())
gbm_ids = set(df_gbm["Patient ID"].values)

#%%
df_gbm_pt = df.loc[df["Patient ID"].isin(gbm_ids)]
print(df_gbm_pt["Patient ID"].nunique())

# Isolate non-GBM entries in GBM patients
df_gbm_rel = df_gbm_pt.loc[~df_gbm_pt["Histologic Type ICD-O-3"].isin([9440, 9441, 9442])] # "~" unary operator to invert
df_gbm_rel.to_excel("data/export2_gbm_rel.xlsx")


site_counter = Counter(df_gbm_rel["Site recode ICD-O-3/WHO 2008"])
type_counter = Counter(df_gbm_rel["ICD-O-3 Hist/behav"])

with open("data/counters.dat", "wb+") as file:
    pickle.dump([site_counter, type_counter], file) # Package counters into single object

#%%

site_data = list(zip(*site_counter.most_common(15))) # list(zip(*var)) to transpose
type_data = list(zip(*type_counter.most_common(15))) # list(zip(*var)) to transpose

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5))
ax1.bar(site_data[0], site_data[1])
plt.draw() # Calculate labels, otherwise .get_xticklabels() doesn't return anything
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, ha="right")
ax1.set_ylabel("Number of GBM cases")
ax1.set_xlabel("Tumour site of associated malignancy")

ax2.bar(type_data[0], type_data[1])
plt.draw() # Calculate labels, otherwise .get_xticklabels() doesn't return anything
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, ha="right")
ax2.set_ylabel("Number of GBM cases")
ax2.set_xlabel("Tumour ICD-O-3 histology code of associated malignancy")

plt.savefig("figures/GBM_associations.png", bbox_inches="tight")