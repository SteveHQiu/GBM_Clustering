# Full pipeline 

#%%
import pickle, os
from typing import Iterable
from collections import Counter


import pandas as pd
from pandas import DataFrame, Series
import matplotlib.pyplot as plt


#%% Constants
ROOT_PATH = R"data/SEER RPD 17 Nov 2021.csv"
ROOT_NAME = os.path.splitext(ROOT_PATH)[0]
COUNTERS_PATH = R"data/counters.dat"

#%%
df = pd.read_csv(ROOT_PATH)

#%% Process df and get counts 
df_gbm = df.loc[df["Histologic Type ICD-O-3"].isin([9440, 9441, 9442, 9445])] # df.isin() method used since "in" operator doesn't work in this context
df_gbm.to_csv(F"{ROOT_PATH}_gbm.csv")
print(df_gbm["Patient ID"].nunique()) 
# SEER RPD 17 Nov 2021 should have ~81,885,000 total encatchment
# Also remember that SEER RPD 17 Nov 2021 is cumulative over 19 years 
gbm_ids = set(df_gbm["Patient ID"].values)


df_gbm2 = df.loc[df["SEER Brain and CNS Recode"] == "1.1.2 Glioblastoma"] # Alternative way to get GBM
print(df_gbm2["Patient ID"].nunique())

df_gbm_pt = df.loc[df["Patient ID"].isin(gbm_ids)]
print(df_gbm_pt["Patient ID"].nunique())

# Isolate non-GBM entries in GBM patients
df_gbm_rel = df_gbm_pt.loc[~df_gbm_pt["Histologic Type ICD-O-3"].isin([9440, 9441, 9442, 9445])] # "~" unary operator to invert
df_gbm_rel.to_csv(F"{ROOT_PATH}_gbm_rel.csv")


all_site_cnt: Counter[str] = Counter(df["Site recode ICD-O-3/WHO 2008"])
all_type_cnt: Counter[str] = Counter(df["ICD-O-3 Hist/behav"])
gbm_site_cnt: Counter[str] = Counter(df_gbm_rel["Site recode ICD-O-3/WHO 2008"])
gbm_type_cnt: Counter[str] = Counter(df_gbm_rel["ICD-O-3 Hist/behav"])

with open(COUNTERS_PATH, "wb+") as file:
    pickle.dump([all_site_cnt,
                 all_type_cnt,
                 gbm_site_cnt,
                 gbm_type_cnt,
                 ], file) # Package counters into single object

#%% Load data checkpoint

with open(COUNTERS_PATH, "rb") as file:
    counters: list[Counter] = pickle.load(file)
    [all_site_cnt,
     all_type_cnt,
     gbm_site_cnt,
     gbm_type_cnt,
     ] = counters # Unpack counters

#%% Visualize raw counts of GBM-related cancers

site_data = list(zip(*gbm_site_cnt.most_common(15))) # list(zip(*var)) to transpose
type_data = list(zip(*gbm_type_cnt.most_common(15))) # list(zip(*var)) to transpose

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

plt.savefig("figures/GBM_assoc.png", bbox_inches="tight")

#%% Visualize incidence-normalized GBM-related cancers
norm_sites = Counter()
for site in gbm_site_cnt:
    norm_count = gbm_site_cnt[site]/all_site_cnt[site]*100 # Calculate % of gbm-related incidences to overall incidence
    norm_sites[site] = norm_count


norm_types = Counter()
for ctype in gbm_type_cnt:
    norm_count = gbm_type_cnt[ctype]/all_type_cnt[ctype]*100 # Calculate % of gbm-related incidences to overall incidence
    norm_types[ctype] = norm_count
    
    
site_data = list(zip(*norm_sites.most_common(15))) # list(zip(*var)) to transpose
type_data = list(zip(*norm_types.most_common(15))) # list(zip(*var)) to transpose

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5))
ax1.bar(site_data[0], site_data[1])
plt.draw() # Calculate labels, otherwise .get_xticklabels() doesn't return anything
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, ha="right")
ax1.set_ylabel("Percentage of tumours at site associated with GBM")
ax1.set_xlabel("Tumour site of associated malignancy")

ax2.bar(type_data[0], type_data[1])
plt.draw() # Calculate labels, otherwise .get_xticklabels() doesn't return anything
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, ha="right")
ax2.set_ylabel("Percentage of tumours with same histology code associated with GBM")
ax2.set_xlabel("Tumour ICD-O-3 histology code of associated malignancy")

plt.savefig("figures/GBM_assoc_norm.png", bbox_inches="tight")

#%% 