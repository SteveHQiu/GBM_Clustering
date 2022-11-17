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
COL_SITE = "Site recode ICD-O-3/WHO 2008"
COL_HIST = "Histologic Type ICD-O-3"
COL_TYPE = "ICD-O-3 Hist/behav"
COL_SURV = "Survival months"

#%%
df = pd.read_csv(ROOT_PATH)
df[COL_SURV] = pd.to_numeric(df[COL_SURV], errors="coerce") # Coerce converts non-numerics into NaN


#%% Process df and get counts 
df_gbm = df.loc[df[COL_HIST].isin([9440, 9441, 9442, 9445])] # df.isin() method used since "in" operator doesn't work in this context
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
df_gbm_rel = df_gbm_pt.loc[~df_gbm_pt[COL_HIST].isin([9440, 9441, 9442, 9445])] # "~" unary operator to invert
df_gbm_rel.to_csv(F"{ROOT_PATH}_gbm_rel.csv")


all_site_cnt: Counter[str] = Counter(df[COL_SITE])
all_type_cnt: Counter[str] = Counter(df[COL_TYPE])
gbm_site_cnt: Counter[str] = Counter(df_gbm_rel[COL_SITE])
gbm_type_cnt: Counter[str] = Counter(df_gbm_rel[COL_TYPE])

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

#%% First cases
df_first = df.loc[df["Record number recode"] == 1] # 6,525,399 first cases 
df_missed = df_first.loc[df_first["Sequence number"].str.match(R"(One)|(^1st)")] # To match "One primary only" and "1st of 2"
print(F"{len(df_first)} first cases, {len(df_missed)} was the first known primary")
print(F"{100-100*len(df_missed)/len(df_first)}% Missing prior cases")
# ~11.3% of first cases was missing a prior cancer in SEER dataset
print(df_missed["Sequence number"].unique())

#%% Cumulative risk 

print(F'{df[COL_SURV].isnull().sum()/len(df)*100}% of entries missing survival data')

site_groups = df.groupby([COL_SITE])[COL_SURV].sum().sort_values(ascending=False)
type_groups = df.groupby([COL_TYPE])[COL_SURV].sum().sort_values(ascending=False)
# Note that each entry has its own survival calculated from diagnosis date to death or current cutoff, independent of co-occurring cancers
# Hence longest survival is always the first entry 
print(site_groups)
print(type_groups)

f_site_groups = df_first.groupby([COL_SITE])[COL_SURV].sum().sort_values(ascending=False)
f_type_groups = df_first.groupby([COL_TYPE])[COL_SURV].sum().sort_values(ascending=False)
# Get survival of first tumours only 

#%% Normalize by cumulative risk 

norm_sites = Counter()
for site in gbm_site_cnt:
    norm_count = gbm_site_cnt[site]/site_groups[site] # Divides by total months of survival across entire cohort
    norm_sites[site] = norm_count


norm_types = Counter()
for ctype in gbm_type_cnt:
    norm_count = gbm_type_cnt[ctype]/type_groups[ctype]*100 # Divides by total months of survival across entire cohort
    norm_types[ctype] = norm_count
    
    
site_data = list(zip(*norm_sites.most_common(15))) # list(zip(*var)) to transpose
type_data = list(zip(*norm_types.most_common(15))) # list(zip(*var)) to transpose

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5))
ax1.bar(site_data[0], site_data[1])
plt.draw() # Calculate labels, otherwise .get_xticklabels() doesn't return anything
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, ha="right")
ax1.set_ylabel("Incidence of tumours at site associated with GBM normalized by cumulative survival of associated tumour")
ax1.set_xlabel("Tumour site of associated malignancy")

ax2.bar(type_data[0], type_data[1])
plt.draw() # Calculate labels, otherwise .get_xticklabels() doesn't return anything
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, ha="right")
ax2.set_ylabel("Incidence of tumour histology associated with GBM normalized by cumulative survival of associated tumour")
ax2.set_xlabel("Tumour ICD-O-3 histology code of associated malignancy")

plt.savefig("figures/GBM_assoc_norm_cum.png", bbox_inches="tight")

#%% Normalize by cumulative risk only based on initial cancer 

norm_sites = Counter()
for site in gbm_site_cnt:
    norm_count = gbm_site_cnt[site]/f_site_groups[site] # Divides by total months of survival across entire cohort
    norm_sites[site] = norm_count


norm_types = Counter()
for ctype in gbm_type_cnt:
    norm_count = gbm_type_cnt[ctype]/f_type_groups[ctype]*100 # Divides by total months of survival across entire cohort
    norm_types[ctype] = norm_count
    
    
site_data = list(zip(*norm_sites.most_common(15))) # list(zip(*var)) to transpose
type_data = list(zip(*norm_types.most_common(15))) # list(zip(*var)) to transpose

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5))
ax1.bar(site_data[0], site_data[1])
plt.draw() # Calculate labels, otherwise .get_xticklabels() doesn't return anything
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, ha="right")
ax1.set_ylabel("Incidence of tumours at site associated with GBM normalized by cumulative survival of associated tumour")
ax1.set_xlabel("Tumour site of associated malignancy")

ax2.bar(type_data[0], type_data[1])
plt.draw() # Calculate labels, otherwise .get_xticklabels() doesn't return anything
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, ha="right")
ax2.set_ylabel("Incidence of tumour histology associated with GBM normalized by cumulative survival of associated tumour")
ax2.set_xlabel("Tumour ICD-O-3 histology code of associated malignancy")

plt.savefig("figures/GBM_assoc_norm_cum_first.png", bbox_inches="tight")

#%% Percentage of GBM cases that had subsequent cancer 

print(df_gbm.groupby(["Sequence number"])["Patient ID"].count())

# Group by patient ID, count cases
# For IDs with more than one entry, get max "Record number recode"
# Check if they had GBM, compare GBM "Record number recode" with max