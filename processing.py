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
COL_RAD = "Radiation recode"
COL_ID = "Patient ID"
COL_SEQ = "Record number recode" # Sequentially numbers a person's tumors within each SEER submission. Order is based on date of diagnosis and then sequence #. 
COL_FIRST_PRIM = "First malignant primary indicator" # First MALIGNANT cancer; does not necessarily mean first cancer as could have many prior non-malignant neoplasm
COL_PRIMARY = "Primary by international rules" # Whether its primary or mets
COL_ORD_PRIM = "Sequence number" # Order of the primary AT THE TIME OF DIAGNOSIS, should only have 1-2 primaries when first diagnosed (rarely should you be walking in with 3+ primaries without being diagnosed) FIXME UNSURE HOW IT CHANGES WITH EACH ENTRY, SHOULD DOUBLE CHECK FOR A SPECIFIC PATIENT WITH MULTIPLE ENTRIES
GBM_HIST_CODES = [9440, 9441, 9442, 9445]

#%% Load CSV
df = pd.read_csv(ROOT_PATH)
df[COL_SURV] = pd.to_numeric(df[COL_SURV], errors="coerce") # Coerce converts non-numerics into NaN
#%%
if 1: # Visualize variable space of df
    for col in df.columns:
        print(col)
        print(df[col].unique())
#%% Get different rates of GBM 

def tempFx(df: DataFrame, frac_gen = 0.000556793063442633):
    df_gbm = df.loc[df[COL_HIST].isin(GBM_HIST_CODES)] # df.isin() method used since "in" operator doesn't work in this context
    n_total = df[COL_ID].nunique()
    n_gbm = df_gbm[COL_ID].nunique()
    n_encatchment = 81885000 # SEER RPD 17 Nov 2021 should have ~81,885,000 total encatchment
    n_years = 20 # Also remember that SEER RPD 17 Nov 2021 is cumulative over 20 years 
    print(F'Number of patients: {n_gbm}') 
    print(F'Fraction of encatchment with GBM: {n_gbm / n_encatchment}\nn = {n_gbm}, N = {n_encatchment}') 
    print(F'GBM rate per 100 000 per year: {n_gbm / (n_encatchment / 100000) / n_years}') # Should be roughly 3.19

    df_gbm_first = df_gbm[df_gbm[COL_ORD_PRIM].str.contains('One primary only|1st of 2 or more primaries')]
    n_gbm_first = df_gbm_first[COL_ID].nunique() 
    assert n_gbm_first == len(df_gbm_first) # Each entry should be unique

    if 0: 
        for pt_id in df_gbm_first[df_gbm_first[COL_SEQ] != 1][COL_ID]: # Output number of entries for patients where their first GBM primary was not the first SEER entry, FIXME Not sure why there should be any of these in the first place
            print(len(df[df[COL_ID] == pt_id]))

    df_rem = df[COL_ID].isin(df_gbm_first[COL_ID])
    
    n_rem_entries = df_rem.value_counts()[True] # Number of entries to remove
    print(F"Entries to remove: {n_rem_entries} | Excess: {n_rem_entries - len(df_gbm_first)}")
    df_no_gbm_first = df.drop(df[df_rem].index)
    print(F"Number of entries removed: {len(df) - len(df_no_gbm_first)}")
    n_no_gbm_first = df_no_gbm_first[COL_ID].nunique() 

    df_gbm_sec = df_no_gbm_first.loc[df_no_gbm_first[COL_HIST].isin(GBM_HIST_CODES)]
    n_gbm_sec = df_gbm_sec[COL_ID].nunique()

    if not frac_gen: # If general rate of GBM not given, then define it here
        frac_gen = (n_gbm - n_gbm_sec) / n_encatchment
    frac_sec = n_gbm_sec / n_no_gbm_first

    print(F'Fraction of encatchment with primary GBM: {frac_gen}\nn = {(n_gbm - n_gbm_sec)}, N = {n_encatchment}') 
    print(F"Fraction of secondary GBM: {frac_sec}\nn = {n_gbm_sec}, N = {n_no_gbm_first}")
    print(F">>> Ratio: {frac_sec/frac_gen} <<<")

#%%

for col in df["Sex"].unique():
    df_temp = df[df["Sex"] == col]
    print(F"Sex: {col} ================================")
    tempFx(df_temp)
    
for col in df["Age recode with <1 year olds"].unique():
    df_temp = df[df["Age recode with <1 year olds"] == col]
    print(F"Age: {col} ================================")
    tempFx(df_temp)
#%%
df_temp = df[df["Age recode with <1 year olds"].isin(['75-79 years', '80-84 years', '85+ years'])]
tempFx(df_temp)
#%%
df_temp = df[(df["Sex"] == "Male") & (df["Age recode with <1 year olds"].isin(['75-79 years', '80-84 years', '85+ years']))]
tempFx(df_temp)
#%%
df_temp = df[(df["Sex"] == "Male") & (df["Age recode with <1 year olds"].isin(['75-79 years']))]
tempFx(df_temp)
#%%



#%%
def tempFx2(df: DataFrame, frac_gen = 0.000556793063442633):
    
    'C61.9-Prostate gland'
    
    df_gbm = df.loc[df[COL_HIST].isin(GBM_HIST_CODES)] # df.isin() method used since "in" operator doesn't work in this context

    # df_gbm_first = df_gbm[df_gbm[COL_ORD_PRIM].str.contains('One primary only|1st of 2 or more primaries')]
    # n_gbm_first = df_gbm_first[COL_ID].nunique() 
    # assert n_gbm_first == len(df_gbm_first) # Each entry should be unique

    # df_rem = df[COL_ID].isin(df_gbm_first[COL_ID])
    
    # n_rem_entries = df_rem.value_counts()[True] # Number of entries to remove
    # print(F"Entries to remove: {n_rem_entries} | Excess: {n_rem_entries - len(df_gbm_first)}")
    # df_no_gbm_first = df.drop(df[df_rem].index)
    # print(F"Number of entries removed: {len(df) - len(df_no_gbm_first)}")
    # n_no_gbm_first = df_no_gbm_first[COL_ID].nunique() 
    n_cases = df[COL_ID].nunique()

    n_gbm_sec = df_gbm[COL_ID].nunique()

    frac_sec = n_gbm_sec / n_cases

    print(F"Fraction of secondary GBM: {frac_sec}\nn = {n_gbm_sec}, N = {n_cases}")
    print(F">>> Ratio: {frac_sec/frac_gen} <<<")


df_pros_entries = df[df["Primary Site - labeled"] == 'C61.9-Prostate gland']
df_pros_entries2 = df_pros_entries[df_pros_entries[COL_ORD_PRIM].str.contains('One primary only|1st of 2 or more primaries')]

df_prost = df[df[COL_ID].isin(df_pros_entries2[COL_ID])]

df_temp = df_prost[(df_prost["Sex"] == "Male") & (df_prost["Age recode with <1 year olds"].isin(['75-79 years']))]
tempFx2(df_temp)


#%% Count sites and types for entries 

df_gbm2 = df.loc[df["SEER Brain and CNS Recode"] == "1.1.2 Glioblastoma"] # Alternative way to get GBM
print(df_gbm2["Patient ID"].nunique())

gbm_ids = set(df_gbm["Patient ID"].values)
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

# Can probably compare rates using chi-squared https://www.medcalc.org/calc/rate_comparison.php