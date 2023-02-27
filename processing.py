# Full pipeline 

#%%
import pickle, os
from typing import Iterable
from collections import Counter


import pandas as pd
from pandas import DataFrame, Series
import matplotlib.pyplot as plt

from internals import LOG


#%% Constants
ROOT_PATH = R"C:\Users\steve\Downloads\ZLocal\SEER RPD 17 Nov 2021.csv"
ROOT_PATH = R"data/SEER RPD 17 Nov 2021.csv"
ROOT_NAME = os.path.splitext(ROOT_PATH)[0]
COUNTERS_PATH = R"data/counters.dat"
COL_AGE = "Age recode with <1 year olds"
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

EXPORT_CHECKPOINTS = 0

#%% Load CSV
df = pd.read_csv(ROOT_PATH)
df[COL_SURV] = pd.to_numeric(df[COL_SURV], errors="coerce") # Coerce converts non-numerics into NaN
n_encatchment = 81885000 # SEER RPD 17 Nov 2021 should have ~81,885,000 total encatchment
n_years = 20 # Also remember that SEER RPD 17 Nov 2021 is cumulative over 20 years 
n_total = df[COL_ID].nunique()

if 0: # Visualize variable space of df
    for col in df.columns:
        LOG.info(col)
        LOG.info(df[col].unique())
#%% Get GBM records

df_gbm = df.loc[df[COL_HIST].isin(GBM_HIST_CODES)] # df.isin() method used since "in" operator doesn't work in this context
n_gbm = df_gbm[COL_ID].nunique()

LOG.info(F'Number of patients: {n_gbm}') 
LOG.info(F'Fraction of encatchment with GBM: {n_gbm / n_encatchment}\nn = {n_gbm}, N = {n_encatchment}') 
LOG.info(F'GBM rate per 100 000 per year: {n_gbm / (n_encatchment / 100000) / n_years}') # Should be roughly 3.19

if EXPORT_CHECKPOINTS:
    df_gbm.to_excel(F"{ROOT_NAME}_gbm.xlsx")


#%%

def reportCancerIncidence(df: DataFrame,
           hist_codes = GBM_HIST_CODES,
           prefix = "gbm",
           incidence_20y_first = 0.000556793,
           export_dfs = False,
           ):
    # Get primary and secondary rates for a target cancer in a DataFrame
    # Target cancer by histology codes
    # frac_gen is prevalence of target cancer as first cancer (all occurences minus secondary ocurrences)
    # frac_gen calculated by (n_gbm - n_gbm_sec) / n_encatchment
    
    df_target_cancer = df.loc[df[COL_HIST].isin(hist_codes)] # df.isin() method used since "in" operator doesn't work in this context
    
    # Find first primaries of target cancer and their patients' related entries
    df_first = df_target_cancer[df_target_cancer[COL_ORD_PRIM].str.contains('One primary only|1st of 2 or more primaries')]
    n_first = df_first[COL_ID].nunique() 
    if export_dfs:
        df_first.to_excel(F"{ROOT_NAME}_{prefix}_first.xlsx")
    
    df_gbm_first_rel = df[df[COL_ID].isin(df_first[COL_ID])]
    if export_dfs:
        df_gbm_first_rel.to_excel(F"{ROOT_NAME}_{prefix}_first_rel.xlsx")
    
    # Debug
    if not n_first == len(df_first): # Each entry in first entry should be unique
        LOG.warning(F"Not all first tumour entries are unique\nLength: {len(df_first)} | Unique: {n_first}")

    if 0: 
        for pt_id in df_first[df_first[COL_SEQ] != 1][COL_ID]: # Output number of entries for patients where their first GBM primary was not the first SEER entry, FIXME Not sure why there should be any of these in the first place
            LOG.info(len(df[df[COL_ID] == pt_id]))
    
    # Define entries to be removed to get primaries of target cancer that were not the first primaries 
    ids_entr_rel_to_first = df[COL_ID].isin(df_first[COL_ID])
    df_not_first = df.drop(df[ids_entr_rel_to_first].index)
    LOG.info(F"Entries to remove: {ids_entr_rel_to_first.value_counts()[True]}")
    LOG.info(F"Number of entries removed: {len(df) - len(df_not_first)}")
    n_not_first = df_not_first[COL_ID].nunique()

    # Find second+ primaries of target cancer and their patients' related entries
    df_second = df_not_first.loc[df_not_first[COL_HIST].isin(hist_codes)]
    n_gbm_sec = df_second[COL_ID].nunique()
    if export_dfs:
        df_second.to_excel(F"{ROOT_NAME}_{prefix}_second.xlsx")
    
    df_gbm_sec_rel = df[df[COL_ID].isin(df_second[COL_ID])]
    if export_dfs:
        df_gbm_sec_rel.to_excel(F"{ROOT_NAME}_{prefix}_sec_rel.xlsx")

    # Calculate rates
    if not incidence_20y_first: # If general rate of GBM not given, then define it here
        incidence_20y_first = (n_gbm - n_gbm_sec) / n_encatchment
        LOG.info(F'Fraction of encatchment with primary GBM: {incidence_20y_first}\nn = {(n_gbm - n_gbm_sec)}, N = {n_encatchment}') 
    else:    
        LOG.info(F'Imported 20 year incidence: {incidence_20y_first}') # Need % of encatchment fitting subpopulation of df to calculate 
    
    incidence_20y_sec = n_gbm_sec / n_not_first # Rate of target cancer as non-first cancer within those with cancer over the past 20 years
    ratio = incidence_20y_sec/incidence_20y_first
    
    LOG.info(F"Fraction of secondary GBM: {incidence_20y_sec}\nn = {n_gbm_sec}, N = {n_not_first}")
    LOG.info(F">>> Ratio: {ratio} <<<")
    
    return ratio
#%%
reportCancerIncidence(df, incidence_20y_first=None, export_dfs=True)

#%% Age

cancer_sites = list(df[COL_AGE].unique())
cancer_sites.sort()
site_ratios = []
for col in cancer_sites:
    df_temp = df[df[COL_AGE] == col]
    LOG.info(F"Age: {col} ================================")
    site_ratios.append(reportCancerIncidence(df_temp)) # Add ratio 

fig, ax = plt.subplots()
ax.bar(cancer_sites, site_ratios)
ax.set_ylabel("Ratio of non-first GBM occurence rate")
fig.autofmt_xdate(rotation=45)
    
#%% Site

#%% 

for col in df["Sex"].unique():
    df_temp = df[df["Sex"] == col]
    LOG.info(F"Sex: {col} ================================")
    reportCancerIncidence(df_temp)


#%%
df_temp = df[df["Age recode with <1 year olds"].isin(['75-79 years', '80-84 years', '85+ years'])]
reportCancerIncidence(df_temp)
#%%
df_temp = df[(df["Sex"] == "Male") & (df["Age recode with <1 year olds"].isin(['75-79 years', '80-84 years', '85+ years']))]
reportCancerIncidence(df_temp)
#%%
df_temp = df[(df["Sex"] == "Male") & (df["Age recode with <1 year olds"].isin(['75-79 years']))]
reportCancerIncidence(df_temp)
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
    # LOG.info(F"Entries to remove: {n_rem_entries} | Excess: {n_rem_entries - len(df_gbm_first)}")
    # df_no_gbm_first = df.drop(df[df_rem].index)
    # LOG.info(F"Number of entries removed: {len(df) - len(df_no_gbm_first)}")
    # n_no_gbm_first = df_no_gbm_first[COL_ID].nunique() 
    n_cases = df[COL_ID].nunique()

    n_gbm_sec = df_gbm[COL_ID].nunique()

    frac_sec = n_gbm_sec / n_cases

    LOG.info(F"Fraction of secondary GBM: {frac_sec}\nn = {n_gbm_sec}, N = {n_cases}")
    LOG.info(F">>> Ratio: {frac_sec/frac_gen} <<<")


df_pros_entries = df[df["Primary Site - labeled"] == 'C61.9-Prostate gland']
df_pros_entries2 = df_pros_entries[df_pros_entries[COL_ORD_PRIM].str.contains('One primary only|1st of 2 or more primaries')]

df_prost = df[df[COL_ID].isin(df_pros_entries2[COL_ID])]

df_temp = df_prost[(df_prost["Sex"] == "Male") & (df_prost["Age recode with <1 year olds"].isin(['75-79 years']))]
tempFx2(df_temp)

# ====================== Original pipeline 
#%% Count sites and types for entries 
df_gbm = df.loc[df[COL_HIST].isin(GBM_HIST_CODES)] # df.isin() method used since "in" operator doesn't work in this context
df_gbm = df_gbm.loc[df_gbm["Radiation recode"].isin(['None/Unknown', 'Refused (1988+)'])]
df_gbm = df_gbm.loc[df_gbm["Chemotherapy recode (yes, no/unk)"].isin(['No/Unknown'])]

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
# %%
