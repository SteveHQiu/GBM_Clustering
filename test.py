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

plt.savefig("figures/GBM_associations.png", bbox_inches="tight")