import matplotlib.pyplot as plt
import pandas as pd
import qq_plot # This is a in-house code to make QQ plot

# ------------------- Question 1. Plot a histogram of allele frequency -------------------
# Read in AF.txt file
fn_af = 'AF.txt'
lst_alt_allele_frequency = []
with open(fn_af) as fh:
    line = fh.readline() # Skip header line
    line = fh.readline().strip()
    while line != '':
        lst_alt_allele_frequency.append(float(line.split()[1]))
        line = fh.readline().strip()

fig_af, ax_af = plt.subplots(figsize=(10,5), dpi=150)
ax_af.hist(lst_alt_allele_frequency, bins=50, rwidth=0.8)
ax_af.set_title('Distribution of estimated AF of alternative allele')
ax_af.set_xlabel('Allele frequency')
ax_af.set_ylabel('Count')
fig_af.savefig('Q1_AF_histogram.jpeg')

# Set threshold of allele frequency to plot
threshold_af = 0.05
lst_alt_allele_frequency_non_zeros = []
for i in lst_alt_allele_frequency:
    if i>threshold_af: lst_alt_allele_frequency_non_zeros.append(i)
fig_af_non_zeros, ax_af_non_zeros = plt.subplots(figsize=(10,5), dpi=150)
ax_af_non_zeros.hist(lst_alt_allele_frequency_non_zeros, bins=50, rwidth=0.8)
ax_af_non_zeros.set_title('Distribution of estimated AF of alternative allele (AF>'+str(threshold_af)+')')
ax_af_non_zeros.set_xlabel('Allele frequency')
ax_af_non_zeros.set_ylabel('Count')
fig_af_non_zeros.savefig('Q1_AF_gt_'+ str(threshold_af) +'_histogram.jpeg')

# ------------------- Question 2. HWE -------------------
# Plot QQ plot for HWE p values
qq_plot.qqplot(filename='HWE.txt', output='Q2_HWE_p_values_QQ_plot.jpeg',
               p_value_column_title = 'HWE_p_vals', title='HWE p values')

# ------------------- Question 3. LD -------------------







# ------------------- Question 4. PCA -------------------
fn_pca_variance = 'PCA_variance_ratio_explained_by_each_PC.txt'
fn_pc_vals = 'PCs.txt'
lst_variance_explained_ratio = [] # Ratio of Variance explained by each PC

df_PCs = pd.read_csv(fn_pc_vals, sep='\t')
with open(fn_pca_variance) as fh:
    line = fh.readline() # Skip header
    line = fh.readline().strip()
    while line != '':
        lst_variance_explained_ratio.append(float(line.split()[2])*100)
        line = fh.readline().strip()

fig_pca, ax_pca = plt.subplots(nrows=2, ncols=2, dpi=150, figsize=(12, 9))
# Create a scree plot to show variance explained by each the first 3
number_of_PCs = 3
ax_pca[0,0].plot([x for x in range(1, number_of_PCs+1)], lst_variance_explained_ratio[:number_of_PCs], ls='-', marker='o')
ax_pca[0,0].set_title('Variance explained by the first '+str(number_of_PCs)+' PCs')
ax_pca[0,0].set_xlabel('PC number')
ax_pca[0,0].set_ylabel('Variance explained (%)')

# Plot PC1-PC2, PC1-PC3, PC2-PC3
ax_pca[0,1].plot(df_PCs['PC1'], df_PCs['PC2'], ls='', marker='o')
ax_pca[0,1].set_title('PC1 vs. PC2')
ax_pca[0,1].set_xlabel('PC1')
ax_pca[0,1].set_ylabel('PC2')

ax_pca[1,0].plot(df_PCs['PC1'], df_PCs['PC3'], ls='', marker='o')
ax_pca[1,0].set_title('PC1 vs. PC3')
ax_pca[1,0].set_xlabel('PC1')
ax_pca[1,0].set_ylabel('PC3')

ax_pca[1,1].plot(df_PCs['PC2'], df_PCs['PC3'], ls='', marker='o')
ax_pca[1,1].set_title('PC2 vs. PC3')
ax_pca[1,1].set_xlabel('PC2')
ax_pca[1,1].set_ylabel('PC3')

fig_pca.tight_layout()
fig_pca.savefig('PCA.jpeg')

# ------ More plot of PCA result (for my own use) ------
fig_pca_2, ax_pca_2 = plt.subplots(nrows=2, ncols=2, dpi=150, figsize=(12, 9))
# Create a scree plot to show variance explained by each PC
number_of_PCs = 30
ax_pca_2[0,0].plot([x for x in range(1, number_of_PCs+1)], lst_variance_explained_ratio[:number_of_PCs], ls='--', marker='.')
ax_pca_2[0,0].set_title('Variance explained by the first '+str(number_of_PCs)+' PCs')
ax_pca_2[0,0].set_xlabel('PC number')
ax_pca_2[0,0].set_ylabel('Variance explained (%)')

# Plot PC1-PC2, PC1-PC3, PC2-PC3
# Color dots by two clusters in PC1 vs. PC2: PC1=0, PC2 values does not matter
ax_pca_2[0,1].plot(df_PCs[df_PCs['PC1']>=0]['PC1'], df_PCs[df_PCs['PC1']>=0]['PC2'], ls='', marker='o', label='Cluster 1')
ax_pca_2[0,1].plot(df_PCs[df_PCs['PC1']<0]['PC1'], df_PCs[df_PCs['PC1']<0]['PC2'], ls='', marker='o', label='Cluster 2')
ax_pca_2[0,1].set_title('PC1 vs. PC2')
ax_pca_2[0,1].set_xlabel('PC1')
ax_pca_2[0,1].set_ylabel('PC2')
ax_pca_2[0,1].legend()

ax_pca_2[1,0].plot(df_PCs[df_PCs['PC1']>=0]['PC1'], df_PCs[df_PCs['PC1']>=0]['PC3'], ls='', marker='o', label='Cluster 1')
ax_pca_2[1,0].plot(df_PCs[df_PCs['PC1']<0]['PC1'], df_PCs[df_PCs['PC1']<0]['PC3'], ls='', marker='o', label='Cluster 2')
ax_pca_2[1,0].set_title('PC1 vs. PC3')
ax_pca_2[1,0].set_xlabel('PC1')
ax_pca_2[1,0].set_ylabel('PC3')
ax_pca_2[1,0].legend()

ax_pca_2[1,1].plot(df_PCs[df_PCs['PC1']>=0]['PC2'], df_PCs[df_PCs['PC1']>=0]['PC3'], ls='', marker='o', label='Cluster 1')
ax_pca_2[1,1].plot(df_PCs[df_PCs['PC1']<0]['PC2'], df_PCs[df_PCs['PC1']<0]['PC3'], ls='', marker='o', label='Cluster 2')
ax_pca_2[1,1].set_title('PC2 vs. PC3')
ax_pca_2[1,1].set_xlabel('PC2')
ax_pca_2[1,1].set_ylabel('PC3')
ax_pca_2[1,1].legend()

fig_pca_2.tight_layout()
fig_pca_2.savefig('PCA_v2.jpeg')