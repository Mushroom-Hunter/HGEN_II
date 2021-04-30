import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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
    if (i>threshold_af and i<=0.5) or (1-i<threshold_af and i>0.5):
        lst_alt_allele_frequency_non_zeros.append(i)

fig_af_non_zeros, ax_af_non_zeros = plt.subplots(figsize=(7, 4), dpi=100)
ax_af_non_zeros.hist(lst_alt_allele_frequency_non_zeros, bins=15, rwidth=0.8)
ax_af_non_zeros.set_title('Distribution of estimated AF of alternative allele (MAF>'+str(threshold_af)+')')
ax_af_non_zeros.set_xlabel('Allele frequency')
ax_af_non_zeros.set_ylabel('Count')
fig_af_non_zeros.savefig('Q1_AF_gt_'+ str(threshold_af) +'_histogram.jpeg')

# ------------------- Question 2. HWE -------------------
# Plot QQ plot for HWE p values
qq_plot.qqplot(filename='HWE.txt', output='Q2_HWE_p_values_QQ_plot.jpeg',
               p_value_column_title = 'HWE_p_vals', title='HWE p values')
# Plot HWE p value distribution
df_HWE_p = pd.read_csv('HWE.txt', sep='\t')
fig_HWE_p, ax_HWE_p = plt.subplots(figsize=(8,6), dpi=150)
ax_HWE_p.hist(df_HWE_p['HWE_p_vals'], bins=50, rwidth=0.8)
ax_HWE_p.set_xlabel('p values')
ax_HWE_p.set_ylabel('Counts')
ax_HWE_p.set_title('Distribution of HWE p values')
fig_HWE_p.savefig('Q2_HWE_p_values_distribution.jpeg')

# ------------------- Question 3. LD -------------------
# Create a heatmap from LD scores
df_ld_d = pd.read_csv('LD_D.txt', sep='\t')
df_ld_d_prime = pd.read_csv('LD_Dprime.txt', sep='\t')
df_ld_r2 = pd.read_csv('LD_r2.txt', sep='\t')

# Merge 3 types of LD scores, sort by r2 then plot to compare
df_ld_merged = pd.DataFrame(columns=['SNP1', 'SNP2','LD_D', 'LD_Dprime', 'LD_r2'])
df_ld_merged['LD_D'] = abs(df_ld_d['LD_D'])
df_ld_merged['LD_Dprime'] = abs(df_ld_d_prime['LD_Dprime'])
df_ld_merged['LD_r2'] = df_ld_r2['LD_r2']
df_ld_merged['SNP1'] = df_ld_r2['SNP1']
df_ld_merged['SNP2'] = df_ld_r2['SNP2']

df_ld_merged.sort_values(by='LD_r2', inplace=True)

# Get a list of all SNPs from LD score data
mask = df_ld_d['SNP1']==df_ld_d.iloc[0,0]
all_snps = df_ld_d[mask]['SNP2'].values

# Construct 3 matrices using above dataframes, rows and columns are SNPs, values are LD scores
# Need to specify dtype, other wise pd assume values are object type in this way of construction
df_ld_d_matrix = pd.DataFrame(columns=all_snps, index=all_snps, dtype='float')
df_ld_d_prime_matrix = pd.DataFrame(columns=all_snps, index=all_snps, dtype='float')
df_ld_r2_matrix = pd.DataFrame(columns=all_snps, index=all_snps, dtype='float')

total_num_of_snps = len(all_snps)
for snp_index in range(len(all_snps)):
    mask_snp = df_ld_d['SNP1'] == all_snps[snp_index]
    # Only fill some columns base on available number of values
    # This is a triangle not really a matrix (symmetrical values are ignored)
    num_of_vals = len(df_ld_d[mask_snp]['LD_D'])
    # Use absolute values for D
    df_ld_d_matrix.iloc[total_num_of_snps-num_of_vals:, snp_index] = abs(df_ld_d[mask_snp]['LD_D'].values)
    # Use absolute values for D' since sign does not matter
    df_ld_d_prime_matrix.iloc[total_num_of_snps-num_of_vals:, snp_index] = abs(df_ld_d_prime[mask_snp]['LD_Dprime'].values)
    df_ld_r2_matrix.iloc[total_num_of_snps-num_of_vals:, snp_index] = df_ld_r2[mask_snp]['LD_r2'].values

fig_ld, ax_ld = plt.subplots(nrows=2, ncols=2, figsize=(8, 8), dpi=200)
# Compare D, D' and r2
ax_ld[0, 0].plot([x for x in range(df_ld_merged.shape[0])], df_ld_merged['LD_D'], ls='', marker='.', label='D', alpha=0.6)
ax_ld[0, 0].plot([x for x in range(df_ld_merged.shape[0])], df_ld_merged['LD_Dprime'], ls='', marker='.', alpha=0.6, label='D\'')
ax_ld[0, 0].plot([x for x in range(df_ld_merged.shape[0])], df_ld_merged['LD_r2'], ls='', marker='.', alpha=0.4, label=r'$R^2$')
ax_ld[0, 0].set_title(r"D, D' and $R^2$")
ax_ld[0, 0].legend(bbox_to_anchor=(0.5, -0.05), ncol=3, loc='upper center')

# Plot D
sns.heatmap(df_ld_d_matrix.values, ax=ax_ld[0, 1], cmap='YlOrRd', square=True,
            xticklabels=10, yticklabels=10, cbar_kws={"shrink": 0.8})
ax_ld[0, 1].set_title('D')
ax_ld[0, 1].set_xlabel('SNP1')
ax_ld[0, 1].set_ylabel('SNP2')
# Plot D'
sns.heatmap(df_ld_d_prime_matrix.values, ax=ax_ld[1, 0], cmap='YlOrRd', square=True,
            xticklabels=10, yticklabels=10, cbar_kws={"shrink": 0.8}) # use cbar_kws={"shrink": .8} to change colorbar length
ax_ld[1, 0].set_title('D\'')
ax_ld[1, 0].set_xlabel('SNP1')
ax_ld[1, 0].set_ylabel('SNP2')
# Plot R2
sns.heatmap(df_ld_r2_matrix.values, ax=ax_ld[1, 1], cmap='YlOrRd', square=True,
            xticklabels=10, yticklabels=10, cbar_kws={"shrink": 0.8})
ax_ld[1, 1].set_title(r'$R^2$')
ax_ld[1, 1].set_xlabel('SNP1')
ax_ld[1, 1].set_ylabel('SNP2')

fig_ld.suptitle('LD scores')
fig_ld.tight_layout()
fig_ld.savefig('Q3_LD.jpeg')

# ---- More plot to compare allele frequencies with D, D' and r2 ----
df_af = pd.read_csv(fn_af, sep='\t')
# Put LD scores and allele frequencies distribution
fig_ld_v2, ax_ld_v2 = plt.subplots(nrows=2, figsize=(6, 8), dpi=100)
# Compare D, D' and r2
ax_ld_v2[0].plot([x for x in range(df_ld_merged.shape[0])], df_ld_merged['LD_D'], ls='', marker='.', label='D', alpha=0.6)
ax_ld_v2[0].plot([x for x in range(df_ld_merged.shape[0])], df_ld_merged['LD_Dprime'], ls='', marker='.', alpha=0.6, label='D\'')
ax_ld_v2[0].plot([x for x in range(df_ld_merged.shape[0])], df_ld_merged['LD_r2'], ls='', marker='.', alpha=0.4, label=r'$R^2$')
ax_ld_v2[0].set_title(r"D, D' and $R^2$")
ax_ld_v2[0].legend(bbox_to_anchor=(0.5, -0.05), ncol=3, loc='upper center')

ax_ld_v2[1].hist(lst_alt_allele_frequency_non_zeros, bins=15, rwidth=0.8)
ax_ld_v2[1].set_title('Distribution of estimated AF of alternative allele (MAF>'+str(threshold_af)+')')
ax_ld_v2[1].set_xlabel('Allele frequency')
ax_ld_v2[1].set_ylabel('Counts')

fig_ld_v2.tight_layout()
fig_ld_v2.savefig('Q3_LD_values_and_AF.jpeg')

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
fig_pca.savefig('Q4_PCA.jpeg')

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
fig_pca_2.savefig('Q4_PCA_colored.jpeg')