# This is the final project of 2021 HGEN II
# Author: Wanying Zhu
# Email: wanying.zhu.1@vanderbilt.edu

import multiprocessing  # For multiprocessing purpose
import pandas as pd
from sklearn.decomposition import PCA
import time # To track how long the script takes to run
from scipy.stats import chisquare
# Can also use this one to calculate chi squared p value
# from scipy.stats import chi2
start_time = time.time()
verbose = True  # For debugging, print out variables
output_log = False  # If true, write some log into output files
save_pca_result = True  # Save PCA results for my own plotting, but not required in this assignment.
# fn = 'HGEN8341_finalproject_data.txt'
# df = pd.read_csv(fn, sep='\t', dtype='str')
fn = 'HGEN8341_finalproject_data.txt.gz'
df = pd.read_csv(fn, sep='\t', compression='gzip', dtype='str')

if verbose:
    print('\n# Starting HGEN II final project by Wanying Zhu\n')
    print('- Original dataset shape:', df.shape)
    print(df.head())

# ------------------ Q1. Allele frequency estimates ------------------
# In the original data file, 0 represents reference allele, other numbers (1, 2 or more) represent alternative allele

# Step 1. Separate 0/1 (alleles) to 0 and 1
# Create a new dataframe (df_allele_separated) with 2 column for each allele, eg. rs116742944_a1 rs116742944_a2
# then separate genotype like 0/1 into the two columns
# Original datafarme df will be needed later, so a new dataframe is necessary
if verbose:
    print('\n\n=======================================\nQ1. Allele frequency estimates')

lst_separated_column_names = ['ID', 'status']  # Create column names for new dataframe
for col in df.columns[2:]:  # The first two items are ID and status, already included when creating the list
    lst_separated_column_names.append(col + '_a1')
    lst_separated_column_names.append(col + '_a2')

df_allele_separated = pd.DataFrame(columns=lst_separated_column_names)
# Fill ID and status columns
df_allele_separated[lst_separated_column_names[0]] = df[lst_separated_column_names[0]]
df_allele_separated[lst_separated_column_names[1]] = df[lst_separated_column_names[1]]
for i in range(2, len(lst_separated_column_names), 2):  # Skip ID and status columns
    df_allele_separated[[lst_separated_column_names[i], lst_separated_column_names[i + 1]]] = df[
        lst_separated_column_names[i] \
            .split('_a1')[0]].str.split('/', expand=True)
if verbose:
    print('- Separate alleles into two columns, new dataframe shape is:', df_allele_separated.shape)
    print(df_allele_separated.head())
    print('- Calculate allele frequency of alternative allele')

# Step 2. Calculate alternative allele (ie note '0' alleles) frequency
mask_alt_allele = df_allele_separated != '0'
alternative_allele_counts = mask_alt_allele.iloc[:, 2:].apply(sum)
number_of_individuals = df.shape[0]
lst_alt_allele_frequency = []  # To store calculated allele frequencies
for i in range(0, len(alternative_allele_counts), 2):
    # Add up counts of two alleles together to get total counts
    total_counts = alternative_allele_counts[i] + alternative_allele_counts[i + 1]
    lst_alt_allele_frequency.append(total_counts / (number_of_individuals * 2))

# Step 3. Output allele frequencies into 'AF.txt' file
# Plotting is in a separated file plotting.py
with open('AF.txt', 'w') as fh_out:
    fh_out.write('rsID\talt_allele_frequency\n')
    lst_variant_IDs = df.columns[2:]  # Skip ID and status in the header
    for i in range(len(lst_variant_IDs)):
        line = lst_variant_IDs[i] + '\t' + str(lst_alt_allele_frequency[i]) + '\n'
        fh_out.write(line)

# ------------------ Q2. Hardy‑Weinberg Equilibrium ------------------
# Only calculate HWE for variants with MAF>0.05
df_variants_HWE = pd.read_csv('AF.txt', sep='\t')
threshold_maf = 0.05
mask_drop = (df_variants_HWE['alt_allele_frequency'] <= threshold_maf) | (
        df_variants_HWE['alt_allele_frequency'] >= 1 - threshold_maf)
df_variants_HWE.drop(labels=(df_variants_HWE[mask_drop]).index, inplace=True)  # Drop where MAF<=0.05

if verbose:
    print('\n\n=======================================\nQ2. HWE')
    print('- Only calculate HWE for SNPs with MAF >', threshold_maf)
    print('- Number of dropped SNPs due to small MAF: ', mask_drop.shape[0])
    print('- Number of SNPs with MAF>0,05: ', df_variants_HWE.shape[0])

# Step 1. Calculate minor allele frequency
df_variants_HWE['MAF'] = df_variants_HWE[
    'alt_allele_frequency']  # Create a column of MAF, set initial values to 'alt_allele_frequency' column
mask_AF_gt_05 = df_variants_HWE['alt_allele_frequency'] >= 0.5  # mask of SNPs with allele frequency >0.5

# Only this or df_variants_HWE.loc[mask_AF_gt_05, ['MAF']] works, could be a bug of pandas
# df_variants_HWE.loc[mask_AF_gt_05]['MAF'] does not work, even though it also returns a Series for assignment
df_variants_HWE.loc[mask_AF_gt_05, 'MAF'] = 1 - df_variants_HWE.loc[mask_AF_gt_05]['MAF']
# df_variants_HWE.to_csv('AF_and_alt_AF_MAF_gt_005_only.txt', sep='\t', index=False) # For debugging purpose

# Step 2. Count and calculate frequencies of 0/0, 0/1 (could be 1/0) and 1/1 in original DataFrame df
mask_00 = df.iloc[:, 2:] == '0/0'  # Skip the first two columns (ID and status)
mask_01 = (df.iloc[:, 2:] == '0/1') | (df.iloc[:, 2:] == '1/0')
mask_11 = df.iloc[:, 2:] == '1/1'

# Store observed frequencies of 0/0 (0), 0/1 (1) and 1/1 (2) to this dataframe
df_pair_count_frequency = pd.concat([mask_00.sum(), mask_01.sum(), mask_11.sum()], axis=1)
df_pair_count_frequency.columns = ['0_count_obs', '1_count_obs', '2_count_obs']
'''
# Do not need to calculate frequencies
df_pair_count_frequency['0_freq_obs'] = df_pair_count_frequency['0_count_obs']/number_of_individuals # divided by number of subjects
df_pair_count_frequency['1_freq_obs'] = df_pair_count_frequency['1_count_obs']/number_of_individuals
df_pair_count_frequency['2_freq_obs'] = df_pair_count_frequency['2_count_obs']/number_of_individuals
'''

if verbose:
    print('- Count and frequency of 0 (0/0), 1 (0/1) and 2 (1/1) for all variants:')
    print(df_pair_count_frequency.head())

# Step 3. Calculate HWE p values
# Note: since minor and major allele frequencies follow HWE,
#       we use 1 as degree of freedom even though there are 3 values 0/0, 0/1 and 1/1
# Use df_variants_HWE and df_pair_count_frequency in this step
# - df_variants_HWE: contains MAFs for SNPs with MAF >0.05
# - df_pair_count_frequency: contains observed counts and frequencies of 0/0, 0/1 and 1/1

# Calculate expected frequencies of 0, 1 and 2 using df_pair_count_frequency.
# Note that 0 is ref allele, not necessarily major allele, so should not use 1-MAF
df_variants_HWE['0_count_exp'] = (1 - df_variants_HWE['alt_allele_frequency']) ** 2 * number_of_individuals
df_variants_HWE['1_count_exp'] = 2 * df_variants_HWE['MAF'] * (1 - df_variants_HWE['MAF']) * number_of_individuals
df_variants_HWE['2_count_exp'] = df_variants_HWE['alt_allele_frequency'] ** 2 * number_of_individuals

# Merge df_variants_HWE and df_pair_count_frequency to merge observed and expected counts into the same dataframe
df_merged_exp_obs_freq = df_variants_HWE.merge(df_pair_count_frequency, left_on='rsID', right_index=True)
# df_merged_exp_obs_freq.to_csv('optional_obs_exp_counts_gt_005_only_for_HWE.txt', sep='\t', index=False) # For debugging

if verbose:
    print('- SNPs with MAF>0.5, observed and expected counts:')
    print(df_merged_exp_obs_freq.head())

# Apply chisquare() to calculated p values
# Note: Use ddof=1 to set degree of freedom to 1, otherwise it will be observation-1=3-1=2
#       Second item in the return value of chisquare() is p value
df_merged_exp_obs_freq['HWE_p_vals'] = df_merged_exp_obs_freq.iloc[:, 3:].apply(
    lambda x: chisquare(f_exp=x[:3], f_obs=x[-3:], ddof=1)[1], axis=1)

# Step 4. Output into HWE.txt
df_merged_exp_obs_freq[['rsID', 'HWE_p_vals']].to_csv('HWE.txt', sep='\t', index=False)

# ------------------ Q3. Linkage Disequilibrium ------------------
# Can also read in AF.txt file to create this dataframe
df_variants_af = pd.DataFrame({'rsID': lst_variant_IDs, 'alt_allele_frequency': lst_alt_allele_frequency})
threshold_maf = 0.05
mask_LD_drop = (df_variants_af['alt_allele_frequency'] < threshold_maf) | (
        df_variants_af['alt_allele_frequency'] > 1 - threshold_maf)
df_variants_af.drop(labels=(df_variants_af[mask_LD_drop]).index, inplace=True)  # Drop where MAF<0.05
number_of_SNPS_to_use = 100  # Only calculate LD for the first 100 SNPs (4950 pairs)
if verbose:
    print('\n\n=======================================\nQ3. LD')
    print('- Calculate pairwise D, D\' and R2 of the first', number_of_SNPS_to_use, 'SNPs')


# ----------- Helper functions -----------
# Estimate haplotype frequencies of snp1-snp2 with EM algorithm
# Return haplotype frequencies of 0|0, 0|1, 1|0 and 1|1
# Original dataframe df is also used in this function
def estimate_haplotype_frequency(snp1, snp2):
    # Get genotype counts of the two SNPs (two-locus genotype)
    # Note:
    # 1. There are 3 genotypes for a single SNP: 0 (00), 1(01) and 2(11),
    #    so for two loci there are 3*3=9 possible two-locus genotypes
    # 2. Explanation of naming below: n_00 means SNP1 and SNP2 are both 00. n_12 means SNP1 is 01 and SNP2 is 11
    # 3. n_11 is the only one with ambiguous haplotype. Since both SNPs are 01, we can have two possible haplotypes:
    #    (1) SNP1-SNP2: 0-1 and 1-0; or (2) SNP1-SNP2: 0-0 and 1-1
    # 4. Use EM algorithm to get haplotype frequency
    n_00 = df[(df[snp1] == '0/0') & (df[snp2] == '0/0')].shape[0]  # haplotype (SNP1|SNP2): 00|00
    n_01 = df[(df[snp1] == '0/0') & (df[snp2] == '0/1')].shape[0]  # haplotype (SNP1|SNP2): 00|01
    n_02 = df[(df[snp1] == '0/0') & (df[snp2] == '1/1')].shape[0]  # haplotype (SNP1|SNP2): 01|01
    n_10 = df[(df[snp1] == '0/1') & (df[snp2] == '0/0')].shape[0]  # haplotype (SNP1|SNP2): 00|10
    n_11 = df[(df[snp1] == '0/1') & (df[snp2] == '0/1')].shape[0]  # haplotype (SNP1|SNP2): 01|10 or 00|11
    n_12 = df[(df[snp1] == '0/1') & (df[snp2] == '1/1')].shape[0]  # haplotype (SNP1|SNP2): 01|11
    n_20 = df[(df[snp1] == '1/1') & (df[snp2] == '0/0')].shape[0]  # haplotype (SNP1|SNP2): 10|10
    n_21 = df[(df[snp1] == '1/1') & (df[snp2] == '0/1')].shape[0]  # haplotype (SNP1|SNP2): 10|11
    n_22 = df[(df[snp1] == '1/1') & (df[snp2] == '1/1')].shape[0]  # haplotype (SNP1|SNP2): 11|11

    # ---- E step ----
    # Initial counts
    # Calculate expected counts of in n_11, ie. when SNP1-SNP2 is 01/01 (hets in both loci)
    # !!!! I am not completely sure about this one, since initially all haplotype frequencies are unknown
    # !!!! For initial counts, assume half of n_11 is 0|1:1|0, the other half is 0|0:1|1
    # !!!! Probably won't matter since we are "guessing" anyway (yes! I tested with other values)
    n_01_10 = n_11 * 0.5  # 01|10
    n_00_11 = n_11 * 0.5  # 00|11

    # haplotype counts of SNP1-SNP2
    # There are 2*2=4 types of haplotype in total (00, 01, 10, 11)
    haplotype_count_00 = n_00 * 2 + n_01 + n_10 + n_00_11
    haplotype_count_01 = n_01 + n_02 * 2 + n_01_10 + n_12
    haplotype_count_10 = n_10 + n_01_10 + n_20 * 2 + n_21
    haplotype_count_11 = n_00_11 + n_12 + n_21 + n_22 * 2

    hap_freq_00 = haplotype_count_00 / (number_of_individuals * 2)
    hap_freq_01 = haplotype_count_01 / (number_of_individuals * 2)
    hap_freq_10 = haplotype_count_10 / (number_of_individuals * 2)
    hap_freq_11 = haplotype_count_11 / (number_of_individuals * 2)

    prev_hap_freq_00 = hap_freq_00
    prev_hap_freq_01 = hap_freq_01
    prev_hap_freq_10 = hap_freq_10
    prev_hap_freq_11 = hap_freq_11

    if output_log:
        with open('optional_EM_iterations_values_at_each_step.txt', 'a') as fh:
            fh.write('\n-----------------' + snp1 + ' ' + snp2 + '-----------------\n')
            fh.write(' n_00\tn_01\tn_02\tn_10\tn_11\tn_12\tn_20\tn_21\tn_22\n ')
            for val in [n_00, n_01, n_02, n_10, n_11, n_12, n_20, n_21, n_22]:
                fh.write(str(val) + '\t')
            fh.write('\n Initial expected haplotype counts (0|0, 0|1, 1|0, 1|1): ')
            for val in [haplotype_count_00, haplotype_count_01, haplotype_count_10]:
                fh.write(str(val) + ',')
            fh.write(str(haplotype_count_11) + '\n')

    # ---- M step ----
    # Update haplotype frequencies with new haplotype counts
    # Stop iteration when changes <10^-5, or reach 1 million iterations
    # From testing, estimated haplotype frequencies usually converge within 40 iterations in this assignment
    number_of_iterations_to_converge = 0
    while number_of_iterations_to_converge < 1000000:
        number_of_iterations_to_converge += 1
        # Use updated haplotype frequencies to calculate expected counts of 0|1:1|0 and 0|0:1|1 within n_11
        n_01_10 = n_11 * (2 * hap_freq_01 * hap_freq_10) / (
                2 * hap_freq_01 * hap_freq_10 + 2 * hap_freq_00 * hap_freq_11)  # 01|10
        n_00_11 = n_11 - n_01_10  # 00|11

        # Update haplotype counts with new haplotype frequencies, then update haplotype frequencies (iteration)
        haplotype_count_00 = n_00 * 2 + n_01 + n_10 + n_00_11
        haplotype_count_01 = n_01 + n_02 * 2 + n_01_10 + n_12
        haplotype_count_10 = n_10 + n_01_10 + n_20 * 2 + n_21
        haplotype_count_11 = n_00_11 + n_12 + n_21 + n_22 * 2

        hap_freq_00 = haplotype_count_00 / (number_of_individuals * 2)
        hap_freq_01 = haplotype_count_01 / (number_of_individuals * 2)
        hap_freq_10 = haplotype_count_10 / (number_of_individuals * 2)
        hap_freq_11 = haplotype_count_11 / (number_of_individuals * 2)

        # Stop iteration when converged (changes of any haplotype frequency less than stop threshold)
        stop_threshold = 1e-6
        if abs(prev_hap_freq_00 - hap_freq_00) < stop_threshold or \
                abs(prev_hap_freq_01 - hap_freq_01) < stop_threshold or \
                abs(prev_hap_freq_10 - hap_freq_10) < stop_threshold or \
                abs(prev_hap_freq_11 - hap_freq_11) < stop_threshold:
            if output_log:
                with open('optional_EM_iterations_values_at_each_step.txt', 'a') as fh:
                    fh.write('# ' + str(number_of_iterations_to_converge) + ' iteration (final):\n')
                    fh.write(' haplotype count: ')
                    for val in [haplotype_count_00, haplotype_count_01, haplotype_count_10]:
                        fh.write(str(val) + ', ')
                    fh.write(str(haplotype_count_11) + '\n')
                    fh.write(' haplotype freq: ')
                    for val in [hap_freq_00, hap_freq_01, hap_freq_10]:
                        fh.write(str(val) + ', ')
                    fh.write(str(hap_freq_11) + '\n')
            break

        if output_log:
            with open('optional_EM_iterations_values_at_each_step.txt', 'a') as fh:
                fh.write('# ' + str(number_of_iterations_to_converge) + ' iteration:\n')
                fh.write(' haplotype count: ')
                for val in [haplotype_count_00, haplotype_count_01, haplotype_count_10]:
                    fh.write(str(val) + ', ')
                fh.write(str(haplotype_count_11) + '\n')
                fh.write(' haplotype freq: ')
                for val in [hap_freq_00, hap_freq_01, hap_freq_10]:
                    fh.write(str(val) + ', ')
                fh.write(str(hap_freq_11) + '\n')

        prev_hap_freq_00 = hap_freq_00
        prev_hap_freq_01 = hap_freq_01
        prev_hap_freq_10 = hap_freq_10
        prev_hap_freq_11 = hap_freq_11

    return hap_freq_00, hap_freq_01, hap_freq_10, hap_freq_11, number_of_iterations_to_converge

# This calculate LD scores using:
# - D(AB) = P(AB) – P(A)P(B)
# - D'(AB) = D(AB)/Dmax(AB) (refer to slides about how to define Dmax(AB))
# - R2 = [D(AB)]^2/[p(A)p(a)p(B)p(b)]
# Note: (1) P(AB) is estimated by EM algorithm
#       (2) A and B are major alleles of each marker (snp1 and snp2)
# Parameters: - snp1_index, snp2_index are indices of SNP1 and SNP2 in df_variants_af.
#             - hap_freq_00, hap_freq_01, hap_freq_10, hap_freq_11: haplotype frequencies of 0|, 0|1, 1|0, and 1|1
#             df_variants_af is a dataframe containing rsID, alternative allele frequencies
# Return: LD scores D, D', R2 of SNP1 and SNP2
def get_LD_scores(snp1_index, snp2_index, hap_freq_00, hap_freq_01, hap_freq_10, hap_freq_11):
    # Get major and minor allele frequencies of SNP1 and SNP2 from variable df_variants_af
    if df_variants_af.iloc[snp1_index, 1] > 0.5:  # In this case alt allele is major allele, ie. '1' is major allele
        snp1_major_af = df_variants_af.iloc[snp1_index, 1]
        snp1_major_allele_number = '1'  # Keep track of which number represents major allele
    else:
        snp1_major_af = 1 - df_variants_af.iloc[snp1_index, 1]
        snp1_major_allele_number = '0'
    if df_variants_af.iloc[snp2_index, 1] > 0.5:
        snp2_major_af = df_variants_af.iloc[snp2_index, 1]
        snp2_major_allele_number = '1'
    else:
        snp2_major_af = 1 - df_variants_af.iloc[snp2_index, 1]
        snp2_major_allele_number = '0'

    snp1_maf = 1 - snp1_major_af
    snp2_maf = 1 - snp2_major_af
    # Get haplotype frequency of when both markers are major alleles
    p_AB = eval('hap_freq_' + snp1_major_allele_number + snp2_major_allele_number)
    ld_D = p_AB - snp1_major_af * snp2_major_af
    ld_r2 = ld_D ** 2 / (snp1_major_af * snp1_maf * snp2_major_af * snp2_maf)

    # D'(AB) = D(AB)/Dmax(AB)
    if ld_D > 0:
        # Dmax(AB) = min[P(A)P(b), P(a)P(B)]
        d_max = min([snp1_major_af * snp2_maf, snp1_maf * snp2_major_af])
    else:
        # Dmax(AB) = min[P(A)P(B), P(a)P(b)]
        d_max = min([snp1_major_af * snp2_major_af, snp1_maf * snp2_maf])
    ld_D_prime = ld_D / d_max

    # !!!!!! Not very sure: Return absolute values of D and D'
    # return abs(ld_D), abs(ld_D_prime), ld_r2
    return ld_D, ld_D_prime, ld_r2
# ----------- End of helper functions -----------

# Calculate LD (D, D' and R2) for the first 100 pairs of SNPs
df_variants_af.reset_index(drop=True, inplace=True)  # Reset index since it was modified due to previous drop MAF<0.05
df_variants_af.drop(labels=[i for i in range(number_of_SNPS_to_use, df_variants_af.shape[0])], inplace=True)

if output_log:  # Optional: output number of iterations to converge and haplotype frequencies for each pair of SNPs
    with open('optional_EM_iterations_to_converge_final_haplotype_freq.txt', 'w') as fh:
        fh.write('SNP1\tSNP2\tnumber_of_iterations_to_converge\t0|0_freq\t0|1_freq\t1|0_freq\t1|1_freq\n')
    with open('optional_EM_iterations_values_at_each_step.txt', 'w') as fh:  # Output details of each pair
        fh.write('')

lst_SNPs_to_calculate_LD = df_variants_af['rsID']  # Actually this is a Series, not a list
count = 0  # Track progres

# Create required output files
with open('LD_D.txt', 'w') as fh:
    fh.write('SNP1\tSNP2\tLD_D\n')
with open('LD_Dprime.txt', 'w') as fh:
    fh.write('SNP1\tSNP2\tLD_Dprime\n')
with open('LD_r2.txt', 'w') as fh:
    fh.write('SNP1\tSNP2\tLD_r2\n')


# ========================= Multiprocessing approach (faster) ================================
# Get all pairs of SNPs to pass to EM algorithm
lst_pairs_of_SNPs = []
lst_pairs_of_SNPs_index = [] # Collect indices of SNP pair, this is for LD calculation
for i in range(number_of_SNPS_to_use):
    snp1 = lst_SNPs_to_calculate_LD[i]
    for j in range(i, number_of_SNPS_to_use):
        snp2 = lst_SNPs_to_calculate_LD[j]
        lst_pairs_of_SNPs.append((snp1, snp2))
        lst_pairs_of_SNPs_index.append((i,j))
for pairs in range(0, len(lst_pairs_of_SNPs), 50):  # Process 50 pairs of SNPs with multiprocessing at a time
    # To be safe, use the max number of cores to do multi processing, unless only one core available
    if multiprocessing.cpu_count() == 1:
        number_of_cores_to_use = 1
    else:
        number_of_cores_to_use = multiprocessing.cpu_count() - 1
    # Get haplotype frequencies
    with multiprocessing.Pool(number_of_cores_to_use) as p:
        if pairs < len(lst_pairs_of_SNPs) - 50:
            # haplotype_freq is a list of tuples, each tuple contains haplotype frequencies of 0|0, 0|1, 1|0 and 1|1 of a SNP pair
            pairs_of_SNPs_index = lst_pairs_of_SNPs_index[pairs:pairs + 50]
            pairs_of_SNPs = lst_pairs_of_SNPs[pairs:pairs + 50]
            haplotype_freq = p.starmap(estimate_haplotype_frequency, pairs_of_SNPs)
        else:  # If less than 50 pairs left
            pairs_of_SNPs_index = lst_pairs_of_SNPs_index[pairs:len(lst_pairs_of_SNPs)]
            pairs_of_SNPs = lst_pairs_of_SNPs[pairs:len(lst_pairs_of_SNPs)]
            haplotype_freq = p.starmap(estimate_haplotype_frequency, pairs_of_SNPs)

        # Calculate LD scores
        parameters_to_pass = []
        for k in range(len(haplotype_freq)):
            # Create a list of values to pass to function get_LD_scores. The list has following format:
            # [(snp1_index, snp2_index, hap_freq_00, hap_freq_01, hap_freq_10, hap_freq_11), (...), ...]
            parameters_to_pass.append(pairs_of_SNPs_index[k] + haplotype_freq[k][:-1])
        ld_scores = p.starmap(get_LD_scores, parameters_to_pass)

        # Output into required files
        fh_ld_d = open('LD_D.txt', 'a')
        fh_ld_d_prime = open('LD_Dprime.txt', 'a')
        fh_ld_r2 = open('LD_r2.txt', 'a')
        for val_index in range(len(ld_scores)):
            snp1 = pairs_of_SNPs[val_index][0]
            snp2 = pairs_of_SNPs[val_index][1]
            fh_ld_d.write(snp1 + '\t' + snp2 + '\t' + str(ld_scores[val_index][0]) + '\n')
            fh_ld_d_prime.write(snp1 + '\t' + snp2 + '\t' + str(ld_scores[val_index][1]) + '\n')
            fh_ld_r2.write(snp1 + '\t' + snp2 + '\t' + str(ld_scores[val_index][2]) + '\n')
        fh_ld_d.close()
        fh_ld_d_prime.close()
        fh_ld_r2.close()

    if (pairs + 50) % 500 == 0:
        print('.', pairs + 50, 'pairs of SNPs processed')
    else:
        print('.', end='', flush=True)
print('\nIn total:', len(lst_pairs_of_SNPs), 'pairs of SNPs processed')
# ========================== END Multiprocessing ===============================

# Use multiprocessing above instead of the block of code below. They generate the same result though
'''
for i in range(number_of_SNPS_to_use):
    snp1 = lst_SNPs_to_calculate_LD[i]
    for j in range(i, number_of_SNPS_to_use):
        count += 1
        snp2 = lst_SNPs_to_calculate_LD[j]

        # ----- Call helper function to get haplotype frequencies with EM algorithm -----
        hap_freq_00, hap_freq_01, hap_freq_10, hap_freq_11, number_of_iterations_to_converge = estimate_haplotype_frequency(
            snp1, snp2)
        if verbose:  # Keep console busy
            if count == 1:  # Process first pair of SNPs
                print('- Estimating haplotype frequencies with EM algorithm, then calculate LD scores')
            elif count % 500 == 0:
                print('.', count, 'pairs processed')
            elif count % 50 == 0:
                print('.', end='', flush=True)

        if output_log:
            with open('optional_EM_iterations_to_converge_final_haplotype_freq.txt', 'a') as fh:
                fh.write(snp1 + '\t' + snp2 + '\t' + str(number_of_iterations_to_converge + 1) + '\t' + str(hap_freq_00) \
                         + '\t' + str(hap_freq_01) + '\t' + str(hap_freq_10) + '\t' + str(hap_freq_11) + '\n')

        # ----- Calculate LD scores -----
        # My plan is to use major alleles in LD calculation (since they have higher allele frequencies)
        ld_D, ld_D_prime, ld_r2 = get_LD_scores(i, j, hap_freq_00, hap_freq_01, hap_freq_10, hap_freq_11)

        # Output in to required files
        with open('LD_D.txt', 'a') as fh:
            fh.write(snp1 + '\t' + snp2 + '\t' + str(ld_D) + '\n')
        with open('LD_Dprime.txt', 'a') as fh:
            fh.write(snp1 + '\t' + snp2 + '\t' + str(ld_D_prime) + '\n')
        with open('LD_r2.txt', 'a') as fh:
            fh.write(snp1 + '\t' + snp2 + '\t' + str(ld_r2) + '\n')

        # print('\n-------------', snp1, snp2,'-------------')
        # print(' alt AF (snp1, snp2):', df_variants_af.iloc[i, 1], df_variants_af.iloc[j, 1])
        # print(' Hap freq (00, 01, 10, 11):', hap_freq_00, hap_freq_01, hap_freq_10, hap_freq_11)
        # print(' LD score (D, D\', r2):', ld_D, ld_D_prime, ld_r2)

print('\n Total:', count, 'pairs of SNPs processed of the first', number_of_SNPS_to_use, 'SNPs')
'''

# ------------------ Q4. Principal component analysis ------------------
if verbose:
    print('\n\n=======================================\nQ4. PCA')
# Create a new dataframe (df_additive_genotype) from original dataframe df, convert 0/0, 0/1, 1/1 to 0, 1, 2
df_additive_genotype = df.replace(to_replace='0/0', value=0).copy()
df_additive_genotype.replace(to_replace='0/1', value=1, inplace=True)
df_additive_genotype.replace(to_replace='1/1', value=2, inplace=True)

if verbose:
    print('- Replace 0/0, 0/1, 1/1 with 0, 1, 2 in a new DataFrame (df_additive_genotype):')
    print(df_additive_genotype.head())

# Perform PCA
pca_result = PCA()  # Keep the 150 PCs by default (since sample size=150, smaller than SNP size)
pc_vals = pca_result.fit_transform(df_additive_genotype.iloc[:, 2:])  # Fit a PCA model and calculate PCs
if verbose:
    print('- Perform PCA, % variance explained by the first 10 PCs are:')
    for var in pca_result.explained_variance_ratio_[:10] * 100:
        print('\t{:.2f}%'.format(var), end='')
    print()

if save_pca_result:
    # Save PC values
    df_PCs = pd.DataFrame(pc_vals)
    df_PCs.columns = ['PC' + str(x + 1) for x in df_PCs.columns]
    df_PCs.to_csv('PCs.txt', sep='\t', index=False)

    # Save variance explained by each PC
    with open('PCA_variance_ratio_explained_by_each_PC.txt', 'w') as fh:
        fh.write('PC_number' + '\t' + 'variance_explained' + '\t' + 'variance_explained_ratio' + '\n')
        for i in range(len(pca_result.explained_variance_ratio_)):
            fh.write(str(i + 1) + '\t' + str(pca_result.explained_variance_[i]) + '\t' + str(
                pca_result.explained_variance_ratio_[i]) + '\n')

print('\n# Run finished in {:.2f} minutes'.format((time.time()-start_time)/60))