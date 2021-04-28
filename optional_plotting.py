import matplotlib.pyplot as plt
import pandas as pd
fn = 'saved_results/optional_EM_iterations_to_converge_final_haplotype_freq.txt'
df_iterations = pd.read_csv(fn, sep='\t')
fig_it, ax_it = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), dpi=150, gridspec_kw={'width_ratios': [2, 1]})
ax_it[0].hist(df_iterations['number_of_iterations_to_converge'], rwidth=0.8)
ax_it[0].set_xlabel('Number of iterations to converge')
ax_it[0].set_ylabel('Counts')
ax_it[1].boxplot(df_iterations['number_of_iterations_to_converge'], labels=[''])
ax_it[1].set_ylabel('Number of iterations to converge')
fig_it.suptitle('EM')
fig_it.tight_layout()
fig_it.savefig('saved_results/optional_EM_number_of_iterations_to_converge.jpeg')
