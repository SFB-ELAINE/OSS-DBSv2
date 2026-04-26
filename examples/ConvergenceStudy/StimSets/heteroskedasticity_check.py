#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 18:08:46 2026

@author: forel
"""

import pandas as pd
import numpy as np

# Load the data
best_df = pd.read_csv('Results_PAM_best_overview.csv')
ref_df = pd.read_csv('Results_PAM_hp_material_refinement_overview.csv')

# Flatten all values from both dataframes into 1D arrays
best_values = best_df.values.flatten()
ref_values = ref_df.values.flatten()

# Filter out entries where the reference (best) activation is zero
mask = best_values > 0
best_filtered = best_values[mask]
ref_filtered = ref_values[mask]

# Compute absolute errors for the filtered values
abs_errors = np.abs(best_filtered - ref_filtered)

# Compute the correlation coefficient
# We use pandas Series to handle the correlation easily
correlation = pd.Series(best_filtered).corr(pd.Series(abs_errors))

print(f"Total non-zero samples: {len(best_filtered)}")
print(f"Global Correlation: {correlation}")

# Prepare a small summary for the user
summary_data = {
    "Metric": ["Total non-zero entries", "Global Correlation", "Mean Activation (non-zero)", "Mean Abs Error (non-zero)"],
    "Value": [len(best_filtered), correlation, best_filtered.mean(), abs_errors.mean()]
}
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv('global_correlation_results.csv', index=False)
print(summary_df)