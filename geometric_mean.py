"""
Peptide Identification Combination Script
Combines results from three search engines (Comet, X!Tandem, MS-GF+) using geometric mean.

Formula: p_combined = (p1 × p2 × p3)^(1/3)
Missing values are filled with 1.0 before combination.
"""

import pandas as pd
import numpy as np
from pathlib import Path

def read_and_normalize_results(file_path, engine_name):
    """Read search engine results and normalize p-values"""
    if engine_name == 'comet':
        df = pd.read_csv(file_path, sep='\t', skiprows=1)
        df.columns = df.columns.str.lower()
        df = df.rename(columns={
            'scannum': 'scan',
            'peptide': 'plain_peptide',
            'expectscore': 'e-value'
        })
    elif engine_name == 'xtandem':
        df = pd.read_csv(file_path)
        df = df.rename(columns={
            'scannum': 'scan',
            'peptide': 'plain_peptide',
            'e.value': 'e-value'
        })
    elif engine_name == 'msgf':
        df = pd.read_csv(file_path)
        df = df.rename(columns={
            'ScanNum': 'scan',
            'SpecEValue': 'e-value'
        })
    
    # Remove NaN peptides
    df = df.dropna(subset=['plain_peptide'])
    
    # Calculate n (number of matches per scan)
    df['n'] = df.groupby('scan')['scan'].transform('count')
    
    # Calculate p-value
    df['p.value'] = df['e-value'] / df['n']
    df = df[df['p.value'] <= 1]
    
    # Calculate adjusted p-value: 1 - (1 - p)^n
    df['p_adjusted'] = 1 - (1 - df['p.value']) ** df['n']
    df = df[(df['p_adjusted'] >= 0) & (df['p_adjusted'] <= 1)]
    
    # Keep only best p-value per scan
    df = df.sort_values('p_adjusted').drop_duplicates(subset=['scan'], keep='first')
    
    return df[['scan', 'plain_peptide', 'p_adjusted']]

def combine_results(df1, df2, df3):
    """Combine results from three search engines using geometric mean"""
    # Create PSM (Peptide-Spectrum Match) identifiers
    for df in [df1, df2, df3]:
        df['psm'] = df['scan'].astype(str) + '-' + df['plain_peptide']
    
    # Merge all results
    merged = pd.merge(df1, df2, on='psm', how='outer', suffixes=('_1', '_2'))
    merged = pd.merge(merged, df3, on='psm', how='outer')
    merged.columns = ['psm', 'scan_1', 'peptide_1', 'p1', 'scan_2', 'peptide_2', 'p2', 
                      'scan', 'plain_peptide', 'p3']
    
    # Fill missing values with 1.0
    merged[['p1', 'p2', 'p3']] = merged[['p1', 'p2', 'p3']].fillna(1.0)
    
    # Extract scan number from PSM
    merged['scan'] = merged['psm'].str.split('-').str[0].astype(float)
    merged['peptide'] = merged['psm'].str.split('-').str[1]
    
    # Combine p-values using geometric mean: (p1 * p2 * p3)^(1/3)
    merged['p_combined'] = (merged['p1'] * merged['p2'] * merged['p3']) ** (1/3)
    
    # Keep best result per scan
    result = merged.sort_values('p_combined').drop_duplicates(subset=['scan'], keep='first')
    
    return result[['scan', 'peptide', 'p_combined', 'p1', 'p2', 'p3']]

def main():
    # File paths
    comet_file = "searched_results/comet/pool1.txt"
    xtandem_file = "searched_results/tandem/pool1.csv"
    msgf_file = "searched_results/msgf/pool1.csv"
    
    # Read and normalize results from each engine
    print("Reading Comet results...")
    df_comet = read_and_normalize_results(comet_file, 'comet')
    
    print("Reading X!Tandem results...")
    df_xtandem = read_and_normalize_results(xtandem_file, 'xtandem')
    
    print("Reading MS-GF+ results...")
    df_msgf = read_and_normalize_results(msgf_file, 'msgf')
    
    # Combine results using geometric mean
    print("Combining results using geometric mean...")
    combined = combine_results(df_comet, df_xtandem, df_msgf)
    
    # Save results
    combined.to_csv("combined_results/geometric_mean/pool1_combine.csv", index=False)
    print(f"Results saved. Total PSMs: {len(combined)}")
    print(f"Average combined p-value: {combined['p_combined'].mean():.4f}")
    
    # Show top 10 results
    print("\nTop 10 peptide identifications:")
    print(combined.head(10)[['scan', 'peptide', 'p_combined']].to_string(index=False))

if __name__ == "__main__":
    main()