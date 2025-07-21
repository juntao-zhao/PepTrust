#!/usr/bin/env python3
"""
Calculate precision for peptide identification methods in PepTrust.
Evaluates each scan-peptide identification against ground truth.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def load_ground_truth(ground_truth_file):
    """Load ground truth peptides from the Excel file."""
    print("Loading ground truth...")
    try:
        # Read the Excel file, specifically the "Peptide identifications" sheet
        df = pd.read_excel(ground_truth_file, sheet_name='Peptide identifications')
        
        # Extract unique peptides from the "Sequence" column
        ground_truth_peptides = set(df['Sequence'].dropna().unique())
        print(f"Loaded {len(ground_truth_peptides)} unique peptides from ground truth")
        
        return ground_truth_peptides
    except Exception as e:
        print(f"Error loading ground truth: {e}")
        return set()


def load_search_results(engine, pool_num):
    """Load scan-peptide identifications from search engine results."""
    file_path = f"searched_results/{engine}/pool{pool_num}.csv"
    try:
        df = pd.read_csv(file_path)
        # Return list of (scan, peptide) tuples
        identifications = []
        for _, row in df.iterrows():
            if pd.notna(row['Peptide']):
                identifications.append((row['scan'], row['Peptide']))
        return identifications
    except Exception as e:
        print(f"Error loading {engine} pool{pool_num}: {e}")
        return []


def load_combination_results(method, pool_num):
    """Load scan-peptide identifications from combination method results."""
    # Define file paths based on method
    if method == 'geometric_mean':
        file_path = f"combined_results/geometric_mean/pool{pool_num}_combine.csv"
    elif method == 'msblender':
        file_path = f"combined_results/msblender/pool{pool_num}.csv"
    elif method == 'voting':
        file_path = f"combined_results/voting/pool_{pool_num}.csv"
    
    try:
        df = pd.read_csv(file_path)
        identifications = []
        
        # Special handling for geometric mean files
        if method == 'geometric_mean':
            # The psm column contains "scan-peptide" format
            if 'psm' in df.columns and 'scan' in df.columns:
                for _, row in df.iterrows():
                    if pd.notna(row['psm']) and pd.notna(row['scan']):
                        # Extract peptide from psm column (format: "scan-peptide")
                        psm_parts = str(row['psm']).split('-', 1)
                        if len(psm_parts) == 2:
                            scan = str(row['scan'])  # Use the scan column directly
                            peptide = psm_parts[1]   # Get peptide part after the dash
                            identifications.append((scan, peptide))
            else:
                print(f"Warning: Expected 'psm' and 'scan' columns not found in {file_path}")
                return []
        
        # Handle MSBlender and voting methods
        else:
            # Extract scan-peptide pairs based on available columns
            if 'scan' in df.columns:
                scan_col = 'scan'
            elif 'Spectrum' in df.columns:
                # For MSBlender, extract scan number from Spectrum column
                df['scan'] = df['Spectrum'].str.extract(r'\.(\d+)\.')
                scan_col = 'scan'
            else:
                print(f"Warning: No scan column found in {file_path}")
                return []
            
            # Find peptide column
            if 'plain_peptide' in df.columns:
                peptide_col = 'plain_peptide'
            elif 'Peptide' in df.columns:
                peptide_col = 'Peptide'
            else:
                print(f"Warning: No peptide column found in {file_path}")
                return []
            
            # Extract identifications
            for _, row in df.iterrows():
                if pd.notna(row[peptide_col]) and pd.notna(row[scan_col]):
                    identifications.append((str(row[scan_col]), row[peptide_col]))
        
        return identifications
    except Exception as e:
        print(f"Error loading {method} pool{pool_num}: {e}")
        return []


def calculate_precision(identifications, ground_truth_peptides):
    """
    Calculate precision for scan-peptide identifications.
    Each identification is evaluated independently.
    """
    if len(identifications) == 0:
        return 0.0, 0, 0, 0
    
    true_positives = 0
    false_positives = 0
    
    # Evaluate each scan-peptide identification
    for scan, peptide in identifications:
        if peptide in ground_truth_peptides:
            true_positives += 1
        else:
            false_positives += 1
    
    total_identifications = true_positives + false_positives
    precision = true_positives / total_identifications if total_identifications > 0 else 0.0
    
    return precision, true_positives, false_positives, total_identifications


def main():
    """Main function to calculate precision for all methods and pools."""
    
    # Define methods and pools
    search_engines = ['comet', 'msgf', 'tandem']
    combination_methods = ['msblender', 'voting', 'geometric_mean']
    all_methods = search_engines + combination_methods
    pools = [1, 2, 3, 4, 5]
    
    # Load ground truth
    ground_truth_file = '41592_2017_BFnmeth4153_MOESM250_ESM.xlsx'
    ground_truth_peptides = load_ground_truth(ground_truth_file)
    
    if not ground_truth_peptides:
        print("Failed to load ground truth. Exiting.")
        return
    
    # Initialize results storage
    results = []
    
    # Calculate precision for each method and pool
    print("\nCalculating precision for each method and pool...")
    print("Note: Each scan-peptide identification is evaluated independently")
    
    for pool_num in pools:
        print(f"\nProcessing Pool {pool_num}...")
        
        # Process search engines
        for engine in search_engines:
            identifications = load_search_results(engine, pool_num)
            precision, tp, fp, total = calculate_precision(identifications, ground_truth_peptides)
            
            results.append({
                'Pool': f'pool{pool_num}',
                'Method': engine.upper() if engine == 'msgf' else engine.capitalize(),
                'True_Positives': tp,
                'False_Positives': fp,
                'Total_Identifications': total,
                'Precision': precision
            })
            
            print(f"  {engine}: {total} identifications, {tp} TP, {fp} FP, precision = {precision:.4f}")
        
        # Process combination methods
        for method in combination_methods:
            identifications = load_combination_results(method, pool_num)
            precision, tp, fp, total = calculate_precision(identifications, ground_truth_peptides)
            
            results.append({
                'Pool': f'pool{pool_num}',
                'Method': 'MSBlender' if method == 'msblender' else method.replace('_', ' ').title(),
                'True_Positives': tp,
                'False_Positives': fp,
                'Total_Identifications': total,
                'Precision': precision
            })
            
            print(f"  {method}: {total} identifications, {tp} TP, {fp} FP, precision = {precision:.4f}")
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Save detailed results
    output_file = 'precision/peptide_identification_precision_results.csv'
    results_df.to_csv(output_file, index=False)
    print(f"\nDetailed results saved to: {output_file}")
    
    # Create and save summary statistics
    summary_df = results_df.pivot_table(
        values='Precision', 
        index='Method', 
        columns='Pool', 
        aggfunc='first'
    )
    
    # Add mean and std columns
    summary_df['Mean_Precision'] = summary_df.mean(axis=1)
    summary_df['Std_Precision'] = summary_df.std(axis=1)
    
    # Calculate total TP, FP, and identifications across all pools
    totals_df = results_df.groupby('Method').agg({
        'True_Positives': 'sum',
        'False_Positives': 'sum',
        'Total_Identifications': 'sum'
    })
    
    # Calculate overall precision from totals
    totals_df['Overall_Precision'] = totals_df['True_Positives'] / totals_df['Total_Identifications']
    
    # Merge with summary
    summary_df = pd.concat([summary_df, totals_df[['Overall_Precision']]], axis=1)
    
    summary_file = 'precision/peptide_identification_precision_summary.csv'
    summary_df.to_csv(summary_file)
    print(f"Summary results saved to: {summary_file}")
    
    # Save totals separately
    totals_file = 'precision/peptide_identification_totals.csv'
    totals_df.to_csv(totals_file)
    print(f"Total counts saved to: {totals_file}")
    
    # Print summary
    print("\n" + "="*80)
    print("PRECISION SUMMARY")
    print("="*80)
    print(f"{'Method':<20} {'Mean±Std':<15} {'Overall':<10} {'Total IDs':<12} {'Total TP':<10} {'Total FP':<10}")
    print("-"*80)
    
    for method in summary_df.index:
        mean_prec = summary_df.loc[method, 'Mean_Precision']
        std_prec = summary_df.loc[method, 'Std_Precision']
        overall_prec = summary_df.loc[method, 'Overall_Precision']
        total_ids = totals_df.loc[method, 'Total_Identifications']
        total_tp = totals_df.loc[method, 'True_Positives']
        total_fp = totals_df.loc[method, 'False_Positives']
        
        print(f"{method:<20} {mean_prec:.3f}±{std_prec:.3f}   {overall_prec:.3f}      "
              f"{total_ids:<12} {total_tp:<10} {total_fp:<10}")
    
    # Find best method
    best_method = summary_df['Overall_Precision'].idxmax()
    best_precision = summary_df['Overall_Precision'].max()
    print(f"\nBest performing method (by overall precision): {best_method} ({best_precision:.4f})")
    
    return results_df, summary_df


if __name__ == "__main__":
    results_df, summary_df = main()