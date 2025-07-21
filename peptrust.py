import pandas as pd
import numpy as np
import os
from pathlib import Path
import sys

class PepTrust:
    """
    PepTrust method for peptide identification validation
    """
    
    def __init__(self, num_datasets=5, dataset_prefix='pool'):
        """
        Initialize PepTrust analysis
        
        Args:
            num_datasets: Number of datasets to analyze (e.g., 5 for pool_1 to pool_5)
            dataset_prefix: Prefix for dataset names (e.g., 'pool' for pool_1, pool_2, etc.)
        """
        self.num_datasets = num_datasets
        self.dataset_prefix = dataset_prefix
        self.search_engines = ['comet', 'msgf', 'tandem']
        self.combination_methods = ['geometric_mean', 'msblender', 'voting']
        
        # Create output directory
        self.output_dir = Path('peptrust_results')
        self.output_dir.mkdir(exist_ok=True)
        
    def get_dataset_name(self, dataset_num):
        """Generate dataset name based on prefix and number"""
        return f"{self.dataset_prefix}_{dataset_num}"
    
    def get_file_path(self, dataset_num, file_type, method):
        """
        Get the correct file path for a dataset
        
        Args:
            dataset_num: Dataset number (1, 2, 3, etc.)
            file_type: 'search' or 'combination'
            method: Method name (e.g., 'comet', 'geometric_mean')
        """
        if file_type == 'search':
            return f"searched_results/{method}/{self.dataset_prefix}{dataset_num}.csv"
        else:  # combination
            if method == 'geometric_mean':
                return f"combined_results/{method}/{self.dataset_prefix}{dataset_num}_combine.csv"
            elif method == 'msblender':
                return f"combined_results/{method}/{self.dataset_prefix}{dataset_num}.csv"
            else:  # voting
                return f"combined_results/{method}/{self.dataset_prefix}_{dataset_num}.csv"
    
    def load_dataset(self, dataset_num):
        """
        Load all search and combination results for a specific dataset
        
        Returns:
            dict: Contains 'search' and 'combination' DataFrames, or None if loading fails
        """
        dataset_name = self.get_dataset_name(dataset_num)
        print(f"Loading {dataset_name}...", end='')
        
        data = {
            'search': {},
            'combination': {}
        }
        
        # Load search engine results
        for engine in self.search_engines:
            filepath = self.get_file_path(dataset_num, 'search', engine)
            try:
                df = pd.read_csv(filepath)
                data['search'][engine] = df[['scan', 'Peptide']]
            except Exception as e:
                print(f"\n  Error loading {filepath}: {e}")
                return None
        
        # Load combination method results
        for method in self.combination_methods:
            filepath = self.get_file_path(dataset_num, 'combination', method)
            try:
                df = pd.read_csv(filepath)
                
                if method == 'geometric_mean':
                    # Extract scan and peptide from psm column
                    df['scan'] = df['psm'].str.split('-').str[0].astype(int)
                    df['Peptide'] = df['psm'].str.split('-').str[1]
                    data['combination'][method] = df[['scan', 'Peptide']]
                else:
                    # For voting and msblender
                    data['combination'][method] = df[['scan', 'plain_peptide']].rename(
                        columns={'plain_peptide': 'Peptide'}
                    )
            except Exception as e:
                print(f"\n  Error loading {filepath}: {e}")
                return None
        
        print(" Done!")
        return data
    
    def calculate_trustworthiness(self, data):
        """
        Calculate τ_{Si|Mj} for each search engine under each combination method
        
        Returns:
            dict: Trustworthiness values
        """
        trustworthiness = {}
        
        for comb_method, comb_df in data['combination'].items():
            for search_engine, search_df in data['search'].items():
                # Merge search and combination results on scan
                merged = search_df.merge(comb_df, on='scan', suffixes=('_search', '_comb'))
                
                # Count matches where peptides agree
                matches = (merged['Peptide_search'] == merged['Peptide_comb']).sum()
                total = len(merged)
                
                # Calculate trustworthiness
                tau = matches / total if total > 0 else 0.5
                trustworthiness[(search_engine, comb_method)] = tau
        
        return trustworthiness
    
    def calculate_log_likelihood(self, data, trustworthiness):
        """
        Calculate log-likelihood ℒ(Mj) for each combination method
        
        Returns:
            dict: Log-likelihood values
        """
        log_likelihoods = {}
        
        for comb_method, comb_df in data['combination'].items():
            log_likelihood = 0.0
            
            for search_engine, search_df in data['search'].items():
                tau = trustworthiness[(search_engine, comb_method)]
                
                # Merge search and combination results
                merged = search_df.merge(comb_df, on='scan', suffixes=('_search', '_comb'))
                
                # Calculate log-likelihood contribution
                for _, row in merged.iterrows():
                    if row['Peptide_search'] == row['Peptide_comb']:
                        # Match case: contribute log(τ)
                        if tau > 0:
                            log_likelihood += np.log(tau)
                    else:
                        # No match case: contribute log(1-τ)
                        if (1 - tau) > 0:
                            log_likelihood += np.log(1 - tau)
            
            log_likelihoods[comb_method] = log_likelihood
        
        return log_likelihoods
    
    def save_dataset_results(self, dataset_num, trustworthiness, log_likelihoods):
        """Save results for a single dataset"""
        dataset_name = self.get_dataset_name(dataset_num)
        
        # Save trustworthiness for this dataset
        trust_rows = []
        for (engine, method), tau in trustworthiness.items():
            trust_rows.append({
                'dataset': dataset_name,
                'search_engine': engine,
                'combination_method': method,
                'trustworthiness': tau
            })
        
        trust_df = pd.DataFrame(trust_rows)
        trust_file = self.output_dir / f'{dataset_name}_trustworthiness.csv'
        trust_df.to_csv(trust_file, index=False)
        
        # Save log-likelihood for this dataset
        ll_rows = []
        for method, ll in log_likelihoods.items():
            ll_rows.append({
                'dataset': dataset_name,
                'combination_method': method,
                'log_likelihood': ll
            })
        
        ll_df = pd.DataFrame(ll_rows)
        ll_file = self.output_dir / f'{dataset_name}_log_likelihood.csv'
        ll_df.to_csv(ll_file, index=False)
        
        return trust_df, ll_df
    
    def run_analysis(self):
        """
        Run the complete PepTrust analysis on all datasets
        """
        print(f"\nPepTrust Analysis")
        print(f"Analyzing {self.num_datasets} datasets with prefix '{self.dataset_prefix}'")
        print("="*70)
        
        # Store all individual results
        all_trust_dfs = []
        all_ll_dfs = []
        
        # Process each dataset
        for dataset_num in range(1, self.num_datasets + 1):
            # Load dataset
            data = self.load_dataset(dataset_num)
            if data is None:
                print(f"  Skipping {self.get_dataset_name(dataset_num)} due to loading errors")
                continue
            
            # Calculate trustworthiness
            trustworthiness = self.calculate_trustworthiness(data)
            
            # Calculate log-likelihood
            log_likelihoods = self.calculate_log_likelihood(data, trustworthiness)
            
            # Save individual dataset results
            trust_df, ll_df = self.save_dataset_results(dataset_num, trustworthiness, log_likelihoods)
            all_trust_dfs.append(trust_df)
            all_ll_dfs.append(ll_df)
            
            print(f"  Saved results for {self.get_dataset_name(dataset_num)}")
        
        # Calculate and save summary statistics
        if all_trust_dfs and all_ll_dfs:
            self.save_summary_results(all_trust_dfs, all_ll_dfs)
        
        print("\nAnalysis complete!")
        print(f"Results saved in '{self.output_dir}/' directory")
    
    def save_summary_results(self, all_trust_dfs, all_ll_dfs):
        """Calculate and save summary statistics across all datasets"""
        print("\nCalculating summary statistics...")
        
        # Combine all trustworthiness results
        all_trust = pd.concat(all_trust_dfs, ignore_index=True)
        
        # Calculate mean trustworthiness
        mean_trust = all_trust.groupby(['search_engine', 'combination_method'])['trustworthiness'].agg([
            'mean', 'std', 'count'
        ]).reset_index()
        mean_trust.columns = ['search_engine', 'combination_method', 'mean_trustworthiness', 'std_trustworthiness', 'n_datasets']
        mean_trust.to_csv(self.output_dir / 'summary_trustworthiness.csv', index=False)
        
        # Combine all log-likelihood results
        all_ll = pd.concat(all_ll_dfs, ignore_index=True)
        
        # Calculate total log-likelihood
        total_ll = all_ll.groupby('combination_method')['log_likelihood'].agg([
            'sum', 'mean', 'std', 'count'
        ]).reset_index()
        total_ll.columns = ['combination_method', 'total_log_likelihood', 'mean_log_likelihood', 'std_log_likelihood', 'n_datasets']
        total_ll = total_ll.sort_values('total_log_likelihood', ascending=False)
        total_ll.to_csv(self.output_dir / 'summary_log_likelihood.csv', index=False)
        
        # Find optimal method
        optimal_method = total_ll.iloc[0]['combination_method']
        optimal_ll = total_ll.iloc[0]['total_log_likelihood']
        
        # Save optimal method
        with open(self.output_dir / 'optimal_method.txt', 'w') as f:
            f.write("PepTrust Analysis Results\n")
            f.write("="*50 + "\n\n")
            f.write(f"Optimal Combination Method: {optimal_method}\n")
            f.write(f"Total Log-Likelihood: {optimal_ll:.4f}\n")
            f.write(f"Number of Datasets: {len(all_trust_dfs)}\n")
            f.write(f"Dataset Prefix: {self.dataset_prefix}\n\n")
            f.write("This combination method should be used as the ground truth.\n")
        
        # Display summary
        print("\n" + "="*70)
        print("SUMMARY RESULTS")
        print("="*70)
        
        # Display mean trustworthiness matrix
        print("\nMean Source Trustworthiness τ_{Si|Mj}:")
        print("-"*70)
        print(f"{'Search Engine':<15}", end='')
        for method in self.combination_methods:
            print(f"{method:<20}", end='')
        print()
        print("-"*70)
        
        for engine in self.search_engines:
            print(f"{engine:<15}", end='')
            for method in self.combination_methods:
                trust_val = mean_trust[
                    (mean_trust['search_engine'] == engine) & 
                    (mean_trust['combination_method'] == method)
                ]['mean_trustworthiness'].values
                if len(trust_val) > 0:
                    print(f"{trust_val[0]:<20.4f}", end='')
                else:
                    print(f"{'N/A':<20}", end='')
            print()
        
        # Display total log-likelihoods
        print("\n\nTotal Log-Likelihoods ℒ(Mj):")
        print("-"*50)
        print(f"{'Combination Method':<25} {'Total LL':<20} {'Mean LL':<20}")
        print("-"*50)
        
        for _, row in total_ll.iterrows():
            print(f"{row['combination_method']:<25} {row['total_log_likelihood']:<20.4f} {row['mean_log_likelihood']:<20.4f}")
        
        print("\n" + "="*70)
        print(f"OPTIMAL METHOD: {optimal_method}")
        print(f"Total Log-Likelihood: {optimal_ll:.4f}")
        print("="*70)


def main():
    """Main function with command line argument handling"""
    # Default values
    num_datasets = 5
    dataset_prefix = 'pool'
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        try:
            num_datasets = int(sys.argv[1])
        except ValueError:
            print("Error: First argument must be a number")
            print("Usage: python peptrust_analysis.py [num_datasets] [dataset_prefix]")
            print("Example: python peptrust_analysis.py 10 pool")
            sys.exit(1)
    
    if len(sys.argv) > 2:
        dataset_prefix = sys.argv[2]
    
    # Run analysis
    analyzer = PepTrust(num_datasets=num_datasets, dataset_prefix=dataset_prefix)
    analyzer.run_analysis()


if __name__ == "__main__":
    main()