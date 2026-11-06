#!/usr/bin/env python
"""
Tree species name translator
Usage: ./translate_tree.py dictionary.csv tree.tre -o output.tre
Note: may require installing some packages.

This script translates species names in Newick format phylogenetic trees 
using a translation dictionary provided in CSV format.

The dictionary file must be a comma-separated (CSV) file with at least two columns:
from: Original species names/IDs in the tree
to: New species names to replace the originals

"""

import sys
import subprocess

def check_packages():
    try:
        import pandas as pd
        import Bio
        print("All required packages are installed")
    except ImportError as e:
        missing = str(e).split("'")[1]
        print(f"Installing missing package: {missing}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", missing])

# Check packages first
check_packages()

# Import required packages
import pandas as pd
from Bio import Phylo
from pathlib import Path
import argparse

def clean_species_name(name: str) -> str:
    """Clean species name by removing quotes, extra spaces and normalize whitespace"""
    return ' '.join(name.strip().strip("'").strip('"').split())

def process_tree(tree_file: Path, dict_file: Path, output_file: Path = None) -> None:
    """Process single tree file using translation dictionary"""
    
    # Check if files exist
    if not tree_file.exists() or not dict_file.exists():
        print(f"Error: Input files not found")
        return
    
    try:
        # Read dictionary and clean species names
        df = pd.read_csv(dict_file)
        df['from'] = df['from'].apply(clean_species_name)
        df['to'] = df['to'].apply(clean_species_name)
        translation_dict = dict(zip(df['from'], df['to']))
        
        # Track translations
        translations_made = {}
        original_names = set()
        translated_names = set()
        
        # Read and process tree
        tree = Phylo.read(tree_file, 'newick')
        
        # Translate names
        for leaf in tree.get_terminals():
            orig_name = clean_species_name(leaf.name)
            original_names.add(orig_name)
            
            if orig_name in translation_dict:
                new_name = translation_dict[orig_name]
                translated_names.add(orig_name)
                leaf.name = new_name
                translations_made[orig_name] = new_name
        
        # Print translation report
        # print("\nTranslations made:")
        # for orig, new in translations_made.items():
        #     print(f"'{orig}' -> '{new}'")
        
        # Calculate statistics
        total_in_tree = len(original_names)
        translated_count = len(translated_names)
        untranslated_count = total_in_tree - translated_count
        unique_in_dict = len(df['from'].unique())
        
        # Print statistics
        print(f"\nStatistics for {tree_file.name}:")
        print(f"Species present in tree: {total_in_tree}")
        print(f"Names translated: {translated_count}")
        print(f"Names untranslated: {untranslated_count}")
        print(f"Total unique names in dictionary: {unique_in_dict}")
        
        # Save translated tree if output specified
        if output_file and translated_count > 0:
            Phylo.write(tree, output_file, 'newick', format_branch_length='%0.8f')
            print(f"\nTranslated tree saved to: {output_file}")
            
    except Exception as e:
        print(f"Error processing files: {e}")
        return

def main():
    parser = argparse.ArgumentParser(description='Translate species names in Newick trees')
    parser.add_argument('dict_file', type=Path, help='CSV dictionary file')
    parser.add_argument('tree_files', type=Path, nargs='+', help='Newick tree file(s)')
    parser.add_argument('-o', '--output', type=Path, help='Output file name (optional)')
    
    args = parser.parse_args()
    
    for tree_file in args.tree_files:
        if args.output and len(args.tree_files) > 1:
            output = tree_file.with_suffix('.translated.tre')
        else:
            output = args.output
        
        process_tree(tree_file, args.dict_file, output)

if __name__ == "__main__":
    main()