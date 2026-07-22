import os
import shutil
import pandas as pd
import argparse

def process_data(input_tsv, output_prefix, completeness_threshold, contamination_threshold, source_base_dir="."):
    """
    Reads the first three columns of a tab-separated file (skipping the header),
    checks completeness and contamination values against provided thresholds,
    and copies files (assuming .fa extension in source directory) to a
    directory derived from the output prefix, with sequential numbering
    appended to the output filename.

    Args:
        input_tsv (str): Path to the input TSV file.
        output_prefix (str): Prefix for the output files, including the desired
                             output directory and base name (e.g.,
                             "path/to/output/base_name_").
        completeness_threshold (float): Minimum completeness value to consider.
        contamination_threshold (float): Minimum contamination value to consider.
        source_base_dir (str, optional): Base directory where the .fa files
                                         are located.
                                         Defaults to the current working directory.
    """
    try:
        # Read only the first three columns, assuming the first row is the header
        df = pd.read_csv(input_tsv, sep='\t', header=0, usecols=[0, 9, 10],
                         names=['filename_base', 'completeness', 'contamination'])
    except FileNotFoundError:
        raise FileNotFoundError(f"Input file not found: {input_tsv}")
    except pd.errors.EmptyDataError:
        print(f"Warning: Input file is empty: {input_tsv}")
        return
    except pd.errors.ParserError as e:
        print(f"Error parsing input file: {e}")
        return

    # Convert completeness and contamination columns to numeric (coerce errors to NaN)
    df['completeness'] = pd.to_numeric(df['completeness'], errors='coerce')
    df['contamination'] = pd.to_numeric(df['contamination'], errors='coerce')

    # Extract output directory from the prefix
    output_dir = os.path.dirname(output_prefix)
    os.makedirs(output_dir, exist_ok=True)

    copied_count = 0
    for index, row in df.iterrows():
        # Get the base filename
        base_file_name = str(row['filename_base']).strip()
        full_source_name = f"{base_file_name}.fa"  # Add the .fa extension for source

        completeness = row['completeness']
        contamination = row['contamination']

        # Check for NaN values after conversion to numeric
        if pd.notna(completeness) and pd.notna(contamination):
            if completeness > completeness_threshold and \
               contamination < contamination_threshold:
                source_path = os.path.join(source_base_dir, full_source_name)
                output_filename = f"{os.path.basename(output_prefix)}{copied_count + 1}.fa"
                destination_path = os.path.join(output_dir, output_filename)

                if os.path.exists(source_path):
                    shutil.copy2(source_path, destination_path)
                    print(f"Copied: {full_source_name} from {source_base_dir} to "
                          f"{destination_path}")
                    copied_count += 1
                else:
                    print(f"Warning: File not found in source directory: "
                          f"{full_source_name} (looked in {source_base_dir})")
            else:
                print(f"Skipping {base_file_name} due to filter criteria.")
        else:
            print(f"Warning: Skipping {base_file_name} due to non-numeric "
                  f"completeness or contamination: {base_file_name}, "
                  f"{row['completeness']}, {row['contamination']}")

    print(f"\nTotal of {copied_count} files copied to: {output_dir}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Process a TSV file, filter based on completeness and "
                    "contamination, and copy .fa files with a combined output "
                    "prefix and sequential numbering.")
    parser.add_argument("input_tsv", help="Path to the input TSV file.")
    parser.add_argument("output_prefix",
                        help="Prefix for the output files, including the "
                             "desired output directory and base name (e.g., "
                             "'path/to/output/base_name_').")
    parser.add_argument("--completeness_threshold", type=float, default=50.0,
                        help="Minimum completeness threshold (default: 50.0).")
    parser.add_argument("--contamination_threshold", type=float, default=10.0,
                        help="Minimum contamination threshold (default: 10.0).")
    parser.add_argument("--source_dir", default=".",
                        help="Base directory for source files (default: "
                             "current directory).")

    args = parser.parse_args()

    process_data(args.input_tsv, args.output_prefix, args.completeness_threshold,
                 args.contamination_threshold, args.source_dir)
