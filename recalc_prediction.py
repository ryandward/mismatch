import pandas as pd
import argparse
import logging
import sys

def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

def calculate_y_pred(original, variant, gc_weight, params):
    y_pred = params['intercept']

    for pos, (orig_base, var_base) in enumerate(zip(original, variant)):
        if orig_base != var_base:
            y_pred += params[f'{pos}']
            y_pred += params[f'{orig_base}{var_base}']

    y_pred += gc_weight * gc_content(original)
    return y_pred

def read_parameters(file_path):
    try:
        df = pd.read_csv(file_path)
        params = df.set_index("feature")["weight"].to_dict()
        return params
    except Exception as e:
        logging.error(f"Error reading parameters file: {e}")
        sys.exit(1)

def main(args):
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=getattr(logging, args.verbosity.upper()))

    logging.info(f"Reading parameters from {args.parameters_file}...")
    params = read_parameters(args.parameters_file)

    logging.info(f"Reading input data from {args.input_data_file}...")
    try:
        data = pd.read_csv(args.input_data_file, sep="\t")
    except Exception as e:
        logging.error(f"Error reading input data file: {e}")
        sys.exit(1)

    logging.info("Calculating y_pred for each row...")
    data['y_pred_calc'] = data.apply(lambda row: calculate_y_pred(row['original'], row['variant'], params['GC_content'], params), axis=1)

    logging.info("Displaying results:")
    print(data[['original', 'variant', 'y_pred', 'y_pred_calc']])

    logging.info(f"Saving results to {args.output_file}...")
    try:
        data.to_csv(args.output_file, sep="\t", index=False)
    except Exception as e:
        logging.error(f"Error saving output file: {e}")
        sys.exit(1)

    logging.info("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate CRISPRi knockdown efficacy with mismatched nucleotides.")
    parser.add_argument("input_data_file", help="Path to the input data file (TSV format).")
    parser.add_argument("parameters_file", help="Path to the parameters file (CSV format).")
    parser.add_argument("output_file", help="Path to the output file (TSV format).")
    parser.add_argument("--verbosity", choices=["debug", "info", "warning", "error", "critical"], default="info", help="Set the logging verbosity level (default: info).")
    args = parser.parse_args()

    main(args)
