import pandas as pd
import numpy as np
import argparse

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

# this function generates a header for the output file
def generate_header():
    header = ["original", "variant", "change_description", "y_pred"]
    # then we join the header list with tabs
    print("\t".join(header))
        
def generate_mismatches(spacers, min_score, max_score, step, parameters):
    nucleotides = ["A", "C", "G", "T"]
    params = parameters["weight"].to_dict()

    for spacer in spacers:
        desired_scores = np.arange(min_score, max_score + step, step)

        mismatches = []
        for pos in range(len(spacer)):
            for nt in nucleotides:
                if nt == spacer[pos]:
                    continue
                target_mismatch = spacer[:pos] + nt + spacer[pos+1:]
                y_pred = calculate_y_pred(spacer, target_mismatch, params['GC_content'], params)
                mismatches.append(((pos, nt), y_pred))

        mismatch_list = []
        for score in desired_scores:
            closest_score = None
            closest_mismatch = None
            for mismatch, mismatch_score in mismatches:
                if closest_score is None or abs(mismatch_score - score) < abs(closest_score - score):
                    if mismatch not in [x[0] for x in mismatch_list]:
                        closest_score = mismatch_score
                        closest_mismatch = mismatch
            mismatch_list.append((closest_mismatch, closest_score))

        for i, (mismatch, score) in enumerate(mismatch_list):
            target_mismatch = spacer[:mismatch[0]] + mismatch[1] + spacer[mismatch[0]+1:]
            change_description = f"{spacer[mismatch[0]]}{mismatch[0]+1}{mismatch[1]}"
            print(f"{spacer}\t{target_mismatch}\t{change_description}\t{score:.4f}")
        print()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate mismatches for a list of spacers.")
    parser.add_argument("--spacers_file", required=True, help="Path to the file containing original spacers (one per line).")
    parser.add_argument("--parameters", required=True, help="Path to the parameters file (CSV format).")
    parser.add_argument("--min", type=float, default=0, help="Minimum desired efficacy (default: 0).")
    parser.add_argument("--max", type=float, default=1, help="Maximum desired efficacy (default: 1).")
    parser.add_argument("--step", type=float, default=0.1, help="Step between desired efficacies (default: 0.1).")
    args = parser.parse_args()

    with open(args.spacers_file, "r") as f:
        spacers = [line.strip() for line in f.readlines() if line.strip()]

    parameters = pd.read_csv(args.parameters, index_col="feature")
    generate_header()
    generate_mismatches(spacers, args.min, args.max, args.step, parameters)

