import pandas as pd
import argparse
from math import comb, sqrt, exp, pi
from tqdm import tqdm

# Function to extract the reference and alternate alleles from the ref/alt format
def extract_ref_alt(value):
    # Check if the value is NaN or not a string
    if pd.isna(value) or not isinstance(value, str):
        return 0, 0  # You can decide how to handle missing values. Here, I treat it as 0 ref and 0 alt.
    
    ref, alt = map(int, value.split('/'))
    return ref, alt

# Normal approximation to the binomial distribution for large n
def normal_approximation(n, k, p):
    mean = n * p
    variance = n * p * (1 - p)
    std_dev = sqrt(variance)
    
    # We use the probability density function (PDF) of the normal distribution
    if std_dev == 0:
        return 0  # Avoid division by zero

    # Normal distribution formula for PDF
    return (1 / (std_dev * sqrt(2 * pi))) * exp(-0.5 * ((k - mean) / std_dev) ** 2)
  # Mapping of dosage to genotype labels
dosage_to_genotype = {
    0: 'AAAA',
    1: 'AAAB',
    2: 'AABB',
    3: 'ABBB',
    4: 'BBBB'
}
# Calculate probabilities for each row and select the most probable genotype if over 95%
def calculate_probabilities_for_row(row, columns):
    inferred_genotypes = {}
    for col in columns:
        ref, alt = extract_ref_alt(row[col])
        probabilities = calculate_probabilities(ref, alt)

        if probabilities:
            best_dosage = max(probabilities, key=probabilities.get)
            best_genotype = dosage_to_genotype[best_dosage]
            best_prob = probabilities[best_dosage] * 100
            if best_prob >= 95:
                inferred_genotypes[col] = f"{best_genotype} ({best_prob:.2f}%)"
            else:
                inferred_genotypes[col] = "Low confidence"
        else:
            inferred_genotypes[col] = "No reads"

    return inferred_genotypes
def calculate_probabilities(a, b, threshold=1000):
    e = 0.001  # Error factor
    n = a + b
    if n == 0:
        return {}

    def _clip01(p):
        tiny = 1e-15
        if p <= 0.0:
            return tiny
        if p >= 1.0:
            return 1.0 - tiny
        return p

    # Use binomial probability for small n, normal approximation for large n
    def binomial_probability(n, k, p):
        p = _clip01(p)
        if n > threshold:
            return normal_approximation(n, k, p)
        else:
            return comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

    # Dosage -> ideal REF fraction f
    f_values = {
        0: 1.0,        # AAAA
        1: 3 / 4,      # AAAB
        2: 2 / 4,      # AABB
        3: 1 / 4,      # ABBB
        4: 0.0         # BBBB
    }

    probs = {}
    for dosage, f in f_values.items():
        if dosage == 0:
            # AAAA
            probs[dosage] = binomial_probability(n, a, 1.0 - e)
        elif dosage == 4:
            # BBBB
            probs[dosage] = binomial_probability(n, a, 0.0 + e)
        else:
            # Symmetric ±e mixture around f (no directional bias)
            p_plus = _clip01(f + e)
            p_minus = _clip01(f - e)
            probs[dosage] = 0.5 * binomial_probability(n, a, p_plus) + 0.5 * binomial_probability(n, a, p_minus)

    # Normalize probabilities
    total_prob = sum(probs.values())
    return {k: v / total_prob for k, v in probs.items()} if total_prob != 0 else {}

# Main function to handle command-line arguments and file processing
def main():
    # Argument parser for input and output files
    parser = argparse.ArgumentParser(description='Infer genotypes and probabilities from a CSV file.')
    parser.add_argument('-i', '--input', required=True, type=str, help='Path to the input CSV file.')
    parser.add_argument('-o', '--output', required=True, type=str, help='Path to the output CSV file.')

    args = parser.parse_args()

    # Read the input CSV file
    df = pd.read_csv(args.input)

    # Extract the columns containing genotype information
    columns = df.columns[2:]  # Assuming the first two columns are 'Chromosome' and 'Position'

    # Apply the inference function to each row with tqdm for progress bar
    inferred_data = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing rows"):
        inferred_data.append(calculate_probabilities_for_row(row, columns))

    # Create a DataFrame with the inferred genotypes
    inferred_df = pd.DataFrame(inferred_data, index=df.index)
    
    # Merge with the original 'Chromosome' and 'Position' columns
    output_df = pd.concat([df[['Chromosome', 'Position']], inferred_df], axis=1)
    output_df = output_df.drop_duplicates(subset=['Chromosome', 'Position'], keep='first')
    # Save the output to a CSV file
    output_df.to_csv(args.output, index=False)

    print(f"Output saved to {args.output}")

if __name__ == "__main__":
    main()
