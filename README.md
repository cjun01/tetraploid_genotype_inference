# Tetraploid Genotype Inference from REF/ALT Read Counts

This script infers tetraploid genotype dosage from REF/ALT read counts using a modified naïve binomial likelihood model. It is designed for autotetraploid organisms such as potato, where each biallelic SNP can have five possible dosage classes.

| ALT dosage | Genotype | Expected REF fraction |
|---:|---|---:|
| 0 | AAAA | 1.00 |
| 1 | AAAB | 0.75 |
| 2 | AABB | 0.50 |
| 3 | ABBB | 0.25 |
| 4 | BBBB | 0.00 |

The approach is adapted from the naïve binomial dosage-probability framework described by Yamamoto et al. (2020), with a small error term used to avoid boundary probabilities and to provide limited flexibility around expected allele-balance ratios.

---

## Main features

- Infers tetraploid dosage classes from REF/ALT read counts.
- Uses a small error term, default `error = 0.001`.
- Calculates likelihoods on the log scale to avoid numerical underflow.
- Reports genotype calls only when the maximum normalized likelihood exceeds a user-defined confidence threshold.
- Marks uncertain calls as `Low confidence`.
- Marks missing or malformed count values as `No reads`.

---

## Model description

For each SNP and individual, the script calculates genotype likelihoods for the five tetraploid dosage classes.

For homozygous classes, the expected REF probabilities are error-adjusted:

```text
AAAA: p_REF = 1 - error
BBBB: p_REF = error
```

For heterozygous classes, the script uses a symmetric two-component mixture around the expected REF fraction:

```text
AAAB: p_REF = 0.75 ± error
AABB: p_REF = 0.50 ± error
ABBB: p_REF = 0.25 ± error
```

For example, when `error = 0.001`, the AAAB class is evaluated as:

```text
0.5 × Binomial(p = 0.751) + 0.5 × Binomial(p = 0.749)
```

All likelihoods are calculated using log-binomial probabilities:

```text
log P(k | n, p)
```

where:

```text
k = number of REF reads
n = total read depth = REF + ALT
p = expected REF fraction for a dosage class
```

The likelihoods are then normalized across the five dosage classes:

```text
normalized likelihood = likelihood for one dosage / sum of likelihoods for all dosages
```

The genotype with the highest normalized likelihood is selected. If the maximum normalized likelihood is below the confidence threshold, the call is classified as `Low confidence`.

---

## Why log-binomial likelihood is used

The current version uses exact binomial log-likelihoods for all depths. This avoids numerical underflow while preserving accurate likelihood calculations for both low-depth and high-depth sites, including extreme REF/ALT ratios.

---

## Input format

The input file should be a CSV file.

The first two columns are treated as metadata, usually:

```text
Chromosome, Position
```

All columns from the third column onward are treated as sample genotype-count columns.

Each sample cell should contain REF and ALT read counts in this format:

```text
REF/ALT
```

Example:

```csv
Chromosome,Position,Sample1,Sample2,Sample3
Chr01,10025,39/0,30/10,0/42
Chr01,10580,741/39,20/20,3/55
Chr02,20210,.,15/5,NA
```

Accepted missing values include:

```text
.
NA
NaN
nan
./.
empty cell
```

Malformed or missing values are reported as:

```text
No reads
```

---

## Output format

The output file is a CSV file containing the first two metadata columns plus inferred genotype calls for each sample.

Example:

```csv
Chromosome,Position,Sample1,Sample2,Sample3
Chr01,10025,AAAA (99.99%),AAAB (98.50%),BBBB (99.99%)
Chr01,10580,AAAB (96.15%),AABB (99.20%),BBBB (97.80%)
Chr02,20210,No reads,AAAB (97.60%),Low confidence
```

Each confident call is reported as:

```text
Genotype (probability%)
```

For example:

```text
AAAB (96.15%)
```

If the maximum normalized likelihood is below the confidence threshold, the output is:

```text
Low confidence
```

If no valid REF/ALT read counts are available, the output is:

```text
No reads
```

---

## Usage

Basic usage:

```bash
python genotype_inference_no_normal.py \
  -i allele_count_Ref_to_Alt.csv \
  -o inferred_genotypes.csv
```

With custom error and confidence threshold:

```bash
python genotype_inference_no_normal.py \
  -i allele_count_Ref_to_Alt.csv \
  -o inferred_genotypes.csv \
  --error 0.001 \
  --confidence 0.95
```

Example on an HPC interactive node:

```bash
salloc --time=3:00:00 --mem=32G --cpus-per-task=2
srun --pty bash

conda activate potato_env

cd /home/cflzxc/potato_smart_culling/Python

python genotype_inference_no_normal.py \
  -i ../allele_count_Ref_to_Alt.csv \
  -o inferred_genotypes.csv
```

---

## Arguments

| Argument | Required | Default | Description |
|---|---:|---:|---|
| `-i`, `--input` | yes | none | Input CSV file containing REF/ALT counts |
| `-o`, `--output` | yes | none | Output CSV file for inferred genotypes |
| `--error` | no | `0.001` | Small error term used in genotype likelihoods |
| `--confidence` | no | `0.95` | Minimum normalized likelihood required to report a genotype |

---

## Example interpretation

For a count such as:

```text
741/39
```

The REF fraction is:

```text
741 / (741 + 39) = 0.95
```

With `error = 0.001`, the model is strict for homozygous classes. Although 95% REF may visually appear close to AAAA, the model expects AAAA to be near 99.9% REF. Therefore, this count may be classified as AAAB if the AAAB likelihood is higher than the AAAA likelihood.

This reflects a key feature of the simplified binomial model: low-frequency minor reads may be interpreted as evidence for heterozygous dosage rather than sequencing, mapping, or allele-balance noise.

---

## Important notes and limitations

This script is a simplified dosage caller. It does not explicitly model:

- allele bias,
- overdispersion,
- SNP-specific sequencing error,
- mapping artifacts,
- paralogous alignment,
- population-level genotype priors.

Tools such as `updog` can model some of these effects and may classify low-minor-allele cases differently.

Therefore, this script is best used as:

- a fast tetraploid dosage caller,
- a transparent binomial baseline,
- an independent comparison or validation method,
- a QC tool for evaluating genotype-calling robustness.

Discordant calls between this script and more complex methods such as `updog` should not automatically be interpreted as errors in either method unless independent validation data are available.

---

## Suggested Materials and Methods wording

A possible description for a manuscript is:

> We adapted the naïve binomial probability approach of Yamamoto et al. (2020) to infer tetraploid genotype dosage from REF/ALT read counts extracted from the VCF. For each SNP and individual, likelihoods were calculated across the five possible dosage classes using dosage-specific expected allele fractions. A small error term (e = 0.001) was included to avoid boundary probabilities and to allow limited deviation from ideal allele balance. The dosage class with the highest normalized likelihood was retained, whereas calls below the confidence threshold were classified as low confidence.

---

## Reference

Yamamoto E, Shirasawa K, Kimura T, Monden Y, Tanaka M, Isobe S. 2020. Genetic mapping in autohexaploid sweet potato with low-coverage NGS-based genotyping data. *G3: Genes, Genomes, Genetics* 10(8):2661–2670. doi:10.1534/g3.120.401433
