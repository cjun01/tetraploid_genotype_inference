# Tetraploid Dosage Inference Assumptions

This document describes the assumptions under which the current genotype-calling script is statistically valid **as written**.

## Purpose

The script infers the most likely dosage state at a **biallelic tetraploid locus** using observed reference (`REF`) and alternate (`ALT`) read counts. It classifies each locus into one of five possible dosage classes:

- `AAAA`
- `AAAB`
- `AABB`
- `ABBB`
- `BBBB`

The script is intended as a **simplified likelihood-based dosage classifier**, not as a full empirical model of sequencing bias, allele-specific distortion, or overdispersion.

## Core Assumptions

### 1. Biallelic tetraploid locus
Each locus is assumed to be:

- **biallelic**, with alleles `A` and `B`
- **tetraploid**, with exactly four allele copies

Therefore, the true genotype must belong to exactly one of the following five dosage classes:

\[
\text{AAAA},\ \text{AAAB},\ \text{AABB},\ \text{ABBB},\ \text{BBBB}
\]

---

### 2. Observed data are allele counts
For each sample at each locus, the observed data are:

- `a`: number of reads supporting the reference allele
- `b`: number of reads supporting the alternate allele
- `n = a + b`: total read depth

If `n = 0`, no genotype inference is made.

---

### 3. Expected reference allele fractions are determined by dosage
The expected reference-read fractions for the five dosage classes are assumed to be:

| Dosage | Genotype | Expected REF fraction |
|---|---|---:|
| 0 | AAAA | 1.00 |
| 1 | AAAB | 0.75 |
| 2 | AABB | 0.50 |
| 3 | ABBB | 0.25 |
| 4 | BBBB | 0.00 |

These values define the characteristic read-count pattern associated with each hidden dosage state.

---

### 4. Conditional independence of reads
Given the true dosage class, reads are assumed to arise independently with a common probability of supporting the reference allele.

Under this assumption, the reference read count `a` is modeled conditionally on total depth `n`.

---

### 5. Likelihood model for homozygous classes
For homozygous dosage classes, the script assumes a binomial model with a small fixed error term:

\[
e = 0.001
\]

Thus:

\[
a \mid \text{AAAA} \sim \text{Binomial}(n, 1 - e)
\]

\[
a \mid \text{BBBB} \sim \text{Binomial}(n, e)
\]

This allows a very small probability of observing the non-expected allele due to sequencing or calling error.

---

### 6. Likelihood model for heterozygous classes
For heterozygous dosage classes (`AAAB`, `AABB`, `ABBB`), the script assumes a **symmetric two-component mixture** around the expected reference fraction `f`.

For each heterozygous class:

\[
P(a \mid d) = \tfrac{1}{2}\,\text{Binomial}(n, f + e) + \tfrac{1}{2}\,\text{Binomial}(n, f - e)
\]

where:

- `f = 0.75` for `AAAB`
- `f = 0.50` for `AABB`
- `f = 0.25` for `ABBB`

This matches the implementation exactly and should be interpreted as a simplified symmetric deviation model around the ideal dosage-specific allele fraction.

---

### 7. Equal prior probability of dosage classes
Before observing read counts, all five dosage classes are assumed to be equally likely:

\[
P(\text{AAAA}) = P(\text{AAAB}) = P(\text{AABB}) = P(\text{ABBB}) = P(\text{BBBB})
\]

This assumption is important because the script converts class-specific likelihoods into probabilities by normalizing them across the five dosage states.

---

### 8. Posterior inference by normalized likelihoods
For each dosage class `d`, the script computes a class-specific likelihood:

\[
L_d = P(a \mid d)
\]

Under the equal-prior assumption, the posterior probability of dosage class `d` is obtained by:

\[
P(d \mid a, n) = \frac{L_d}{\sum_j L_j}
\]

Thus, the normalized values returned by the script can be interpreted as posterior probabilities under the model assumptions listed above.

---

### 9. Normal approximation at very high depth
For read depths up to the defined threshold (`n <= 1000`), the script uses the exact binomial probability.

For very high read depth (`n > 1000`), the script replaces the exact binomial calculation with a normal approximation having matching mean and variance.

Therefore:

- for `n <= 1000`, the inference is exact under the stated model
- for `n > 1000`, the inference is approximate under the stated model

---

### 10. Independent locus-wise inference
Each locus is inferred independently.

The script does **not**:

- borrow information across loci
- borrow information across samples
- model linkage disequilibrium
- model population structure
- estimate locus-specific bias parameters

---

## Interpretation of output probabilities

Under the assumptions above, the reported probabilities are valid **model-based posterior probabilities**.

For example, a call such as:

- `AABB (97.3%)`

means:

> Under the specified five-state dosage model, equal genotype priors, and the likelihood rules implemented in the script, the posterior probability of `AABB` is 0.973.

It does **not** necessarily mean that the probability is perfectly calibrated to all sources of real sequencing noise or biological bias.

---

## Confidence threshold

The script reports a genotype label only when the highest posterior probability is at least 95%.

So the rule is:

- if `max posterior >= 0.95`, report the best genotype
- otherwise, report `Low confidence`

This threshold is valid as a decision rule under the model assumptions above.

---

## Practical scope

This script is statistically defensible when used as a **simplified dosage classifier** for:

- biallelic tetraploid loci
- reference/alternate read count data
- locus-wise independent genotype inference
- a fixed small error term
- equal prior probabilities across the five dosage states

It should **not** be described as a fully calibrated empirical sequencing-error model unless additional validation is provided.

---

## Suggested methods text

### Short version

Genotype inference was performed for biallelic tetraploid loci by classifying each sample into one of five dosage states: `AAAA`, `AAAB`, `AABB`, `ABBB`, or `BBBB`. For each locus, reference and alternate read counts were modeled conditionally on dosage using binomial likelihoods with a small fixed error term (`e = 0.001`). Heterozygous classes were represented using a symmetric two-component mixture around their expected allele fractions. Equal prior probabilities were assumed across dosage classes, and posterior probabilities were obtained by normalizing class-specific likelihoods. For depths greater than 1000, a normal approximation to the binomial likelihood was used.

### Full version

For each biallelic tetraploid locus, genotype inference was treated as classification among five possible dosage states: `AAAA`, `AAAB`, `AABB`, `ABBB`, and `BBBB`. For each sample at a locus, the observed data consisted of the number of reference-supporting reads (`a`) and alternate-supporting reads (`b`), with total depth defined as `n = a + b`. Expected reference allele fractions were assumed to be `1.0`, `0.75`, `0.5`, `0.25`, and `0.0` for the five dosage states, respectively. Read counts were modeled conditionally on dosage using binomial likelihoods. Homozygous classes were modeled using a small fixed error term (`e = 0.001`), such that `AAAA` and `BBBB` corresponded to reference-read probabilities of `1 - e` and `e`, respectively. Heterozygous classes were modeled using a symmetric two-component mixture around the expected reference fraction, evaluated at `f - e` and `f + e`. Equal prior probabilities were assumed for all five dosage classes, and posterior probabilities were obtained by normalizing class-specific likelihoods across the five states. For depths greater than 1000, a normal approximation to the binomial likelihood was used.

---

## Limitation statement

A precise way to describe the method is:

> This script implements a simplified Bayesian-style likelihood framework for dosage inference in biallelic tetraploids, with equal genotype priors and a fixed small error term.

A cautious limitation statement is:

> The method is intended as a simplified model-based classifier rather than a full empirical model of allele-specific bias, overdispersion, or locus-specific sequencing error.

