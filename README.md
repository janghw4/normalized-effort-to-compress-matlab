# Normalized Effort-to-Compress for MATLAB

Pure MATLAB implementation of Effort-to-Compress (ETC). No MEX compilation or additional toolboxes are required.

## Usage

Add this folder to the MATLAB path:

```matlab
addpath('/path/to/normalized-effort-to-compress-matlab')

compute_etc([0 1 0 1])   % 1: one substitution gives [2 2]
compute_etc([7 7 7])     % 0: already constant
compute_etc([0 1 2])     % 2
example_etc
```

`compute_etc` accepts row or column vectors of finite integer labels, including negative labels. It returns the number of pair substitutions needed to reach a constant sequence. The API returns zero for empty input.

## Normalization

```matlab
value = normalize_etc(etc, lower_bound, reference_mean);
```

This computes `(etc - lower_bound) / (reference_mean - lower_bound)`. Supply the certified global lower bound and matching global ensemble mean for nETC. For cETC, supply the corresponding anchors conditioned on the ordered symbol-count vector. Reference means must be estimated separately.

The reference mean must exceed the lower bound. Values above one are retained.

## Substitution rules

Pairs are counted without overlap in a left-to-right scan. Ties are resolved by the smallest second label, then the smallest first label. Each substitution uses a fresh label larger than all existing labels. The algorithm checks for an exactly constant sequence before counting pairs. When the maximum pair frequency is one, it adds the current length minus one to the count and returns.

Original labels are mapped to ascending ranks. Their numerical order and the count associated with each label are preserved. Pair counts use a sparse matrix.

## Tests

```matlab
addpath('tests')
run_tests
```

Tests compare with an explicit pair-by-pair MATLAB implementation and saved results from the [C implementation](https://github.com/janghw4/normalized-effort-to-compress). The tests cover all binary sequences through length 12, all ternary sequences through length 8, random sequences, sparse labels, row and column inputs, and invalid inputs.

Verified on MATLAB R2024b Update 6: 18,255 oracle comparisons and 120 C comparisons passed.

The C fixture can be regenerated with Python and a compiled C library:

```sh
python tests/generate_c_cases.py /path/to/libetc.so
```

## References

- N. Nagaraj, K. Balasubramanian, and S. Dey, A new complexity measure for time series analysis and classification, European Physical Journal Special Topics 222, 847-860 (2013).
- K. Balasubramanian, [Effort-To-Compress (ETC) measure using Non-Sequential Recursive Pair Substitution (NSRPS)](https://www.mathworks.com/matlabcentral/fileexchange/67665-effort-to-compress-etc-measure-using-non-sequential-recursive-pair-substitution-nsrps), MATLAB Central File Exchange (2018).

MIT license. This repository implements the same substitution convention as the companion C repository. The earlier File Exchange implementation is available from the reference above.
