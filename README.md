# Normalized Effort-to-Compress for MATLAB

MATLAB implementation of Effort-to-Compress (ETC), nETC, and cETC. 

## Usage

Add this folder to the MATLAB path:

```matlab
addpath('/path/to/normalized-effort-to-compress-matlab')

compute_etc([0 1 0 1])   % 1: one substitution gives [2 2]
compute_etc([7 7 7])     % 0: already constant
compute_etc([0 1 2])     % 2
```

`compute_etc` accepts row or column vectors of finite integer labels, including negative labels. It returns the number of pair substitutions needed to reach a constant sequence. The API returns zero for empty input.

## Normalization

```matlab
sequence = [0 0 0 1 0 1 0 0 1 0];
n = nETC(sequence);
c = cETC(sequence);

% Optional: 1,000 reference samples and a reproducible seed.
[n, global_info] = nETC(sequence, 1000, 42);
[c, composition_info] = cETC(sequence, 1000, 42);
```

Each function computes raw ETC, its certified lower bound, and a Monte Carlo reference mean. The result is `(ETC - lower_bound) / (reference_mean - lower_bound)`. The default is 200 reference samples. With a seed, the functions use a local random stream and leave MATLAB's global random state unchanged. Without a seed, they use the global stream.

`nETC` draws uniformly from sequences with the same length and observed alphabet size, requiring every symbol to appear. It uses rejection sampling when acceptance is high and a completion-count sampler otherwise. `cETC` uniformly shuffles the observed sequence, preserving the count associated with each numerical label.

For length `L` and alphabet size `k`, the global lower bound is `k - 1 + ceil(log2(L/(k*M)))`, where `M` is the largest divisor of `L` not exceeding `floor(L/k)`. The composition-conditioned bound replaces `M` with the greatest common divisor of the observed symbol counts. Both are evaluated by integer doubling to avoid rounding at powers of two.

The optional second output includes `raw_etc`, `lower_bound`, `reference_mean`, `reference_sd`, `reference_mcse`, the individual `reference_etc` values, and the ordered `symbols` and `counts`.

Normalization requires `2 <= k <= L-2`. Constant sequences and the boundaries `k = L-1` and `k = L` have zero reference denominators. If a finite sample also gives a zero denominator, the function reports an error so the sample count can be increased. Values above one are retained.

## Substitution rules

Pairs are counted without overlap in a left-to-right scan. Ties are resolved by the smallest second label, then the smallest first label. Each substitution uses a fresh label larger than all existing labels. The algorithm checks for an exactly constant sequence before counting pairs. When the maximum pair frequency is one, it adds the current length minus one to the count and returns.

Original labels are mapped to ascending ranks. Their numerical order and the count associated with each label are preserved. Pair counts use a sparse matrix.

## References

- N. Nagaraj, K. Balasubramanian, and S. Dey, A new complexity measure for time series analysis and classification, European Physical Journal Special Topics 222, 847-860 (2013).
- K. Balasubramanian, [Effort-To-Compress (ETC) measure using Non-Sequential Recursive Pair Substitution (NSRPS)](https://www.mathworks.com/matlabcentral/fileexchange/67665-effort-to-compress-etc-measure-using-non-sequential-recursive-pair-substitution-nsrps), MATLAB Central File Exchange (2018).

MIT license. This repository implements the same substitution convention as the companion C repository. The earlier File Exchange implementation is available from the reference above.
