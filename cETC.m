function [value, info] = cETC(sequence, numSamples, seed)
%CETC Composition-conditioned ETC, including both reference anchors.
%   VALUE = CETC(SEQUENCE) uses 200 shuffled references with the same
%   ordered symbol counts. [VALUE, INFO] = CETC(SEQUENCE, NUMSAMPLES, SEED)
%   sets the sample count and an optional local random seed.

if nargin < 2 || isempty(numSamples), numSamples = 200; end
if nargin < 3, seed = []; end
validateattributes(numSamples, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, mfilename, 'numSamples');
if isempty(seed)
    stream = RandStream.getGlobalStream;
else
    validateattributes(seed, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', '>=', 0, '<=', 2^32-1}, mfilename, 'seed');
    stream = RandStream('mt19937ar', 'Seed', double(seed));
end
rawETC = compute_etc(sequence);
[symbols, ~, ranks] = unique(sequence(:), 'sorted');
L = numel(sequence);
k = numel(symbols);
if k < 2 || k > L - 2
    error('cETC:UndefinedNormalization', ...
        'Normalization requires at least two symbols and k <= L-2.');
end
counts = accumarray(ranks, 1, [k 1]);
commonDivisor = counts(1);
for index = 2:k
    commonDivisor = gcd(commonDivisor, counts(index));
end
lowerBound = k - 1;
capacity = k * commonDivisor;
while capacity < L
    capacity = 2 * capacity;
    lowerBound = lowerBound + 1;
end

% Shuffle positions, preserving each count's association with its label.
rankedSequence = double(ranks(:).');
referenceETC = zeros(double(numSamples), 1);
for sample = 1:numSamples
    reference = rankedSequence(randperm(stream, L));
    referenceETC(sample) = compute_etc(reference);
end
referenceMean = mean(referenceETC);
if referenceMean <= lowerBound
    error('cETC:ZeroEstimatedDenominator', ...
        'The sampled mean equals the lower bound. Increase numSamples.');
end
value = (rawETC - lowerBound) / (referenceMean - lowerBound);
info = struct('raw_etc', rawETC, 'lower_bound', lowerBound, ...
    'reference_mean', referenceMean, 'reference_sd', std(referenceETC), ...
    'reference_mcse', std(referenceETC) / sqrt(double(numSamples)), ...
    'reference_etc', referenceETC, 'num_samples', double(numSamples), ...
    'length', L, 'alphabet_size', k, 'symbols', symbols, 'counts', counts, ...
    'sampling_method', 'uniform permutation', 'seed', seed);
end
