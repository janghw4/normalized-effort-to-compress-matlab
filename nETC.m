function [value, info] = nETC(sequence, numSamples, seed)
%NETC Global-ensemble normalized ETC, including both reference anchors.
%   VALUE = NETC(SEQUENCE) uses 200 uniform full-alphabet references.
%   [VALUE, INFO] = NETC(SEQUENCE, NUMSAMPLES, SEED) sets the sample count
%   and an optional local random seed. INFO contains the estimated anchors.

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
    error('nETC:UndefinedNormalization', ...
        'Normalization requires at least two symbols and k <= L-2.');
end
counts = accumarray(ranks, 1, [k 1]);

% Largest divisor of L no greater than floor(L/k).
terminal = 1;
limit = floor(L / k);
for divisor = 1:floor(sqrt(L))
    if mod(L, divisor) == 0
        paired = L / divisor;
        if divisor <= limit, terminal = max(terminal, divisor); end
        if paired <= limit, terminal = max(terminal, paired); end
    end
end
lowerBound = k - 1;
capacity = k * terminal;
while capacity < L
    capacity = 2 * capacity;
    lowerBound = lowerBound + 1;
end

% Rejection is efficient when the missing-symbol union bound is small.
useRejection = log(k) + L * log1p(-1/k) <= log(0.25);
if useRejection
    method = 'uniform rejection';
else
    method = 'completion count';
    logD = completion_counts(L, k);
end
referenceETC = zeros(double(numSamples), 1);
for sample = 1:numSamples
    if useRejection
        reference = randi(stream, k, 1, L);
        while numel(unique(reference)) ~= k
            reference = randi(stream, k, 1, L);
        end
    else
        reference = draw_full_alphabet(L, k, logD, stream);
    end
    referenceETC(sample) = compute_etc(reference);
end
referenceMean = mean(referenceETC);
if referenceMean <= lowerBound
    error('nETC:ZeroEstimatedDenominator', ...
        'The sampled mean equals the lower bound. Increase numSamples.');
end
value = (rawETC - lowerBound) / (referenceMean - lowerBound);
info = struct('raw_etc', rawETC, 'lower_bound', lowerBound, ...
    'reference_mean', referenceMean, 'reference_sd', std(referenceETC), ...
    'reference_mcse', std(referenceETC) / sqrt(double(numSamples)), ...
    'reference_etc', referenceETC, 'num_samples', double(numSamples), ...
    'length', L, 'alphabet_size', k, 'symbols', symbols, 'counts', counts, ...
    'sampling_method', method, 'seed', seed);
end

function logD = completion_counts(L, k)
% D(r,u): completions with r positions left and u labels still unseen.
% Logarithms keep the integer-count recurrence within floating-point range.
logD = -inf(L + 1, k + 1);
logD(1, 1) = 0;
for remaining = 1:L
    for unseen = 0:min(k, remaining)
        old = -inf;
        fresh = -inf;
        if unseen < k
            old = log(k - unseen) + logD(remaining, unseen + 1);
        end
        if unseen > 0
            fresh = log(unseen) + logD(remaining, unseen);
        end
        largest = max(old, fresh);
        logD(remaining + 1, unseen + 1) = ...
            largest + log(exp(old - largest) + exp(fresh - largest));
    end
end
end

function reference = draw_full_alphabet(L, k, logD, stream)
reference = zeros(1, L);
unseenLabels = 1:k;
seenLabels = zeros(1, k);
unseen = k;
for position = 1:L
    remaining = L - position + 1;
    if unseen == 0
        chooseNew = false;
    elseif unseen == k || unseen == remaining
        chooseNew = true;
    else
        probability = exp(log(unseen) + logD(remaining, unseen) ...
            - logD(remaining + 1, unseen + 1));
        chooseNew = rand(stream) < probability;
    end
    if chooseNew
        index = randi(stream, unseen);
        label = unseenLabels(index);
        unseenLabels(index) = unseenLabels(unseen);
        unseen = unseen - 1;
        seenLabels(k - unseen) = label;
    else
        label = seenLabels(randi(stream, k - unseen));
    end
    reference(position) = label;
end
end
