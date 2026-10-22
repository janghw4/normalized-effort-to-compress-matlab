function result = run_tests
%RUN_TESTS Compare with an explicit pair-by-pair oracle and saved C results.
root = fileparts(fileparts(mfilename('fullpath')));
previousPath = path;
restorePath = onCleanup(@() path(previousPath)); %#ok<NASGU>
addpath(root, '-begin');
originalRng = rng;
restoreRng = onCleanup(@() rng(originalRng)); %#ok<NASGU>
rng(20261022, 'twister');

comparisons = 0;
for alphabet = [2 3]
    if alphabet == 2, maxLength = 12; else, maxLength = 8; end
    for length = 1:maxLength
        for code = 0:alphabet^length-1
            sequence = zeros(1, length);
            residual = code;
            for position = length:-1:1
                sequence(position) = mod(residual, alphabet);
                residual = floor(residual / alphabet);
            end
            assert(compute_etc(sequence) == oracle(sequence), ...
                'Oracle mismatch for %s', mat2str(sequence));
            comparisons = comparisons + 1;
        end
    end
end
for alphabet = [2 3 5 10 20]
    for length = [20 100 500]
        for trial = 1:5
            sequence = randi(alphabet, 1, length) - 1;
            expected = oracle(sequence);
            assert(compute_etc(sequence) == expected);
            assert(compute_etc(sequence(:)) == expected);
            assert(compute_etc(sequence * 100000 - 10000000000) == expected);
            comparisons = comparisons + 3;
        end
    end
end

assert(compute_etc([]) == 0);
assert(compute_etc(ones(1, 100)) == 0);
assert(compute_etc([0 1 0 1]) == 1);
assert(compute_etc(logical([0 1 0 1])) == 1);
assert(compute_etc(uint64([0 1 0 1]) + bitshift(uint64(1), 63)) == 1);
assert(compute_etc([intmin('int64'), intmax('int64')]) == 1);
invalid = {[0 0.5], [0 NaN], [0 Inf], [1 2; 3 4], '0101', {0, 1}, [1 1i]};
for index = 1:numel(invalid)
    must_error(@() compute_etc(invalid{index}), 'compute_etc:InvalidInput');
end
assert(normalize_etc(5, 1, 3) == 2);
assert(normalize_etc(1, 1, 3) == 0);
assert(normalize_etc(3, 1, 3) == 1);
must_error(@() normalize_etc(1, 1, 1), 'normalize_etc:InvalidAnchors');
must_error(@() normalize_etc(0, 1, 3), 'normalize_etc:InvalidAnchors');

fixture = readtable(fullfile(root, 'tests', 'c_reference_cases.csv'), 'TextType', 'string');
for index = 1:height(fixture)
    sequence = sscanf(char(fixture.sequence(index)), '%d').';
    assert(compute_etc(sequence) == fixture.etc(index), 'C mismatch at row %d', index);
end
result = struct('status', 'pass', 'oracle_comparisons', comparisons, ...
    'c_comparisons', height(fixture), 'matlab_version', version);
fprintf('%s\n', jsonencode(result));
end

function steps = oracle(sequence)
% Count each candidate with its own scan; perform every substitution.
steps = 0;
while numel(unique(sequence)) > 1
    candidates = unique([sequence(1:end-1).', sequence(2:end).'], 'rows');
    candidates = sortrows(candidates, [2 1]);
    counts = zeros(size(candidates, 1), 1);
    for index = 1:size(candidates, 1)
        position = 1;
        while position < numel(sequence)
            if isequal(sequence(position:position+1), candidates(index, :))
                counts(index) = counts(index) + 1;
                position = position + 2;
            else
                position = position + 1;
            end
        end
    end
    [~, best] = max(counts);
    pair = candidates(best, :);
    fresh = max(sequence) + 1;
    output = [];
    position = 1;
    while position <= numel(sequence)
        if position < numel(sequence) && isequal(sequence(position:position+1), pair)
            output(end+1) = fresh; %#ok<AGROW>
            position = position + 2;
        else
            output(end+1) = sequence(position); %#ok<AGROW>
            position = position + 1;
        end
    end
    sequence = output;
    steps = steps + 1;
end
end

function must_error(action, identifier)
try
    action();
catch exception
    assert(strcmp(exception.identifier, identifier), 'Unexpected error: %s', exception.identifier);
    return
end
error('run_tests:MissingError', 'Expected %s', identifier);
end
