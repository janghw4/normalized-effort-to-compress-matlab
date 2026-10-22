function steps = compute_etc(sequence)
%COMPUTE_ETC Count pair substitutions until a symbolic sequence is constant.
%   STEPS = COMPUTE_ETC(SEQUENCE) accepts a row or column vector of finite
%   integer labels. Numerical label order determines tie breaking.
%   Empty and constant inputs return zero.

if ~(isnumeric(sequence) || islogical(sequence)) || ~isreal(sequence) || ...
        ~ismatrix(sequence) || (~isvector(sequence) && ~isempty(sequence)) || ...
        any(~isfinite(sequence(:))) || any(sequence(:) ~= fix(sequence(:)))
    error('compute_etc:InvalidInput', ...
        'Input must be a real vector of finite integer labels.');
end
steps = 0;
if isempty(sequence)
    return
end

% Rank the original labels without changing their numerical order.
% UNIQUE operates on the input type, preserving distinct int64/uint64 labels.
[labels, ~, ranks] = unique(sequence(:), 'sorted');
current = double(ranks(:).');
nextLabel = numel(labels) + 1;
currentLength = numel(current);
buffer = zeros(size(current));

while currentLength > 1
    if all(current(1:currentLength) == current(1))
        return
    end

    left = zeros(1, currentLength - 1);
    right = zeros(1, currentLength - 1);
    pairCount = 0;
    position = 1;
    while position < currentLength
        a = current(position);
        b = current(position + 1);
        pairCount = pairCount + 1;
        left(pairCount) = a;
        right(pairCount) = b;
        % Equal-symbol runs contribute disjoint pairs from their left edge.
        if a == b && position < currentLength - 1 && current(position + 2) == a
            position = position + 1;
        end
        position = position + 1;
    end
    counts = sparse(left(1:pairCount), right(1:pairCount), ...
        ones(1, pairCount), nextLabel - 1, nextLabel - 1);
    % FIND visits columns first: smallest second label, then first label.
    [firstLabels, secondLabels, frequencies] = find(counts);
    [frequency, selected] = max(frequencies);
    if frequency == 1
        steps = steps + currentLength - 1;
        return
    end
    a = firstLabels(selected);
    b = secondLabels(selected);

    position = 1;
    written = 0;
    while position <= currentLength
        written = written + 1;
        if position < currentLength && current(position) == a && current(position + 1) == b
            buffer(written) = nextLabel;
            position = position + 2;
        else
            buffer(written) = current(position);
            position = position + 1;
        end
    end
    previous = current;
    current = buffer;
    buffer = previous;
    currentLength = written;
    nextLabel = nextLabel + 1;
    steps = steps + 1;
end
end
