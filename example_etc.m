function example_etc
%EXAMPLE_ETC Run a few symbolic sequences and an illustrative normalization.
sequences = {[0 1 0 1], [7 7 7], [0 1 2]};
for index = 1:numel(sequences)
    fprintf('%s: ETC = %d\n', mat2str(sequences{index}), compute_etc(sequences{index}));
end
% Illustrative anchors; replace these with matched reference estimates.
fprintf('Normalized value = %.1f\n', normalize_etc(5, 1, 3));
end
