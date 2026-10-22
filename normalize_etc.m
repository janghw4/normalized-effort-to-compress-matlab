function value = normalize_etc(etc, lowerBound, referenceMean)
%NORMALIZE_ETC Apply supplied global or composition-conditioned anchors.
%   VALUE = (ETC - LOWERBOUND) / (REFERENCEMEAN - LOWERBOUND).
%   Use global anchors for nETC and labeled-composition anchors for cETC.
%   Estimate the reference mean separately. Values above one are retained.

validateattributes(etc, {'numeric'}, {'real', 'finite', 'scalar'}, mfilename, 'etc');
validateattributes(lowerBound, {'numeric'}, {'real', 'finite', 'scalar'}, mfilename, 'lowerBound');
validateattributes(referenceMean, {'numeric'}, {'real', 'finite', 'scalar'}, mfilename, 'referenceMean');
etc = double(etc);
lowerBound = double(lowerBound);
referenceMean = double(referenceMean);
if etc < lowerBound || referenceMean <= lowerBound
    error('normalize_etc:InvalidAnchors', ...
        'Require ETC >= lower bound and reference mean > lower bound.');
end
value = (etc - lowerBound) / (referenceMean - lowerBound);
end
