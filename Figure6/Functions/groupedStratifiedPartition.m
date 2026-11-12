function cvp = groupedStratifiedPartition(y, group, requestedFolds, seed)
%GROUPEDSTRATIFIEDPARTITION Stratified K-fold split that keeps groups intact.
%
% Every row sharing a group ID is assigned to the same fold. A group
% therefore appears entirely in training or entirely in test. Groups are
% stratified by outcome, so the outcome label must be constant within each
% group. The returned object supports training() and test().

y = double(y(:));
group = group(:);
if numel(group) ~= numel(y)
    error('mechlearn:GroupLengthMismatch', ...
        'Expected one group ID per observation; found %d IDs and %d outcomes.', ...
        numel(group), numel(y));
end
if any(~isfinite(group))
    error('mechlearn:NonfiniteGroup', 'Group IDs must be finite.');
end

[groups, ~, groupIndex] = unique(group, 'stable');
yGroup = zeros(numel(groups), 1);
for g = 1:numel(groups)
    labels = unique(y(groupIndex == g));
    if numel(labels) ~= 1
        error('mechlearn:MixedGroupOutcome', ...
            'Group %g has more than one outcome label.', groups(g));
    end
    yGroup(g) = labels;
end

specimenPartition = stratifiedPartition( ...
    yGroup, requestedFolds, seed);
foldOfSpecimen = zeros(numel(groups), 1);
for f = 1:specimenPartition.NumTestSets
    foldOfSpecimen(test(specimenPartition, f)) = f;
end
if any(foldOfSpecimen == 0)
    error('mechlearn:UnassignedGroup', ...
        'At least one specimen was not assigned to a fold.');
end

foldOfRow = foldOfSpecimen(groupIndex);
cvp = cvpartition(CustomPartition=foldOfRow);
end
