function assertDisjointGroups(trainGroup, testGroup)
%ASSERTDISJOINTGROUPS Fail if any group ID is in both training and test.

shared = intersect(trainGroup(:), testGroup(:));
if ~isempty(shared)
    error('exhaustive:SpecimenInBothFolds', ...
        ['%d specimen group(s) appear in both training and test. ' ...
         'First shared ID: %g.'], numel(shared), shared(1));
end
end
