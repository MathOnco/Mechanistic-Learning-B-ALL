function run_flow_abundance_nested_10fold(outcomes, models)
%RUN_FLOW_ABUNDANCE_NESTED_10FOLD Exhaustive nested search over x1..x4.

if nargin < 1, outcomes = []; end
if nargin < 2, models = []; end
run_exhaustive_nested_10fold(outcomes, models, "flow");
end
