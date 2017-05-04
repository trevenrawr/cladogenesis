distances = [];
for ii = 1:length(no_heuristic) - 1
    if no_heuristic(ii) > no_heuristic(ii + 1)
        distances(length(distances) + 1, :) = [no_heuristic(ii), no_heuristic(ii) - no_heuristic(ii + 1)];
    end
end