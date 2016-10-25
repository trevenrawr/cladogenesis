% Takes in a full set of species, where each column holds the following:
% 1: id, 2: mass, 3: m_min, 4: death, 5: parent
function plot_clade_spawnrate( m )
    [n, ~] = size(m);

    % Assuming each ratchet mass is unique to that lineage
    x_mins = unique(m(:, 3));
    
    n_bins = 1000;
    bin_size = n / n_bins;
    edges = linspace(0, n, n_bins);

    figure;
    for ii = 1:length(x_mins)
        clade = find(m(:, 3) == x_mins(ii));
        id_clade = m(clade, 1);  % IDs, but also birthdates!
        
        [y, ~] = histcounts(id_clade, edges);

        plot(edges(1:end-1), y);
        hold on;
    end
    xlabel('Model Time');
    ylabel(sprintf('Speciation rate (num born per %0.4g model steps)', bin_size));
end