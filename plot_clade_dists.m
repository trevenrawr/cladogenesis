% Takes in a full set of species, where each column holds the following:
% 1: id, 2: mass, 3: m_min, 4: death, 5: parent
function plot_clade_dists( m )
    % Assuming each ratchet mass is unique to that lineage
    x_mins = unique(m(:, 3));
    
    % Pick the smallest clade size to plot
    cutoff = 500;

    figure;
    for ii = 1:length(x_mins)
        clade = m(:, 3) == x_mins(ii);
        if sum(clade) > cutoff
            m_clade = m(clade, 2);

            edges = logspace(0, 10, 70);
            [y, edges] = histcounts(m_clade, edges);
            scatter(edges(1:end-1), y, 'filled', 'd');

            hold on;
        end
    end
    
    set(gca, 'xscale', 'log');
    set(gca, 'yscale', 'log');
    xlim([0 10^8]);
    xlabel('Species mass, g');
    ylabel('Species count');
    title(sprintf('Clades with over %u species', cutoff));
end