function plot_dist( x )
    edges = logspace(0, 10, 70);
    [y, edges] = histcounts(x, edges, 'Normalization', 'probability');
    scatter(edges(1:end-1), y, 'filled', 'd');
    set(gca, 'xscale', 'log');
    xlim([0 10^8]);
        xlabel('Species mass, g');
        ylabel('Proportion');
end