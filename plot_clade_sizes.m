% Takes in a full set of species, where each column holds the following:
% 1: id, 2: mass, 3: m_min, 4: death, 5: parent
function clade_sizes = plot_clade_sizes( m )
    % Assuming each ratchet mass is unique to that lineage
    x_mins = unique(m(:, 3));
    
    clade_sizes = zeros(length(x_mins), 1);
    for ii = 1:length(x_mins)
        clade_sizes(ii) = length(find(m(:, 3) == x_mins(ii)));
    end

    figure;
    scatter(x_mins, clade_sizes, 'filled', 'd');
    set(gca, 'xscale', 'log');
    set(gca, 'yscale', 'log');
    xlim([0 10^8]);
    xlabel('x_{min} for clade, g');
    ylabel('Number of species in clade');
    title(sprintf('Sizes of %u clades', length(x_mins)));
end