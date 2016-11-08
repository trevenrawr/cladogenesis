% Takes in a full set of species, where each column holds the following:
% 1: id, 2: mass, 3: m_min, 4: death, 5: parent
function plot_clade_n_extant( m )
    [n, ~] = size(m);

    % Assuming each ratchet mass is unique to that lineage
    x_mins = unique(m(:, 3));
    
    n_bins = 1000;
    bin_size = floor((n-1)/2 / n_bins);
    edges = 1:bin_size:(n-1)/2;

    figure;
    % For every clade
    for ii = 1:length(x_mins)
        clade = find(m(:, 3) == x_mins(ii));
        ls_clade = m(clade, [1, 4]);  % Birthdates and death dates
        
        % Keep track of how many are alive at each "edge"
        extants = zeros(length(edges), 1, 'uint16');
        for jj = 1:length(edges)
            % Check every species in the clade
            for kk = 1:length(clade)
                birth = floor(ls_clade(kk, 1) / 2) + 1;

                death = ls_clade(kk, 2);
                
                % If the birth is before the edge and the death is after
                if birth <= edges(jj) && death >= edges(jj)
                    extants(jj) = extants(jj) + 1;
                end
            end
        end
        plot(extants);
        hold on;
    end
    
    xlabel('Model Time');
    ylabel('Number of extant species in clade');
end