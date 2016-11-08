% Takes in a full set of species, where each column holds the following:
% 1: id, 2: mass, 3: m_min, 4: death, 5: parent
function plot_n_extant( m )
    [n, ~] = size(m);

    extants = zeros((n-1)/2, 1, 'uint16');
    for ii = 1:n
        birth = floor(m(ii, 1) / 2) + 1;
        
        death = m(ii, 4);
        if death > length(extants)
            death = length(extants);
        end
        
        extants(birth:death) = extants(birth:death) + 1;
    end
    
    figure;
    plot(extants);
    xlabel('Model Time');
    ylabel('Number of extant species in clade');
end