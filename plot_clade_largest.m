% Takes in a full set of species, where each column holds the following:
% 1: id, 2: mass, 3: m_min, 4: death, 5: parent
function plot_clade_largest( m )
    % Assuming each ratchet mass is unique to that lineage
    x_mins = unique(m(:, 3));

    figure;
    for ii = 1:length(x_mins)
        clade = find(m(:, 3) == x_mins(ii));
        m_clade = m(clade, 2);

        % Write over smaller sizes, so we plot only the largest species see yet
        for jj = 2:length(clade)
            if m_clade(jj) < m_clade(jj - 1)
                m_clade(jj) = m_clade(jj - 1);
            end
        end

        plot(clade, m_clade);
        hold on;
    end
    set(gca, 'yscale', 'log');
    xlabel('Model time');
    ylabel('Largest mass in clade');
end