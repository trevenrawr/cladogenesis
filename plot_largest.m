%default values
min = 1.8;
x_0 = 40;
n = 5000;

%1: id, 2: mass, 3: m_min, 4: death, 5: parent
m = csvread(sprintf('all_%g_%d_%d.csv', min, x_0, n));
[l, ~] = size(m);

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

figure;
plot_dist(m(:, 2));