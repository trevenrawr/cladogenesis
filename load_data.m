%default values
model = 1;
min = 1.8;
x_0 = 40;
n = 5000;

model = 3;

m = csvread(sprintf('all_m%d_%g_%d_%d.csv', model, min, x_0, n));

plot_clade_largest(m);
plot_clade_spawnrate(m);
plot_clade_dists(m);
plot_clade_sizes(m);
plot_clade_n_extant(m);

m_ext = csvread(sprintf('extant_m%d_%g_%d_%d.csv', model, min, x_0, n));

figure;
plot_dist(m_ext(:, 2));
title('Extant species');