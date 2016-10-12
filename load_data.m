%default values
min = 1.8;
x_0 = 40;
n = 5000;

m = csvread(sprintf('all_%g_%d_%d.csv', min, x_0, n));

figure;
plot_dist(m(:, 2));

plot_clade_largest(m);
plot_clade_dists(m);
plot_clade_sizes(m);