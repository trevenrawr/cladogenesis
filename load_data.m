%default values
min = 1.8;
x_0 = 40;
n = 5000;

to_load = [4 40 400 4000 40000 400000];

curr = 1;
for x_0 = to_load
    fid = fopen(sprintf('extant_%g_%d_%d.csv', min, x_0, n));
    % Just load the mass data, thus the * in %*u
    m(curr) = textscan(fid, '%f %*u %*u', 'Delimiter', ',');
    fclose(fid);
    curr = curr + 1;
end

figure;
for ii = 1:length(m)
    plot_dist(m{ii});
    hold on;
end
legend(num2str(to_load'));