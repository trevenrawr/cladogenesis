figure;
cdfplot(log10(m4));
hold on;
cdfplot(log10(m40));

[~, p1, k1] = kstest2(m4, m40)
[~, p2, k2] = wkstest(m4, m40)