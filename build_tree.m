% [id, mass, m_min, death, parent]
A = csvread('all_1.8_40_5000.csv');

[n, ~] = size(A);

tree = sparse(A(2:end, 1), A(2:end, 5), ones(n-1, 1), n, n);