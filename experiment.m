format longEng
load minnesota;
A = Problem.A;

rho = 3.232397; % matrix norm of `minnesota`, cfr. https://sparse.tamu.edu/Gleich/minnesota
alpha = 0.85 / rho;
L_max = 30;
tol = 1e-4;
pcg_tol = 1e-5;

runs = 30;

n = size(A, 1);

w = 1011;
u = 1011;
v = 1015;

I = speye(n);
one = ones(n, 1);
zero = zeros(n, 1);

x = (I - alpha * A) \ one; % Katz vector

% Without one edge

A_eps = A;
A_eps(u, v) = 0;
A_eps(v, u) = 0;

time_i_edge = 0;
time_ii_edge = 0;
time_iii_edge = 0;

L_edge = 0;
it_ii_edge = 0;
it_iii_edge = 0;

for i = 1:runs
    a = tic;
    [x_paper_edge, L] = katz_edge(A, x, alpha, L_max, tol, u, v);
    time_i_edge = time_i_edge + toc(a);
    L_edge = L_edge + L;

    a = tic;
    [x_pcg_zero_edge, ~, ~, it] = pcg(I - alpha * A_eps, one, pcg_tol, [], [], [], zero);
    time_ii_edge = time_ii_edge + toc(a);
    it_ii_edge = it_ii_edge + it;

    a = tic;
    [x_pcg_katz_edge, ~, ~, it] = pcg(I - alpha * A_eps, one, pcg_tol, [], [], [], x);
    time_iii_edge = time_iii_edge + toc(a);
    it_iii_edge = it_iii_edge + it;
end

time_i_edge = time_i_edge / runs;
time_ii_edge = time_ii_edge / runs;
time_iii_edge = time_iii_edge / runs;

L_edge = L_edge / runs;
it_ii_edge = it_ii_edge / runs;
it_iii_edge = it_iii_edge / runs;

fprintf('\n===== Without edge {%d, %d} =====\n', u, v);
fprintf('%-16s %-20s %-25s\n', 'Method', 'Iter/L', 'Time (s)');
fprintf('%-15s L = %-5.2f   %25.16e\n', '(i)', L_edge, time_i_edge);
fprintf('%-15s %7.2f     %25.16e\n', '(ii)', it_ii_edge, time_ii_edge);
fprintf('%-15s %7.2f     %25.16e\n', '(iii)', it_iii_edge, time_iii_edge);

rel_diff = norm(x_paper_edge - x_pcg_katz_edge) / norm(x_pcg_katz_edge);
fprintf('\nRelative difference between `katz_edge` and `pcg` (with guess = x): ');
fprintf('%.16e\n', rel_diff);

% Without one node

A_n = A;
A_n(:, w) = zero;
A_n(w, :) = zero';

time_i_node = 0;
time_ii_node = 0;
time_iii_node = 0;

L_node = 0;
it_ii_node = 0;
it_iii_node = 0;

for i = 1:runs
    a = tic;
    [x_paper_node, L] = katz_node(A, x, alpha, L_max, tol, w);
    time_i_node = time_i_node + toc(a);
    L_node = L_node + L;

    a = tic;
    [x_pcg_zero_node, ~, ~, it] = pcg(I - alpha * A_n, one, pcg_tol, [], [], [], zero);
    time_ii_node = time_ii_node + toc(a);
    it_ii_node = it_ii_node + it;

    a = tic;
    [x_pcg_katz_node, ~, ~, it] = pcg(I - alpha * A_n, one, pcg_tol, [], [], [], x);
    time_iii_node = time_iii_node + toc(a);
    it_iii_node = it_iii_node + it;
end

time_i_node = time_i_node / runs;
time_ii_node = time_ii_node / runs;
time_iii_node = time_iii_node / runs;

L_node = L_node / runs;
it_ii_node = it_ii_node / runs;
it_iii_node = it_iii_node / runs;

fprintf('\n===== Without node %d =====\n', w);
fprintf('%-16s %-20s %-25s\n', 'Method', 'Iter/L', 'Time (s)');
fprintf('%-15s L = %-5.2f   %25.16e\n', '(i)', L_node, time_i_node);
fprintf('%-15s %7.2f     %25.16e\n', '(ii)', it_ii_node, time_ii_node);
fprintf('%-15s %7.2f     %25.16e\n', '(iii)', it_iii_node, time_iii_node);

rel_diff = norm(x_paper_node - x_pcg_katz_node) / norm(x_pcg_katz_node);
fprintf('\nRelative difference between `katz_node` and `pcg` (with guess = x): ');
fprintf('%.16e\n', rel_diff);
