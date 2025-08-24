function numerical_benchmarks(N_list, t)
% NUMERICAL_BENCHMARKS  Benchmark dense vs sparse eigenvalue solvers
%   N_list - vector of system sizes (e.g., [100, 1000, 5000, 10000])
%   t      - hopping amplitude (used to construct Hamiltonian)
%
% Outputs timing table to console.

if nargin < 2
    t = 1.0;
end

% Preallocate results
dense_time  = zeros(size(N_list));
sparse_time = zeros(size(N_list));

for idx = 1:length(N_list)
    N = N_list(idx);
    % Build tridiagonal Hamiltonian
    e = ones(N,1);
    H = spdiags([t*e, zeros(N,1), t*e], [-1,0,1], N, N);
    % Dense solver timing
    fullH = full(H);
    tic;
    eig(fullH);
    dense_time(idx) = toc;
    % Sparse solver timing (compute all eigenvalues via eigs)
    opts.tol   = 1e-12;
    opts.maxit = 10*N;
    tic;
    eigs(H, N, 'smallestabs', opts);
    sparse_time(idx) = toc;
end

% Display table
fprintf('%8s %12s %12s %10s\n', 'N', 'Dense (s)', 'Sparse (s)', 'Speedup');
for i = 1:length(N_list)
    sp = dense_time(i)/sparse_time(i);
    fprintf('%8d %12.4f %12.4f %10.1fx\n', N_list(i), dense_time(i), sparse_time(i), sp);
end

end
