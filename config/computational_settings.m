function settings = computational_settings()
% COMPUTATIONAL_SETTINGS  Returns solver and numerical configuration
%   Provides rigorous, reproducible settings for QFI and eigenvalue
%   computations as used in the manuscript and supplementary material.
%
%   Fields:
%     eig_tol       - convergence tolerance for eigensolver
%     max_iters     - maximum iterations for iterative methods
%     sparse_method - whether to use sparse Lanczos (true/false)
%     random_seed   - seed for random number generators
%     precision     - 'double' or 'single'
%     parallel      - whether to enable MATLAB parallel pool
%     pool_size     - number of workers if parallel is true
%
% Example:
%   settings = computational_settings();

% Eigenproblem settings
settings.eig_tol       = 1e-12;       % target relative tolerance
settings.max_iters     = 1000;        % max Lanczos or Arnoldi iterations
settings.sparse_method = true;        % use sparse solver for N >= 200
settings.precision     = 'double';    % double precision arithmetic

% Randomness
settings.random_seed   = 20250820;    % reproducible seed based on date
rng(settings.random_seed);

% Parallel configuration
settings.parallel      = true;        % enable parallel computation
settings.pool_size     = min(8, feature('numcores')); % up to 8 workers

% If parallel requested, try to start pool
if settings.parallel
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= settings.pool_size
        try
            % force use of local cluster
            parpool('local', settings.pool_size);
        catch ME
            warning('Parallel pool unavailable (%s). Falling back to serial execution.', ME.message);
            settings.parallel = false;
        end
    end
end

end
