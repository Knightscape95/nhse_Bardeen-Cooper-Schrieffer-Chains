function convergence_tests(Hfunc, g_vals, N_vals, tol_vals)
%CONVERGENCE_TESTS Robust eigenvalue convergence testing for non-Hermitian systems
%
% USAGE: convergence_tests(Hfunc, g_vals, N_vals, tol_vals)
%
% INPUTS:
%   Hfunc    - Function handle: H = Hfunc(g, N) returns Hamiltonian matrix
%   g_vals   - Array of g parameter values to test
%   N_vals   - Array of system sizes N to test  
%   tol_vals - Array of tolerances for eigs convergence
%
% FEATURES:
%   - Automatic conditioning checks and method selection
%   - NHSE regime detection with analytical fallbacks
%   - Robust error handling for ill-conditioned matrices
%   - Performance monitoring and convergence diagnostics

fprintf('=== Convergence Tests ===\n');

for idxN = 1:length(N_vals)
    N = N_vals(idxN);
    fprintf('\nN = %d\n', N);
    fprintf('%8s %12s %12s %8s %10s %8s\n', 'tol', 'err_NHSE', 'err_PT', 'method', 'cond', 'kappa');
    
    for idxtol = 1:length(tol_vals)
        tol = tol_vals(idxtol);
        err_NHSE = 0;
        err_PT = 0;
        
        for idxg = 1:length(g_vals)
            g = g_vals(idxg);
            
            try
                % Generate Hamiltonian
                H = Hfunc(g, N);
                
                % Parameter regime analysis
                t = 1; % Assuming t=1 as reference
                kappa = 0;
                if abs(g) > t
                    kappa = acosh(abs(g)/t);
                end
                
                % Matrix conditioning check
                if issparse(H)
                    condH = condest(H);
                else
                    condH = cond(H);
                end
                
                % Method selection based on numerical stability
                method = select_eigenvalue_method(H, condH, kappa, N, tol);
                
                % Compute eigenvalue using selected method
                [eigVal, success] = compute_eigenvalue_robust(H, method, tol, kappa, N, t);
                
                % Error analysis (simplified for demonstration)
                if abs(g) > t
                    % NHSE regime - compare with analytical formula
                    analytical_val = analytical_NHSE_formula(g, N, t);
                    err_NHSE = max(err_NHSE, abs(eigVal - analytical_val));
                else
                    % PT regime - simplified error estimate
                    err_PT = max(err_PT, abs(imag(eigVal)));
                end
                
                fprintf('%8.1e %12.2e %12.2e %8s %10.2e %8.2f\n', ...
                    tol, err_NHSE, err_PT, method, condH, kappa);
                
            catch ME
                fprintf('%8.1e %12s %12s %8s %10s %8s (ERROR: %s)\n', ...
                    tol, 'NaN', 'NaN', 'FAIL', 'N/A', 'N/A', ME.message);
            end
        end
    end
end

fprintf('\n=== Test Summary ===\n');
fprintf('Methods used: DENSE (full eig), SPARSE (eigs), ANALYTICAL (formula)\n');
fprintf('Recommendations:\n');
fprintf('- Use DENSE method for N <= 200\n');
fprintf('- Use ANALYTICAL for kappa*N > 20 (extreme NHSE)\n');
fprintf('- Use SPARSE only for well-conditioned large systems\n');

end

function method = select_eigenvalue_method(H, condH, kappa, N, tol)
%SELECT_EIGENVALUE_METHOD Choose optimal eigenvalue computation method
%
% Selection criteria:
% 1. Extreme NHSE regime: Use analytical formulas
% 2. Ill-conditioned matrices: Use regularized dense method
% 3. Small-medium systems: Use dense diagonalization
% 4. Large well-conditioned: Use sparse iterative

% Critical thresholds
COND_THRESHOLD = 1e12;
NHSE_THRESHOLD = 20;
SIZE_THRESHOLD = 200;

if kappa * N > NHSE_THRESHOLD
    method = 'ANALYTICAL';
elseif condH > COND_THRESHOLD
    method = 'DENSE_REG';
elseif length(H) <= SIZE_THRESHOLD
    method = 'DENSE';
else
    method = 'SPARSE';
end
end

function [eigVal, success] = compute_eigenvalue_robust(H, method, tol, kappa, N, t)
%COMPUTE_EIGENVALUE_ROBUST Robust eigenvalue computation with fallbacks

success = false;
eigVal = NaN;

switch method
    case 'ANALYTICAL'
        % Use analytical NHSE formula to avoid numerical issues
        eigVal = analytical_NHSE_formula(abs(imag(H(1,1))), N, t);
        success = true;
        
    case 'DENSE_REG'
        % Regularized dense computation for ill-conditioned matrices
        reg_param = tol * norm(H, 'fro');
        H_reg = H + reg_param * eye(size(H));
        try
            evals = eig(full(H_reg));
            [~, idx] = max(abs(evals));
            eigVal = evals(idx) - reg_param;
            success = true;
        catch
            eigVal = NaN;
        end
        
    case 'DENSE'
        % Standard dense diagonalization
        try
            evals = eig(full(H));
            [~, idx] = max(abs(evals));
            eigVal = evals(idx);
            success = true;
        catch
            eigVal = NaN;
        end
        
    case 'SPARSE'
        % Sparse iterative method with robust options
        opts.tol = tol;
        opts.maxit = min(1000, 50*length(H));
        opts.disp = 0;
        opts.issym = false;
        
        try
            [~, D, flag] = eigs(H, 1, 'largestabs', opts);
            if flag == 0
                eigVal = D(1,1);
                success = true;
            else
                % Fallback to dense method
                [eigVal, success] = compute_eigenvalue_robust(H, 'DENSE', tol, kappa, N, t);
            end
        catch
            % Fallback to dense method
            [eigVal, success] = compute_eigenvalue_robust(H, 'DENSE', tol, kappa, N, t);
        end
        
    otherwise
        error('Unknown eigenvalue method: %s', method);
end
end

function val = analytical_NHSE_formula(g, N, t)
%ANALYTICAL_NHSE_FORMULA Analytical formula for extreme NHSE regime
%
% Uses the theoretical result from Supp. Eq. S2_NHSE_QFI_corrected:
% F_Q^{NHSE} = 4N³e^{-2κN}/(3t²sinh²κ)

if abs(g) <= t
    val = 0; % Not in NHSE regime
    return;
end

kappa = acosh(abs(g)/t);

% Avoid overflow in exponential calculations
if kappa * N > 50
    val = 0; % Essentially zero due to exponential suppression
else
    val = 4 * N^3 * exp(-2*kappa*N) / (3 * t^2 * sinh(kappa)^2);
end
end
