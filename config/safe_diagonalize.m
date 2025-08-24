function results = safe_diagonalize(H, params, model)
    fprintf('========================================\n');
    fprintf('NON-HERMITIAN BDG SPARSE DIAGONALIZATION\n');
    fprintf('========================================\n');
    fprintf('Model: %s, N = %d\n', model, params.N);
    fprintf('Parameters: t = %.3e, γ = %.3e, Δ = %.3e\n', ...
        params.t, params.gamma, params.Delta);

    % ---- Dimension handling ----
    dim = size(H,1);
    N = params.N;
    if dim ~= 2*N
        warning('Mismatch: Hamiltonian size = %d, expected = 2N = %d', dim, 2*N);
    end

    % How many eigenvalues to compute?
    k = min(dim, 2*N);   % never exceed matrix dimension

    % Method selection
    sparsity = nnz(H) / numel(H);
    fprintf('Using pre-constructed Hamiltonian.\n');
    fprintf('Matrix sparsity: %.2f%% (%d/%d non-zeros)\n', ...
        sparsity*100, nnz(H), numel(H));

    if sparsity < 0.2 && dim > 200
        method = 'eigs_shift_invert';
    else
        method = 'eig_full';
    end
    fprintf('Selected method: %s\n', method);

    % ---- Eigen-decomposition ----
    try
        switch method
            case 'eigs_shift_invert'
                sigma = params.t;  % shift near hopping
                fprintf('shift_invert: computing %d eigenvalues (limited from ''all''), σ = %.3e\n', k, sigma);
                [V,D] = eigs(H, k, sigma);
                d = diag(D);
            otherwise
                [V,D] = eig(full(H));
                d = diag(D);
        end
    catch ME
        warning('Sparse diagonalization failed (%s). Falling back to dense eig.', char(ME.message));
        [V,D] = eig(full(H));
        d = diag(D);
    end

    % ---- Biorthogonal normalization ----
    try
        Vl = inv(V); % left eigenvectors
        err0 = norm(Vl*V - eye(size(H)));
        fprintf('Initial biorthogonal error: %.2e\n', err0);

        % Biorthogonal Gram-Schmidt
        [V,Vl] = biorthogonalize(V,Vl);
        err1 = norm(Vl*V - eye(size(H)));
        fprintf('Final biorthogonal error: %.2e (improvement: %.1fx)\n', ...
            err1, err0/max(err1,eps));
    catch ME
        warning('Biorthogonal normalization failed: %s', char(ME.message));
        Vl = inv(V); % fallback
    end

    % ---- PT-symmetric diagnostics ----
    if strcmp(model,'pt_symmetric')
        try
            P = kron(flip(eye(N)), eye(2));   % parity in BdG space
            comm_err = norm(P*conj(H) - H*P);
            fprintf('PT commutator norm ||[PT,H]|| = %.3e\n', comm_err);
        catch ME
            warning('PT-symmetry check failed: %s', char(ME.message));
        end
    end

    % ---- Save results ----
    results.V = V;      % right eigenvectors
    results.L = Vl;     % left eigenvectors
    results.d = d;      % eigenvalues
    results.method = method;
end

% ---------- Helper ----------
function [Vr, Vl] = biorthogonalize(Vr, Vl)
    % Biorthogonal Gram-Schmidt
    n = size(Vr,2);
    for j = 1:n
        for i = 1:j-1
            coeff = Vl(:,i)'*Vr(:,j);
            Vr(:,j) = Vr(:,j) - Vr(:,i)*coeff;
            Vl(:,j) = Vl(:,j) - Vl(:,i)*conj(coeff);
        end
        norm_factor = sqrt(abs(Vl(:,j)'*Vr(:,j)));
        if norm_factor > 1e-12
            Vr(:,j) = Vr(:,j)/norm_factor;
            Vl(:,j) = Vl(:,j)/conj(norm_factor);
        end
    end
end
