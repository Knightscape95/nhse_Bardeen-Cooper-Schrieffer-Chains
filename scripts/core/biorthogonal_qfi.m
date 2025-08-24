function F = biorthogonal_qfi(system, param, opts)
% BIORTHOGONAL_QFI  Quantum Fisher information for non-Hermitian systems
% using biorthogonal eigenvectors and first-order biorthogonal perturbation theory.
%
%   F = biorthogonal_qfi(system, param)
%   F = biorthogonal_qfi(system, param, opts)
%
% REQUIRED FIELDS in 'system':
%   system.psi_R     : right eigenvectors  (N x N)
%   system.psi_L     : left  eigenvectors  (N x N)
%   system.eigenvals : eigenvalues         (N x 1) or (1 x N)
%   system.H_params  : struct with fields needed by construct_hamiltonian
%
% OPTIONAL:
%   system.V         : precomputed perturbation operator dH/d(param) (N x N or 2N x 2N)
%
% PARAM:
%   param            : char/string, parameter name in H_params
%
% OPTS (optional struct fields):
%   opts.state       : integer, which eigenstate column (default: 1)
%   opts.fd_step     : finite-difference step for dH (default: 1e-6)
%   opts.tol         : small tolerance to avoid dividing by tiny gaps (default: 1e-12)
%
% OUTPUT:
%   F                : scalar QFI for the selected eigenstate

    %% ---------------- Input validation ----------------
    if ~(isstruct(system))
        error('MATLAB:biorthogonal_qfi:invalidInput', 'Input "system" must be a struct.');
    end
    reqFields = {'psi_R','psi_L','eigenvals','H_params'};
    for f = reqFields
        if ~isfield(system, f{1})
            error('MATLAB:biorthogonal_qfi:invalidInput', ...
                'System struct must contain field "%s".', f{1});
        end
    end
    if ~(ischar(param) || isstring(param))
        error('MATLAB:biorthogonal_qfi:invalidInput', 'Parameter name must be a char or string.');
    end
    param = char(param);

    psi_R = system.psi_R;
    psi_L = system.psi_L;
    E     = system.eigenvals(:);  % ensure column

    [N1, M1] = size(psi_R);
    [N2, M2] = size(psi_L);
    if N1 ~= N2 || M1 ~= M2 || N1 ~= numel(E)
        error('MATLAB:biorthogonal_qfi:invalidInput', ...
            'Dimensions must match: psi_R, psi_L are N×N and length(eigenvals)=N.');
    end
    N = N1;

    %% ---------------- Options ----------------
    if nargin < 3 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'state'),   opts.state = 1;      end
    if ~isfield(opts,'fd_step'), opts.fd_step = 1e-6; end
    if ~isfield(opts,'tol'),     opts.tol = 1e-12;    end
    n = opts.state;
    if ~(isscalar(n) && n>=1 && n<=N && n==floor(n))
        error('MATLAB:biorthogonal_qfi:invalidInput', 'opts.state must be an integer in [1..N].');
    end

    %% ---------------- Robust Biorthonormalization ----------------
    [U,S,V] = svd(psi_L'*psi_R);
    psi_R = psi_R * V * diag(1./sqrt(diag(S)));
    psi_L = psi_L * U * diag(1./sqrt(diag(S)));

    %% ---------------- Build V = dH/d(param) ----------------
    if isfield(system, 'V') && ~isempty(system.V)
        V = system.V;
        if ~isequal(size(V,1),N) || ~isequal(size(V,2),N)
            error('MATLAB:biorthogonal_qfi:invalidInput', 'Provided V must be N×N to match eigenvector dimensions.');
        end
    else
        if ~isfield(system.H_params, param)
            error('MATLAB:biorthogonal_qfi:invalidInput', ...
                'Parameter "%s" not found in system.H_params.', param);
        end
        hp = system.H_params;
        h  = opts.fd_step;

        hp_plus  = hp;  hp_plus.(param)  = hp.(param) + h;
        hp_minus = hp;  hp_minus.(param) = hp.(param) - h;

        Hplus  = construct_hamiltonian(hp_plus);
        Hminus = construct_hamiltonian(hp_minus);

        if ~isequal(size(Hplus), [N N]) || ~isequal(size(Hminus), [N N])
            error('MATLAB:biorthogonal_qfi:invalidInput', ...
                 'Inconsistent sizes from construct_hamiltonian during finite-difference.');
        end
        V = (Hplus - Hminus) / (2*h);
    end

    %% ---------------- d|psi_R^n>, d<psi_L^n| via perturbation theory ----------------
    dpsi_R_n = zeros(N,1);
    dpsi_L_n = zeros(N,1);
    En = E(n);

    % Precompute matrix elements
    M = psi_L' * (V * psi_R);

    for m = 1:N
        if m == n, continue; end
        gap = En - E(m);
        if abs(gap) < opts.tol
            gap = sign(gap) * opts.tol;  % regularize small gap
        end
        coeff_R = - M(m,n) / gap;
        dpsi_R_n = dpsi_R_n + psi_R(:,m) * coeff_R;

        coeff_L = - M(n,m) / gap;
        dpsi_L_n = dpsi_L_n + psi_L(:,m) * conj(coeff_L);
    end

    %% ---------------- QFI for the selected state ----------------
    psiR_n = psi_R(:,n);
    psiL_n = psi_L(:,n);

    term1 = real( (dpsi_L_n') * dpsi_R_n );
    term2 = abs( (psiL_n') * dpsi_R_n )^2;

    F = 4 * (term1 - term2);

    % Ensure real non-negative
    F = real(F);
    F = max(F,0);
end
