function experimental_setups()
    % =========================================
    % experimental_setups.m
    % Runs NHSE and PT-symmetric cases with safe diagonalization
    % =========================================

    %% NHSE CASE
    disp('========================================');
    disp('NHSE CASE');
    disp('========================================');

    params = struct();
    params.N = 60;
    params.t = 1.0e7 * 2*pi;       % hopping
    params.gamma = 1.0e8 * pi;     % asymmetry
    params.Delta = 1.0e6;          % pairing
    params.model_type = 'nhse';

    fprintf('Model: %s, N = %d\n', params.model_type, params.N);
    fprintf('Parameters: t = %.3e, γ = %.3e, Δ = %.3e\n', params.t, params.gamma, params.Delta);

    H_NHSE = construct_hamiltonian(params);
    results_nhse = safe_diagonalize(H_NHSE, params, 'nhse');
        % === Biorthogonal overlap diagnostics ===
    if isfield(results_nhse, 'L') && isfield(results_nhse, 'V')
        O = results_nhse.L' * results_nhse.V;
        overlap_error = norm(O - eye(size(O)), 'fro');
        fprintf('Biorthogonal overlap ||O - I|| = %.3e\n', overlap_error);

        % Store in results
        results_nhse.overlap_matrix = O;
        results_nhse.overlap_error = overlap_error;
    end

    %% PT-SYMMETRIC CASE
    disp('========================================');
    disp('PT-SYMMETRIC CASE');
    disp('========================================');

    params = struct();
    params.N = 60;
    params.t = 1.0e7 * 2*pi;       % hopping
    params.gamma = 1.0e7 * 2*pi;   % gain/loss
    params.Delta = 1.0e6;          % pairing
    params.model_type = 'pt_symmetric';

    fprintf('Model: %s, N = %d\n', params.model_type, params.N);
    fprintf('Parameters: t = %.3e, γ = %.3e, Δ = %.3e\n', params.t, params.gamma, params.Delta);

    H_PT = construct_hamiltonian(params);
    results_pt   = safe_diagonalize(H_PT, params, 'pt_symmetric');
        % === Biorthogonal overlap diagnostics ===
    if isfield(results_nhse, 'L') && isfield(results_nhse, 'V')
        O = results_nhse.L' * results_nhse.V;
        overlap_error = norm(O - eye(size(O)), 'fro');
        fprintf('Biorthogonal overlap ||O - I|| = %.3e\n', overlap_error);

        % Store in results
        results_nhse.overlap_matrix = O;
        results_nhse.overlap_error = overlap_error;
    end

    save('experimental_results.mat','results_nhse','results_pt');
    disp('Results saved to experimental_results.mat');
end

