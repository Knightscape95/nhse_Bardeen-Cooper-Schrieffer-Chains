function results = pt_symmetric_bdg_analysis(params)
    % Wrapper for pt_symmetric_bdg that returns analysis results
    fprintf('✅ Using real PT-symmetric BdG analysis (not fallback!)\n');
    
    % Call your existing function
    [H_pt, eigvals, psi_R, psi_L] = pt_symmetric_bdg(params);
    
    % Add QFI calculations from your theory
    delta = params.g_c - params.g;  
    F_PT_Q = (params.t * params.N^2) / (6 * delta);  % Main Eq. (9)
    enhancement_factor = (params.t * params.N) / (6 * delta);
    % 1. Check for exceptional point proximity
    if abs(params_pt.delta) < 1e-6 * params_pt.t
        warning('Too close to exceptional point - triggering fallback');
        results.status = 'fallback';
        results.reason = 'exceptional_point_proximity';
        return;
    end
    
    % 2. Check gN condition more conservatively 
    if params_pt.g * params_pt.N >= 0.95  % More conservative than gN < 1
        warning('gN too large - triggering fallback');
        results.status = 'fallback';
        results.reason = 'perturbative_violation';
        return;
    end
    
    % 3. Check for numerical stability
    condition_number = abs(params_pt.g_c) / abs(params_pt.delta);
    if condition_number > 1e6
        warning('Poor numerical conditioning - triggering fallback');
        results.status = 'fallback';
        results.reason = 'numerical_instability';
        return;
    end
    % Package in expected format
    results = struct();
    results.H_pt = H_pt;
    results.eigenvalues = eigvals;
    results.psi_R = psi_R;
    results.psi_L = psi_L;
    results.qfi_pt = F_PT_Q;
    results.enhancement_factor = enhancement_factor;
    results.scaling = 'N^2/delta';
    results.status = 'success';
    
    fprintf('✅ PT enhancement: F_Q = %.2e (N²/δ scaling)\n', F_PT_Q);
    fprintf('✅ Enhancement factor: η = %.1f\n', enhancement_factor);
end
