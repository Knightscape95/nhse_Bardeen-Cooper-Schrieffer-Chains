function params = params_nhse()
% NHSE Configuration - Non-Hermitian Skin Effect Model
% Based on Eq. (2): H_NHSE = Σ_j [(t+γ)c†_j c_{j+1} + (t-γ)c†_{j+1} c_j]

    params = struct();
    
    % Model identification
    params.model = 'nhse';
    params.description = 'Non-Hermitian Skin Effect with asymmetric hoppings';
    
    % System parameters (SI units: rad/s)
    params.N = 60;                    % System size (sites)
    params.t = 2*pi*10e6;            % Hopping amplitude: t/2π = 10 MHz
    params.Delta = 2*pi*1e6;         % BdG pairing: Δ/2π = 1 MHz
    
    % NHSE-specific parameters - CORRECTED
    params.gamma = 1.2 * params.t;   % Non-reciprocity: γ > t for strong NHSE
    
    % Additional parameters for multiparameter estimation
    params.mu = 0.1 * params.t;      % Chemical potential: |μ| < 2t (metallic regime)
    params.phi = 0.05;               % Peierls phase: |φ| ≪ 1 (perturbative regime)
    
    % CORRECTED: Proper localization parameter calculation
    gamma_over_t = abs(params.gamma) / params.t;
    if gamma_over_t > 1
        % Strong NHSE regime: κ = arccosh(|γ|/t)
        params.kappa = acosh(gamma_over_t);
        params.localization_length = 1/params.kappa;
        params.nhse_regime = 'strong';
    else
        % Weak non-reciprocity: use alternative formula
        params.localization_length = 1/log((params.t + abs(params.gamma))/(params.t - abs(params.gamma)));
        params.kappa = 1/params.localization_length;
        params.nhse_regime = 'weak';
    end
    
    params.expected_qfi_suppression = exp(-2*params.kappa*params.N);
    
    % Experimental constraints
    params.coherence_time = 100e-6;  % T₂* = 100 μs (superconducting circuits)
    params.decoherence_rate = 2*pi*50e3; % γ_φ/2π ~ 50 kHz
    
    % CORRECTED: Validity checks with warnings instead of errors
    if abs(params.mu) >= 2*params.t
        warning('μ = %.2f*t violates |μ| < 2t condition for metallic regime', abs(params.mu)/params.t);
    end
    
    if abs(params.phi) >= 1
        warning('φ = %.2f violates |φ| ≪ 1 condition for perturbative validity', abs(params.phi));
    end
    
    if params.N * params.kappa <= 1
        warning('κN = %.2f should be > 1 for clear NHSE regime (current regime: %s)', ...
                params.N * params.kappa, params.nhse_regime);
    end
    
    fprintf('NHSE Configuration Loaded:\n');
    fprintf('γ/t = %.2f (%s NHSE regime)\n', gamma_over_t, params.nhse_regime);
    fprintf('κ = %.3f, κN = %.2f\n', params.kappa, params.N * params.kappa);
    fprintf('Expected QFI suppression: %.2e\n', params.expected_qfi_suppression);
    
end
