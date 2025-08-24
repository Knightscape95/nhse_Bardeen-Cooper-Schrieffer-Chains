function params = params_pt()
% PT-Symmetric Configuration - Balanced Gain/Loss Model  

    params = struct();
    
    % Model identification  
    params.model = 'pt_symmetric';
    params.description = 'PT-symmetric with balanced gain/loss';
    
    % System parameters (SI units: rad/s)
    params.N = 60;                    % System size (sites)
    params.t = 2*pi*10e6;            % Hopping amplitude: t/2π = 10 MHz  
    params.Delta = 2*pi*1e6;         % BdG pairing: Δ/2π = 1 MHz
    
    % PT-breaking threshold (exact finite-size formula)
    params.g_c = 2*params.t * cos(pi/(params.N + 1));
    params.g_c_asymptotic = 2*params.t * (1 - pi^2/(2*params.N^2));
    
    % PT gain/loss parameter (unbroken phase) - CORRECTED ORDER
    params.safety_margin = 0.05;     % 5% safety margin from g_c
    params.delta = params.safety_margin * params.g_c;
    params.g = params.g_c - params.delta;
    
    % CRITICAL FIX: Ensure gN < 1 condition is satisfied
    if params.g * params.N >= 1
        warning('gN = %.1f ≥ 1 violates perturbative condition. Adjusting g to satisfy gN < 1.', params.g * params.N);
        params.g = 0.9 / params.N;  % Set gN = 0.9 < 1
        params.delta = params.g_c - params.g;
        params.safety_margin = params.delta / params.g_c;
        fprintf('Adjusted: g = %.3e, δ = %.3e, safety_margin = %.3f\n', params.g, params.delta, params.safety_margin);
    end
    
    % Additional parameters for multiparameter estimation
    params.mu = 0.1 * params.t;      % Chemical potential: |μ| < 2t
    params.phi = 0.05;               % Peierls phase: |φ| ≪ 1
    
    % Theoretical predictions (using corrected g)
    params.expected_qfi_enhancement = params.t * params.N^2 / (6 * params.delta);
    params.enhancement_factor = params.t * params.N / (6 * params.delta);
    params.heisenberg_advantage = params.enhancement_factor / sqrt(params.N);
    
    % Multiparameter QFI matrix elements (weak-pairing limit)
    params.F_mu_mu = params.N^2 / (4 * params.Delta^2);
    params.F_phi_phi = 3 * params.N^2 * params.t^4 / 2;
    params.F_g_g = params.N^2 * params.Delta^2 / (4 * params.t^2);
    params.F_mu_phi = -params.N^2 * params.Delta * params.t^2 / 4;
    
    % Experimental constraints
    params.coherence_time = 100e-6;  % T₂* = 100 μs
    params.decoherence_rate = 2*pi*50e3; % γ_φ/2π ~ 50 kHz  
    params.pt_breaking_noise = 2*pi*1e3; % γ_PT/2π ~ 1 kHz
    
    % Stability requirements
    params.spectral_stability = params.delta > params.decoherence_rate;
    params.adiabatic_condition = params.delta > 1/(params.t * params.coherence_time);
    
    % Validity checks with warnings instead of errors
    if abs(params.mu) >= 2*params.t
        warning('|μ| = %.2f*t ≥ 2t violates metallic regime condition', abs(params.mu)/params.t);        
    end
    
    if abs(params.phi) >= 1
        warning('|φ| = %.2f ≥ 1 violates perturbative condition', abs(params.phi));
    end
    
    if params.Delta >= params.t
        warning('Δ = %.2f*t ≥ t violates weak-pairing condition', params.Delta/params.t);
    end
    
    if ~params.spectral_stability
        warning('δ = %.2e < decoherence rate %.2e may affect stability', params.delta, params.decoherence_rate);
    end
    
    % Performance metrics
    fprintf('PT Configuration Summary:\n');
    fprintf('g_c = %.3e rad/s (%.2f MHz)\n', params.g_c, params.g_c/(2*pi*1e6));
    fprintf('g = %.3e rad/s (%.2f MHz)\n', params.g, params.g/(2*pi*1e6));
    fprintf('δ = %.3e rad/s (%.2f kHz)\n', params.delta, params.delta/(2*pi*1e3));
    fprintf('gN = %.3f (must be < 1)\n', params.g * params.N);
    fprintf('Enhancement factor η = %.1f\n', params.enhancement_factor);
    fprintf('Heisenberg advantage = %.1f × √N\n', params.heisenberg_advantage);
    
end
