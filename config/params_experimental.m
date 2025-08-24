function params = params_experimental()
% Experimental Configuration - Realistic superconducting circuit parameters
% Based on manuscript Sec. VI: Platform requirements and validation strategy

params = struct();

% Model identification
params.model = 'experimental';
params.description = 'Realistic superconducting circuit parameters';

% Superconducting circuit parameters (from manuscript)
params.t_freq = 1e6; % t/2π = 10 MHz (inductive coupling)
params.Delta_freq = 1e6; % Δ/2π = 1 MHz (pairing)
params.g_freq = 1e5; % g/2π = 1 MHz (pump-driven)

% Convert to angular frequencies
params.t = 2*pi * params.t_freq;
params.Delta = 2*pi * params.Delta_freq;
params.g = 2*pi * params.g_freq;

% System sizes for different applications
params.N_proof_of_principle = 8; % Small-system validation
params.N_scaling_verification = 50; % Scaling studies
params.N_modular_limit = 100; % Modular chip scalability

% FIX: Calculate PT-breaking threshold and operating parameters
N = params.N_scaling_verification; % Use scaling verification size as default
params.g_c = 2*params.t * cos(pi/(N + 1)); % PT-breaking threshold
params.safety_margin = 0.05; % 5% safety margin from g_c
params.delta = params.safety_margin * params.g_c; % Detuning from threshold
params.g_operating = params.g_c - params.delta; % Actual operating point

% Ensure gN < 1 constraint for perturbative validity
if params.g_operating * N >= 1
    warning('gN = %.1f ≥ 1 violates perturbative condition. Adjusting g to satisfy gN < 1.', params.g_operating * N);
    params.g_operating = 0.9 / N;
    params.delta = params.g_c - params.g_operating;
    params.safety_margin = params.delta / params.g_c;
end

% Noise and decoherence (realistic values)
params.gamma_phi = 2*pi * 50e3; % Dephasing: γ_φ/2π ~ 50 kHz
params.gamma_minus = 2*pi * 10e3; % Amplitude damping: γ_-/2π ~ 10 kHz
params.gamma_pt = 2*pi * 1e3; % PT-breaking noise: γ_PT/2π ~ 1 kHz
params.T2_star = 100e-6; % Coherence time: T₂* = 100 μs

% Control and readout
params.readout_fidelity = 0.99; % Single-shot readout fidelity
params.gate_fidelity = 0.999; % Two-qubit gate fidelity
params.initialization_fidelity = 0.995; % State preparation fidelity

% Performance predictions (from manuscript Eq. 13-15)
nu = 1e6; % Number of measurements

% Enhancement factors over classical sensing
params.eta_mu = (2*params.t/params.Delta) * sqrt(N); % ≈ 20√N
params.eta_phi = (params.t^2 * sqrt(3*N)) / 2; % ≈ t²√(3N/2)
params.eta_g = (params.Delta/params.t) * sqrt(N); % ≈ (Δ/t)√N

% Minimum achievable sensitivities (normalized)
weak_pairing_factor = 6 - params.Delta^4;
params.delta_mu_achievable = 2*params.Delta*sqrt(6) / (N * sqrt(nu * weak_pairing_factor));
params.delta_phi_achievable = 2 / (N * params.t^2 * sqrt(nu * weak_pairing_factor));
params.delta_g_achievable = 2*params.t / (N * params.Delta * sqrt(nu));

% Systematic error contributions (from manuscript Table I)
params.calibration_drift = 1e-4; % Δg_cal/g ~ 10⁻⁴
params.temperature_fluctuations = 1e-3; % ΔT/T ~ 10⁻³
params.magnetic_noise = 1e-6; % ΔB ~ 1 μT (Tesla)
params.readout_nonlinearity = 1e-3; % ε_NL ~ 10⁻³

% Resource requirements (NOW USING CORRECT params.delta)
params.initialization_time = 10 / params.t; % T_prep = 10/t
params.ramping_time = 1 / (params.t * params.delta); % T_ramp = ℏ/(tδ)
params.interaction_time = sqrt(nu/params.t) * sqrt(N^2/params.delta); % T_int
params.readout_time_per_site = N / params.t; % T_read = N/t_read

% Validation against manuscript benchmarks
fprintf('Experimental Configuration Validation:\n');
fprintf('g_c = %.3e rad/s (%.2f MHz)\n', params.g_c, params.g_c/(2*pi*1e6));
fprintf('g_operating = %.3e rad/s (%.2f MHz)\n', params.g_operating, params.g_operating/(2*pi*1e6));
fprintf('δ = %.3e rad/s (%.2f kHz)\n', params.delta, params.delta/(2*pi*1e3));
fprintf('gN = %.3f (must be < 1)\n', params.g_operating * N);
fprintf('Enhancement factors: η_μ = %.1f, η_φ = %.1f, η_g = %.1f\n', ...
    params.eta_mu, params.eta_phi, params.eta_g);
fprintf('Achievable sensitivities: δμ = %.2e, δφ = %.2e, δg = %.2e\n', ...
    params.delta_mu_achievable, params.delta_phi_achievable, params.delta_g_achievable);

end
