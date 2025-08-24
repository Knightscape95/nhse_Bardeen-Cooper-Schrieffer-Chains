function params = params_validation()

% Validation Configuration - Small system for testing

% Based on manuscript Sec. VI.B: N = 4-8 for proof-of-principle

params = struct();

% Model identification
params.model = 'validation';
params.description = 'Small-system validation and benchmarking';

% Validation parameters
params.N_test = [4, 6, 8, 10]; % Test system sizes
params.t = 1.0; % Normalized units (t = 1)
params.Delta = 0.1; % Weak pairing: Δ/t = 0.1

% NHSE validation
params.gamma_nhse = [0.5, 0.8, 0.95] * params.t;
params.kappa_values = arrayfun(@(g) acosh(abs(g)/params.t), params.gamma_nhse);

% PT validation
params.N_pt = 10;
params.g_c_validation = 2*params.t * cos(pi/(params.N_pt + 1));
params.delta_values = [1e-4, 5e-4, 1e-3, 5e-3];
params.g_pt_values = params.g_c_validation - params.delta_values;

% Numerical precision requirements
params.tolerance = 1e-12; % Numerical accuracy target
params.max_iterations = 50; % Convergence limit
params.relative_error_threshold = 1e-10;

% Benchmarking parameters
params.timing_sizes = [100, 1000, 5000, 10000];
params.memory_threshold = 1e9; % 1 GB memory limit
params.parallel_cores = 4; % For distributed computation

% --- Additional test for gN violation in PT model (Optional) ---
params.test_gn_violation = true;
params.N_gn_test = 8;
params.t_gn_test = 1.0; % Consistent with normalized units
params.g_c_gn_test = 2 * params.t_gn_test * cos(pi/(params.N_gn_test + 1));
params.safety_margin_gn_test = 0.01; % 1% safety margin
params.delta_gn_test = params.safety_margin_gn_test * params.g_c_gn_test;
params.g_gn_test = params.g_c_gn_test - params.delta_gn_test;

fprintf('GN Violation Test Configuration:\n');
fprintf('N = %d, g_c = %.3f, g = %.3f, δ = %.3f\n', params.N_gn_test, params.g_c_gn_test, params.g_gn_test, params.delta_gn_test);
fprintf('gN = %.3f (violates gN < 1 constraint)\n', params.g_gn_test * params.N_gn_test);
fprintf('Expected enhancement η = %.1f\n', params.t_gn_test * params.N_gn_test / (6 * params.delta_gn_test));

end
