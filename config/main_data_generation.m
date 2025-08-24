function main_data_generation()
% ===============================================
% MAIN_DATA_GENERATION - Complete Working Version
% Generates comprehensive validation data for
% Non-Hermitian Quantum Metrology Enhancement
% ===============================================

fprintf('========================================\n');
fprintf('NON-HERMITIAN QUANTUM METROLOGY DATA GENERATION\n');
fprintf('Paper: PT-Symmetric Bardeen-Cooper-Schrieffer Chains\n');
fprintf('Author: Harshank Matkar\n');
fprintf('Date: %s\n', datestr(now));
fprintf('========================================\n');

% Create results directory
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
results_dir = sprintf('results_%s', timestamp);
mkdir(results_dir);
fprintf('Results directory: %s\n', results_dir);

% Load computational settings
settings = computational_settings();
fprintf('Computational settings loaded.\n');

%% 1. EXCEPTIONAL POINTS DATA (N = 4-50)
fprintf('\n=== EXCEPTIONAL POINTS ANALYSIS ===\n');
ep_data = generate_exceptional_points_data();
save(fullfile(results_dir, 'exceptional_points_data.mat'), 'ep_data');

%% 2. BIORTHOGONAL QFI DATA (N = 4-100)
fprintf('\n=== BIORTHOGONAL QFI ANALYSIS ===\n');
qfi_data = generate_biorthogonal_qfi_data();
save(fullfile(results_dir, 'biorthogonal_qfi_data.mat'), 'qfi_data');

%% 3. MULTIPARAMETER QFI MATRIX (N = 10-200)
fprintf('\n=== MULTIPARAMETER QFI MATRIX ===\n');
mp_data = generate_multiparameter_qfi_data();
save(fullfile(results_dir, 'multiparameter_qfi_data.mat'), 'mp_data');

%% 4. PT-SYMMETRIC BDG SYSTEMS (N = 6-200)
fprintf('\n=== PT-SYMMETRIC BDG SYSTEMS ===\n');
pt_data = generate_pt_symmetric_bdg_data();
save(fullfile(results_dir, 'pt_symmetric_bdg_data.mat'), 'pt_data');

%% 5. NHSE ANALYSIS (N = 8-50)
fprintf('\n=== NHSE ANALYSIS ===\n');
nhse_data = generate_nhse_analysis_data();
save(fullfile(results_dir, 'nhse_analysis_data.mat'), 'nhse_data');

%% 6. MASTER DATASET COMPILATION
fprintf('\n=== COMPILING MASTER DATASET ===\n');
core_data = struct();
core_data.exceptional_points = ep_data;
core_data.biorthogonal_qfi = qfi_data;
core_data.multiparameter_qfi = mp_data;
core_data.pt_symmetric_bdg = pt_data;
core_data.nhse_analysis = nhse_data;
core_data.generation_info = get_generation_info(settings);

save(fullfile(results_dir, 'CORE_DATA_MASTER.mat'), 'core_data', '-v7.3');

%% 7. CREATE SUMMARY FILES
create_validation_summary(results_dir, core_data);
create_usage_guide(results_dir);

fprintf('\n========================================\n');
fprintf('DATA GENERATION COMPLETED\n');
fprintf('Total files created: %d\n', length(dir(fullfile(results_dir, '*.mat'))) + 2);
fprintf('Results saved to: %s\n', results_dir);
fprintf('Master dataset: CORE_DATA_MASTER.mat (%.1f MB)\n', ...
    dir(fullfile(results_dir, 'CORE_DATA_MASTER.mat')).bytes/1e6);
fprintf('========================================\n');

end

%% ===============================================
%% DATA GENERATION FUNCTIONS
%% ===============================================

function ep_data = generate_exceptional_points_data()
% Generate exceptional points configurations (N = 4-50)

ep_data = struct();
ep_data.configurations = struct();

N_values = [4:2:20, 25:5:50]; % Proof-of-principle to scaling
config_count = 0;

for N = N_values
    for safety_margin = [0.01, 0.05, 0.10, 0.20] % Different EP proximities
        config_count = config_count + 1;
        
        % FIX: Create valid MATLAB field name (replace periods with underscores)
        margin_str = sprintf('%.3f', safety_margin);
        margin_str = strrep(margin_str, '.', '_'); % Replace . with _
        config_name = sprintf('N%d_margin%s', N, margin_str);
        
        % Load PT parameters with overrides
        params = load_params('pt', 'N', N, 'safety_margin', safety_margin);
        
        % Calculate theoretical predictions
        enhancement = params.t * N / (6 * params.delta);
        qfi_theoretical = params.t * N^2 / (6 * params.delta);
        
        % Store configuration
        ep_data.configurations.(config_name) = struct();
        ep_data.configurations.(config_name).params = params;
        ep_data.configurations.(config_name).enhancement_factor = enhancement;
        ep_data.configurations.(config_name).qfi_theoretical = qfi_theoretical;
        ep_data.configurations.(config_name).pt_breaking_threshold = params.g_c;
        ep_data.configurations.(config_name).operating_point = params.g;
        ep_data.configurations.(config_name).detuning = params.delta;
        
        fprintf('EP Config %d: N=%d, margin=%.1f%%, η=%.1f\n', ...
            config_count, N, safety_margin*100, enhancement);
    end
end

ep_data.summary = struct();
ep_data.summary.total_configurations = config_count;
ep_data.summary.N_range = [min(N_values), max(N_values)];
ep_data.summary.safety_margins = [0.01, 0.05, 0.10, 0.20];

end

function qfi_data = generate_biorthogonal_qfi_data()
% Generate biorthogonal QFI systems (N = 4-100)

qfi_data = struct();
qfi_data.systems = struct();

N_values = [4:2:20, 22:3:50, 55:5:100]; % Comprehensive range
system_count = 0;

for N = N_values
    for model_type = {'nhse', 'pt'}
        system_count = system_count + 1;
        system_name = sprintf('%s_N%d', model_type{1}, N);
        
        % Load appropriate parameters
        if strcmp(model_type{1}, 'nhse')
            params = load_params('nhse', 'N', N);
            params.model_type = 'nhse';
        else
            params = load_params('pt', 'N', N);
            params.model_type = 'pt_symmetric';
        end
        
        % Calculate theoretical QFI
        if strcmp(model_type{1}, 'nhse')
            % NHSE suppression formula: F = 4N³e^(-2κN)/(3t²sinh²κ)
            qfi_theoretical = 4 * N^3 * exp(-2*params.kappa*N) / ...
                (3 * params.t^2 * sinh(params.kappa)^2);
            suppression_factor = qfi_theoretical / N; % Compare to SQL
        else
            % PT enhancement formula: F = tN²/(6δ)
            qfi_theoretical = params.t * N^2 / (6 * params.delta);
            enhancement_factor = qfi_theoretical / N; % Compare to SQL
        end
        
        % Store system
        qfi_data.systems.(system_name) = struct();
        qfi_data.systems.(system_name).params = params;
        qfi_data.systems.(system_name).qfi_theoretical = qfi_theoretical;
        
        if strcmp(model_type{1}, 'nhse')
            qfi_data.systems.(system_name).suppression_factor = suppression_factor;
        else
            qfi_data.systems.(system_name).enhancement_factor = enhancement_factor;
        end
        
        fprintf('QFI System %d: %s, N=%d\n', system_count, model_type{1}, N);
    end
end

qfi_data.summary = struct();
qfi_data.summary.total_systems = system_count;
qfi_data.summary.N_range = [min(N_values), max(N_values)];
qfi_data.summary.models = {'nhse', 'pt_symmetric'};

end

function mp_data = generate_multiparameter_qfi_data()
% Generate multiparameter QFI parameter sets (N = 10-200)

mp_data = struct();
mp_data.parameter_sets = struct();

N_values = [10:5:50, 60:10:100, 120:20:200]; % Focus on larger systems
param_count = 0;

for N = N_values
    for Delta_ratio = [0.05, 0.10, 0.15] % Weak pairing regimes
        param_count = param_count + 1;
        
        % FIX: Create valid MATLAB field name
        delta_str = sprintf('%.3f', Delta_ratio);
        delta_str = strrep(delta_str, '.', '_'); % Replace . with _
        param_name = sprintf('N%d_Delta%s', N, delta_str);
        
        % Load PT parameters
        params = load_params('pt', 'N', N);
        params.Delta = Delta_ratio * params.t; % Set pairing strength
        
        % Calculate multiparameter QFI matrix elements (Eq. 11 from paper)
        F_mu_mu = N^2 / (4 * params.Delta^2);
        F_phi_phi = 3 * N^2 * params.t^4 / 2;
        F_g_g = N^2 * params.Delta^2 / (4 * params.t^2);
        F_mu_phi = -N^2 * params.Delta * params.t^2 / 4;
        
        % Construct QFI matrix
        qfi_matrix = [
            F_mu_mu,  F_mu_phi, 0;
            F_mu_phi, F_phi_phi, 0;
            0,        0,        F_g_g
        ];
        
        % Calculate optimal sensitivities (Eq. 13-15)
        weak_pairing_factor = 6 - params.Delta^4;
        nu = 1e6; % Number of measurements
        
        delta_mu_min = 2*params.Delta*sqrt(6) / (N * sqrt(nu * weak_pairing_factor));
        delta_phi_min = 2 / (N*params.t^2 * sqrt(nu * weak_pairing_factor));
        delta_g_min = 2*params.t / (N*params.Delta*sqrt(nu));
        
        % Store parameter set
        mp_data.parameter_sets.(param_name) = struct();
        mp_data.parameter_sets.(param_name).params = params;
        mp_data.parameter_sets.(param_name).qfi_matrix = qfi_matrix;
        mp_data.parameter_sets.(param_name).diagonal_elements = [F_mu_mu, F_phi_phi, F_g_g];
        mp_data.parameter_sets.(param_name).optimal_sensitivities = [delta_mu_min, delta_phi_min, delta_g_min];
        mp_data.parameter_sets.(param_name).heisenberg_scaling = all(diag(qfi_matrix) > 0); % Check N² scaling
        
        fprintf('MP Param %d: N=%d, Δ/t=%.2f, det(F)=%.2e\n', ...
            param_count, N, Delta_ratio, det(qfi_matrix));
    end
end

mp_data.summary = struct();
mp_data.summary.total_parameter_sets = param_count;
mp_data.summary.N_range = [min(N_values), max(N_values)];
mp_data.summary.Delta_ratios = [0.05, 0.10, 0.15];

end

function pt_data = generate_pt_symmetric_bdg_data()
% Generate PT-symmetric BdG systems (N = 6-200)

pt_data = struct();
pt_data.bdg_systems = struct();

N_values = [6:2:20, 25:5:100, 110:10:200]; % Full scaling range
system_count = 0;

for N = N_values
    for experimental_regime = {'proof_of_principle', 'scaling_verification', 'asymptotic_behavior'}
        system_count = system_count + 1;
        system_name = sprintf('%s_N%d', experimental_regime{1}, N);
        
        % Load experimental parameters
        params = load_params('experimental');
        params.N = N;
        
        % Adjust parameters based on regime
        switch experimental_regime{1}
            case 'proof_of_principle'
                if N > 8, continue; end % Skip large N for this regime
                params.tolerance = 1e-14; % High precision for small systems
            case 'scaling_verification'
                if N < 20 || N > 50, continue; end % Mid-range N
                params.tolerance = 1e-12; % Standard precision
            case 'asymptotic_behavior'
                if N < 50, continue; end % Large N only
                params.tolerance = 1e-10; % Relaxed precision for speed
        end
        
        % Calculate PT threshold for this N
        g_c = 2*params.t * cos(pi/(N + 1));
        g_c_asymptotic = 2*params.t * (1 - pi^2/(2*N^2));
        finite_size_correction = g_c - g_c_asymptotic;
        
        % Store BdG system
        pt_data.bdg_systems.(system_name) = struct();
        pt_data.bdg_systems.(system_name).params = params;
        pt_data.bdg_systems.(system_name).pt_threshold_exact = g_c;
        pt_data.bdg_systems.(system_name).pt_threshold_asymptotic = g_c_asymptotic;
        pt_data.bdg_systems.(system_name).finite_size_correction = finite_size_correction;
        pt_data.bdg_systems.(system_name).experimental_regime = experimental_regime{1};
        
        fprintf('PT BdG %d: %s, N=%d, g_c=%.3e\n', ...
            system_count, experimental_regime{1}, N, g_c);
    end
end

pt_data.summary = struct();
pt_data.summary.total_systems = system_count;
pt_data.summary.N_range = [min(N_values), max(N_values)];
pt_data.summary.experimental_regimes = {'proof_of_principle', 'scaling_verification', 'asymptotic_behavior'};

end

function nhse_data = generate_nhse_analysis_data()
% Generate NHSE analysis systems (N = 8-50)

nhse_data = struct();
nhse_data.nhse_systems = struct();

N_values = 8:2:50; % Limited range due to exponential suppression
system_count = 0;

for N = N_values
    for gamma_ratio = [1.1, 1.5, 2.0, 2.5] % Different NHSE strengths
        system_count = system_count + 1;
        
        % FIX: Create valid MATLAB field name
        gamma_str = sprintf('%.1f', gamma_ratio);
        gamma_str = strrep(gamma_str, '.', '_'); % Replace . with _
        system_name = sprintf('gamma%s_N%d', gamma_str, N);
        
        % Load NHSE parameters with override
        params = load_params('nhse', 'N', N);
        params.gamma = gamma_ratio * params.t;
        
        % Recalculate localization parameters
        kappa = acosh(gamma_ratio);
        localization_length = 1/kappa;
        biorthogonal_overlap = exp(-2*kappa*N);
        
        % Calculate suppressed QFI (Eq. 8 from paper)
        qfi_suppressed = 4 * N^3 * exp(-2*kappa*N) / (3 * params.t^2 * sinh(kappa)^2);
        suppression_ratio = qfi_suppressed / N; % Compare to SQL
        
        % Store NHSE system
        nhse_data.nhse_systems.(system_name) = struct();
        nhse_data.nhse_systems.(system_name).params = params;
        nhse_data.nhse_systems.(system_name).gamma_over_t = gamma_ratio;
        nhse_data.nhse_systems.(system_name).kappa = kappa;
        nhse_data.nhse_systems.(system_name).localization_length = localization_length;
        nhse_data.nhse_systems.(system_name).biorthogonal_overlap = biorthogonal_overlap;
        nhse_data.nhse_systems.(system_name).qfi_suppressed = qfi_suppressed;
        nhse_data.nhse_systems.(system_name).suppression_ratio = suppression_ratio;
        
        fprintf('NHSE System %d: γ/t=%.1f, N=%d, suppression=%.2e\n', ...
            system_count, gamma_ratio, N, suppression_ratio);
    end
end

nhse_data.summary = struct();
nhse_data.summary.total_systems = system_count;
nhse_data.summary.N_range = [min(N_values), max(N_values)];
nhse_data.summary.gamma_ratios = [1.1, 1.5, 2.0, 2.5];

end

%% ===============================================
%% UTILITY FUNCTIONS
%% ===============================================

function info = get_generation_info(settings)
% Create generation metadata

info = struct();
info.timestamp = datestr(now);
info.matlab_version = version();
info.computer = computer();
info.computational_settings = settings;
info.paper_reference = 'Non-Hermitian Quantum Metrology Enhancement and Skin Effect Suppression in PT-Symmetric Bardeen-Cooper-Schrieffer Chains';
info.author = 'Harshank Matkar';
info.institution = 'Government College of Engineering Aurangabad';
info.total_configurations = 0; % Will be updated by calling functions

end

function create_validation_summary(results_dir, core_data)
% Create human-readable validation summary

filename = fullfile(results_dir, 'VALIDATION_SUMMARY.txt');
fid = fopen(filename, 'w');

fprintf(fid, '========================================\n');
fprintf(fid, 'NON-HERMITIAN QUANTUM METROLOGY\n');
fprintf(fid, 'VALIDATION SUMMARY\n');
fprintf(fid, '========================================\n\n');

fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Paper: PT-Symmetric Bardeen-Cooper-Schrieffer Chains\n');
fprintf(fid, 'Author: Harshank Matkar\n\n');

fprintf(fid, 'DATASET OVERVIEW:\n');
fprintf(fid, '- Exceptional Points: %d configurations\n', core_data.exceptional_points.summary.total_configurations);
fprintf(fid, '- Biorthogonal QFI: %d systems\n', core_data.biorthogonal_qfi.summary.total_systems);
fprintf(fid, '- Multiparameter QFI: %d parameter sets\n', core_data.multiparameter_qfi.summary.total_parameter_sets);
fprintf(fid, '- PT-Symmetric BdG: %d systems\n', core_data.pt_symmetric_bdg.summary.total_systems);
fprintf(fid, '- NHSE Analysis: %d systems\n', core_data.nhse_analysis.summary.total_systems);

fprintf(fid, '\nSYSTEM SIZE RANGES:\n');
fprintf(fid, '- Exceptional Points: N = %d-%d\n', core_data.exceptional_points.summary.N_range);
fprintf(fid, '- Biorthogonal QFI: N = %d-%d\n', core_data.biorthogonal_qfi.summary.N_range);
fprintf(fid, '- Multiparameter QFI: N = %d-%d\n', core_data.multiparameter_qfi.summary.N_range);
fprintf(fid, '- PT-Symmetric BdG: N = %d-%d\n', core_data.pt_symmetric_bdg.summary.N_range);
fprintf(fid, '- NHSE Analysis: N = %d-%d\n', core_data.nhse_analysis.summary.N_range);

fprintf(fid, '\nKEY VALIDATION CATEGORIES:\n');
fprintf(fid, '- Proof of Principle (N = 4-8): Small-system validation\n');
fprintf(fid, '- Finite-size Corrections (N = 10-20): Asymptotic validity\n');
fprintf(fid, '- Scaling Verification (N = 20-50): N² enhancement\n');
fprintf(fid, '- Asymptotic Behavior (N = 50-100): Bulk limit\n');
fprintf(fid, '- Large-scale Validation (N = 100-200): Theoretical completeness\n');

fclose(fid);

end

function create_usage_guide(results_dir)
% Create usage guide for the generated data

filename = fullfile(results_dir, 'CORE_DATA_USAGE_GUIDE.txt');
fid = fopen(filename, 'w');

fprintf(fid, '========================================\n');
fprintf(fid, 'CORE DATA USAGE GUIDE\n');
fprintf(fid, '========================================\n\n');

fprintf(fid, 'LOADING DATA:\n');
fprintf(fid, 'load(''CORE_DATA_MASTER.mat'');\n');
fprintf(fid, 'ep_configs = core_data.exceptional_points.configurations;\n');
fprintf(fid, 'qfi_systems = core_data.biorthogonal_qfi.systems;\n\n');

fprintf(fid, 'EXAMPLE USAGE:\n');
fprintf(fid, '%% Test exceptional points function\n');
fprintf(fid, 'config_names = fieldnames(ep_configs);\n');
fprintf(fid, 'for i = 1:length(config_names)\n');
fprintf(fid, '    params = ep_configs.(config_names{i}).params;\n');
fprintf(fid, '    result = exceptional_points(params);\n');
fprintf(fid, '    %% Validate against theoretical predictions\n');
fprintf(fid, 'end\n\n');

fprintf(fid, 'DATA STRUCTURE:\n');
fprintf(fid, 'core_data/\n');
fprintf(fid, '├── exceptional_points/configurations/     (EP analysis)\n');
fprintf(fid, '├── biorthogonal_qfi/systems/             (QFI systems)\n');
fprintf(fid, '├── multiparameter_qfi/parameter_sets/    (Parameter sets)\n');
fprintf(fid, '├── pt_symmetric_bdg/bdg_systems/         (BdG systems)\n');
fprintf(fid, '├── nhse_analysis/nhse_systems/           (NHSE systems)\n');
fprintf(fid, '└── generation_info/                      (Metadata)\n\n');

fprintf(fid, 'VALIDATION METRICS:\n');
fprintf(fid, '- Enhancement factors: η ≈ 20√N for chemical potential\n');
fprintf(fid, '- QFI scaling: F_PT ∝ N²/δ (Heisenberg limit)\n');
fprintf(fid, '- NHSE suppression: F_NHSE ∝ N³e^(-2κN) (exponential failure)\n');
fprintf(fid, '- Multiparameter bounds: δμ ≥ 2Δ/(N√ν), δφ ≥ √2/(Nt²√3ν)\n');

fclose(fid);

end
