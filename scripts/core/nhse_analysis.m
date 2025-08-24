function results = nhse_analysis()
% NHSE_ANALYSIS  Compute NHSE metrics using precomputed system parameters
%
% Output:
%   results - struct containing nhse_001...nhse_004 with fields:
%       kappa, overlap2, FQ_NHSE, details, regime

% Load core MAT file
S = load('CORE_DATA_MASTER.mat');  % adjust path as needed
systems = S.core_data.nhse_analysis.nhse_systems;
sys_names = fieldnames(systems);

results = struct();

for i = 1:numel(sys_names)
    sys = systems.(sys_names{i});
    
    % Extract required parameters
    params = sys.params;  % N, t, g
    
    N = params.N;
    t = params.t;
    g = params.g;
    
    % Compute localization
    kappa = acosh(abs(g)/t);
    
    % Biorthogonal overlap squared (Eq. S2)
    exp2kN = exp(-2*kappa*N);
    exp2k = exp(-2*kappa);
    num = (1 - exp2k)^2;
    den = (1 - exp2kN)*(exp(2*kappa*N) - 1);
    overlap2 = num / den;
    
    % QFI suppression
    sinhk = sinh(kappa);
    FQ_NHSE = 4*N^3 * exp2kN / (3*t^2 * sinhk^2);
    
    % Save results
    details = struct();
    details.sum_geom = (1 - exp2kN)/(1 - exp2k);
    details.overlap_num = num;
    details.overlap_den = den;
    details.sinhk = sinhk;
    
    results.(sys_names{i}) = struct(...
        'kappa', kappa, ...
        'overlap2', overlap2, ...
        'FQ_NHSE', FQ_NHSE, ...
        'details', details, ...
        'regime', sys.regime ...
    );
end
end
