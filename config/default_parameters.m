function params = default_parameters()
% DEFAULT_PARAMETERS  Returns a struct of default simulation parameters
%   Provides consistent values referenced in the manuscript and
%   supplementary material for NHSE, PT, and QFI tests.
%
%   Fields:
%     N_list       - array of system sizes for benchmarking
%     t            - hopping amplitude (Hz units)
%     Delta        - pairing potential (same units as t)
%     mu           - chemical potential (same units as t)
%     gamma_NHSE   - non-Hermitian strength for NHSE model (> t)
%     g_c_PT       - PT-breaking threshold for PT model
%     delta_EP     - detuning from threshold for EP enhancement
%     dparam       - small finite-difference step
%     num_shots    - number of measurements (for QFI bounds)
%
% Example:
%   params = default_parameters();

% System sizes
params.N_list     = [20, 40, 60, 80, 100];

% Hamiltonian parameters
params.t          = 2*pi*10e6;    % 10 MHz hopping
params.Delta      = 2*pi*1e6;     % 1 MHz pairing
params.mu         = 0;            % zero detuning baseline

% NHSE model
params.gamma_NHSE = 1.5 * params.t;   % choose gamma > t

% PT model
% Compute critical threshold g_c = 2t*cos(pi/(N+1)) for largest N
Nmax = max(params.N_list);
params.g_c_PT     = 2 * params.t * cos(pi/(Nmax+1));
params.delta_EP   = 1e-3 * params.t; % small detuning below threshold

% Numerical differentiation step
params.dparam     = 1e-6 * params.t;

% QFI bound sampling
params.num_shots  = 1e5;         % number of repeated measurements

end
