function [ep_results, T] = exceptional_points(system_key, returnTable)
% EXCEPTIONAL_POINTS  Compute PT-breaking threshold and EP splitting coefficients
%   Uses parameters from CORE_DATA_MASTER.mat (bdg_systems).
%
% USAGE:
%   ep_results = exceptional_points();                   % default: bdg_001_minimal_N6
%   ep_results = exceptional_points('bdg_003_medium_N12'); % single system
%   ep_results = exceptional_points('all');              % all systems
%   [ep_results, T] = exceptional_points('all', true);   % also return comparison table
%
% OUTPUT (single system):
%   ep_results.g_c_exact   - exact threshold
%   ep_results.g_c_asympt  - asymptotic expansion
%   ep_results.alpha       - splitting prefactor
%   ep_results.details     - intermediate values
%
% OUTPUT (all systems):
%   ep_results.<system_key> = struct(...)
%   T = MATLAB table with columns [System, N, g_c_exact, g_c_asympt, alpha]

    % ---------------------------------------------------------------------
    % Load master data
    % ---------------------------------------------------------------------
    datafile = '/MATLAB Drive/ptyy/npj_quantum_matlab/results_20250821_163854/CORE_DATA_MASTER.mat';
    S = load(datafile);

    systems = fieldnames(S.core_data.pt_symmetric_bdg.bdg_systems);

    % ---------------------------------------------------------------------
    % Handle defaults
    % ---------------------------------------------------------------------
    if nargin < 1 || isempty(system_key)
        system_key = 'bdg_001_minimal_N6'; % default
        fprintf('[exceptional_points] No system specified. Using default: %s\n', system_key);
    end
    if nargin < 2
        returnTable = false;
    end

    % ---------------------------------------------------------------------
    % Handle "all" mode
    % ---------------------------------------------------------------------
    if strcmpi(system_key, 'all')
        ep_results = struct();
        rows = [];
        for i = 1:numel(systems)
            key = systems{i};
            res = compute_ep(S.core_data.pt_symmetric_bdg.bdg_systems.(key), key);
            ep_results.(key) = res;
            rows = [rows; {key, res.details.N, res.g_c_exact, res.g_c_asympt, res.alpha}]; %#ok<AGROW>
        end

        if returnTable
            T = cell2table(rows, 'VariableNames', ...
                {'System','N','g_c_exact','g_c_asympt','alpha'});
            disp(T);
        else
            T = [];
        end
        return;
    end

    % ---------------------------------------------------------------------
    % Validate key
    % ---------------------------------------------------------------------
    if ~ismember(system_key, systems)
        error('exceptional_points:invalidKey', ...
              'System key "%s" not found. Available: %s', ...
              system_key, strjoin(systems, ', '));
    end

    % ---------------------------------------------------------------------
    % Single system case
    % ---------------------------------------------------------------------
    ep_results = compute_ep(S.core_data.pt_symmetric_bdg.bdg_systems.(system_key), system_key);
    T = [];

end


% =====================================================================
% Helper: compute EP results for one system
% =====================================================================
function result = compute_ep(sys, system_key)
    params = sys.params;
    N = params.N;
    t = params.t;

    % Validate
    assert(N > 1, 'params.N must exceed 1');
    assert(t > 0, 'params.t must be > 0');

    % Quantized minimum momentum
    k1 = pi/(N+1);

    % Exact threshold
    g_c_exact = 2 * t * cos(k1);

    % Asymptotic expansion
    g_c_asympt = 2 * t * (1 - (k1^2)/2);

    % EP splitting prefactor alpha
    alpha = sqrt((2 * t^2 * sin(k1)^2) / (N+1));

    % Package
    details = struct( ...
        'system_key', system_key, ...
        'N', N, ...
        'k1', k1, ...
        'cos_k1', cos(k1), ...
        'sin_k1', sin(k1), ...
        'Delta', params.Delta, ...
        'mu', params.mu, ...
        'g', params.g ...
    );

    result = struct( ...
        'g_c_exact',  g_c_exact, ...
        'g_c_asympt', g_c_asympt, ...
        'alpha',      alpha, ...
        'details',    details ...
    );
end
