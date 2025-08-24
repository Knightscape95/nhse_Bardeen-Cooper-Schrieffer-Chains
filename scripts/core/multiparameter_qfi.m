function [F, F_inv, components] = multiparameter_qfi(params, mode)
% MULTIPARAMETER_QFI - Rigorous multiparameter QFI matrix
%
%   [F, Finv, comp] = multiparameter_qfi(params, mode)
%
%   Inputs:
%       params - struct with fields: N, t, Delta, mu, g
%       mode   - 'analytical' (default) or 'numerical'
%
%   If no params provided, loads defaults from CORE_DATA_MASTER.mat
%   in the most recent results_* folder in the project.

    % ---------------------------------------------------------------------
    % Defaults (load from master_core_data if no params supplied)
    % ---------------------------------------------------------------------
    if nargin < 1 || isempty(params)

        % Starting directory = where this function lives
        thisDir = fileparts(mfilename('fullpath'));

        % Candidate parent levels to search
        searchDirs = { ...
            fullfile(thisDir, '..'), ...
            fullfile(thisDir, '..', '..'), ...
            fullfile(thisDir, '..', '..', '..')};

        datafile = '';
        for s = 1:numel(searchDirs)
            d = dir(fullfile(searchDirs{s}, 'results_*'));
            d = d([d.isdir]);
            if ~isempty(d)
                [~, idx] = max([d.datenum]); % pick most recent
                candidate = fullfile(searchDirs{s}, d(idx).name, 'CORE_DATA_MASTER.mat');
                if isfile(candidate)
                    datafile = candidate;
                    break;
                end
            end
        end

        % If still not found, do a global search for safety
        if isempty(datafile)
            d = dir(fullfile(searchDirs{end}, '**', 'CORE_DATA_MASTER.mat'));
            if isempty(d)
                error('multiparameter_qfi:missingData', ...
                    'CORE_DATA_MASTER.mat not found in project root or results_* folders.');
            end
            [~, idx] = max([d.datenum]);
            datafile = fullfile(d(idx).folder, d(idx).name);
        end

        fprintf('[multiparameter_qfi] Loading defaults from: %s\n', datafile);
        S = load(datafile, 'core_data');

        % Extract default params (minimal PT-symmetric system)
        params = S.core_data.pt_symmetric_bdg.bdg_systems.bdg_001_minimal_N6.params;
        fprintf('[multiparameter_qfi] Using default params (N=%d)\n', params.N);
    end

    if nargin < 2 || isempty(mode)
        mode = 'analytical';  % default
    end

    % ---------------------------------------------------------------------
    % Validate inputs
    % ---------------------------------------------------------------------
    if ~(isfield(params,'N') && params.N >= 2)
        error('multiparameter_qfi:invalidN', 'params.N required, N >= 2');
    end
    required = {'t','Delta','mu','g'};
    for k = 1:numel(required)
        if ~isfield(params, required{k})
            error('multiparameter_qfi:missingField', 'params.%s required', required{k});
        end
    end

    % Extract
    N      = params.N;
    t      = params.t;
    Delta  = params.Delta;
    mu     = params.mu; %#ok<NASGU>
    g      = params.g;  %#ok<NASGU>

    % ---------------------------------------------------------------------
    % Construct QFI matrix (from Supplement formulas)
    % ---------------------------------------------------------------------
    F_mu_mu    = N^2/(4*Delta^2);
    F_phi_phi  = (3/2)*N^2*t^4;
    F_gg       = N^2*Delta^2/(4*t^2);
    F_mu_phi   = -N^2*Delta*t^2/4;

    F = [F_mu_mu, F_mu_phi, 0;
         F_mu_phi, F_phi_phi, 0;
         0,        0,        F_gg];

    % ---------------------------------------------------------------------
    % Inversion strategy
    % ---------------------------------------------------------------------
    switch lower(mode)
        case 'analytical'
            det_block = F_mu_mu*F_phi_phi - F_mu_phi^2;
            if abs(det_block) < 1e-14
                warning('Analytical inverse unstable, falling back to numerical');
                F_inv = pinv(F);
            else
                inv_block = (1/det_block) * [F_phi_phi, -F_mu_phi; -F_mu_phi, F_mu_mu];
                F_inv = blkdiag(inv_block, 1/F_gg);
            end
        case 'numerical'
            try
                F_inv = inv(F);
            catch
                F_inv = pinv(F);
            end
        otherwise
            error('multiparameter_qfi:badMode', 'mode must be ''analytical'' or ''numerical''');
    end

    % ---------------------------------------------------------------------
    % Consistency check
    % ---------------------------------------------------------------------
    I = eye(size(F));
    err_rel = norm(F*F_inv - I, 'fro') / norm(I, 'fro');

    % ---------------------------------------------------------------------
    % Package components
    % ---------------------------------------------------------------------
    components = struct();
    components.F_mu_mu   = F_mu_mu;
    components.F_phi_phi = F_phi_phi;
    components.F_gg      = F_gg;
    components.F_mu_phi  = F_mu_phi;
    components.err_rel   = err_rel;
    components.mode      = mode;
    components.labels    = {'mu','phi','g'};

end
