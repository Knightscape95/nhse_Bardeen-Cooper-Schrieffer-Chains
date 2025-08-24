function params = load_params(model_type, varargin)

% Universal parameter loader for all configurations
% Usage: params = load_params('nhse') or load_params('pt') etc.
% CORRECTED: Robust path handling and enforcing gN < 1 constraint after overrides

try
    config_path = fileparts(mfilename('fullpath'));
catch
    config_path = pwd; % fallback to current directory
end

% Add config folder to path if not already present
if ~contains(path, config_path)
    addpath(config_path);
end

% Load specified configuration
switch lower(model_type)
    case {'nhse', 'skin_effect'}
        params = params_nhse();
    case {'pt', 'pt_symmetric', 'exceptional_point'}
        params = params_pt();
    case {'validation', 'test', 'benchmark'}
        params = params_validation();
    case {'experimental', 'realistic', 'superconducting'}
        params = params_experimental();
    otherwise
        error('Unknown model type: %s. Available: nhse, pt, validation, experimental', model_type);
end

% Apply any parameter overrides
if nargin > 1
    for i = 1:2:length(varargin)
        if isfield(params, varargin{i})
            old_val = params.(varargin{i});
            params.(varargin{i}) = varargin{i+1};
            fprintf('Override: %s = %g (was %g)\n', varargin{i}, varargin{i+1}, old_val);

            % Recalculate dependent parameters if needed
            if strcmp(varargin{i}, 'N') && strcmp(params.model, 'pt_symmetric')
                % Recalculate PT threshold for new N
                params.g_c = 2*params.t * cos(pi/(params.N + 1));
                params.delta = params.safety_margin * params.g_c;
                params.g = params.g_c - params.delta;

                % Enforce gN < 1
                if params.g * params.N >= 1
                    warning('gN = %.1f ≥ 1 violates perturbative condition. Adjusting g to satisfy gN < 1.', params.g * params.N);
                    params.g = 0.9 / params.N;
                    params.delta = params.g_c - params.g;
                    params.safety_margin = params.delta / params.g_c;
                end
                fprintf('Recalculated: g_c = %.3e, g = %.3e, delta = %.3e\n', params.g_c, params.g, params.delta);
            end

            if strcmp(varargin{i}, 'safety_margin') && strcmp(params.model, 'pt_symmetric')
                % Recalculate delta and g for new safety margin
                params.delta = params.safety_margin * params.g_c;
                params.g = params.g_c - params.delta;

                % Enforce gN < 1
                if params.g * params.N >= 1
                    warning('gN = %.1f ≥ 1 violates perturbative condition. Adjusting g to satisfy gN < 1.', params.g * params.N);
                    params.g = 0.9 / params.N;
                    params.delta = params.g_c - params.g;
                    params.safety_margin = params.delta / params.g_c;
                end

                fprintf('Recalculated: delta = %.3e, g = %.3e, safety_margin = %.3f\n', params.delta, params.g, params.safety_margin);
            end
        else
            warning('Parameter %s not found in configuration', varargin{i});
        end
    end
end

% Add metadata
params.config_loaded = datetime('now');
params.matlab_version = version();
params.config_path = config_path;

end
