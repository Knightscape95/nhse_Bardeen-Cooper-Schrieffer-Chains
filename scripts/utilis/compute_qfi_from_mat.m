function [QFI_nhse, QFI_pt] = compute_qfi_from_mat(matfile)
% COMPUTE_QFI_FROM_MAT Load experimental results and compute biorthogonal QFI
% 
% Inputs:
%   matfile - path to .mat file containing results_nhse and results_pt
%
% Outputs:
%   QFI_nhse - QFI for NHSE case
%   QFI_pt   - QFI for PT-symmetric case

    % Load the mat file
    S = load("/MATLAB Drive/ptyy/experimental_results.mat");

    % --- NHSE ---
    if ~isfield(S, 'results_nhse')
        error('MAT file must contain results_nhse struct');
    end
    nhse = S.results_nhse;

    % Fill missing fields
    if ~isfield(nhse, 'H')
        if isfield(S, 'H1'), nhse.H = S.H1; else, error('H1 not found for NHSE'); end
    end
    if ~isfield(nhse, 'L')
        nhse.L = inv(nhse.V);  % compute left eigenvectors
    end
    if ~isfield(nhse, 'd')
        nhse.d = diag(nhse.V' * nhse.H * nhse.V); % fallback, just in case
    end

    % Compute QFI for NHSE
    fprintf('Computing QFI for NHSE...\n');
    QFI_nhse = biorthogonal_qfi(nhse.H, nhse.H, nhse.d, nhse.V, nhse.L); % H passed twice if params not needed

    % --- PT-SYMMETRIC ---
    if ~isfield(S, 'results_pt')
        error('MAT file must contain results_pt struct');
    end
    pt = S.results_pt;

    % Fill missing fields
    if ~isfield(pt, 'H')
        if isfield(S, 'H2'), pt.H = S.H2; else, error('H2 not found for PT'); end
    end
    if ~isfield(pt, 'L')
        pt.L = inv(pt.V);
    end
    if ~isfield(pt, 'd')
        pt.d = diag(pt.V' * pt.H * pt.V);
    end

    % Compute QFI for PT-symmetric
    fprintf('Computing QFI for PT-symmetric...\n');
    QFI_pt = biorthogonal_qfi(pt.H, pt.H, pt.d, pt.V, pt.L);

    fprintf('Done. QFI computed for both NHSE and PT-symmetric cases.\n');
end
