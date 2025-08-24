function summarize_nhse_results(results)
% SUMMARIZE_NHSE_RESULTS  Summarize and visualize NHSE analysis results
%
% Input:
%   results - struct returned by nhse_analysis()
%
% Example usage:
%   results = nhse_analysis();
%   summarize_nhse_results(results);

% Extract system names
system_names = fieldnames(results);

% Preallocate arrays
Nsystems = numel(system_names);
kappa_vals = zeros(Nsystems,1);
overlap_vals = zeros(Nsystems,1);
FQ_vals = zeros(Nsystems,1);
regimes = strings(Nsystems,1);

% Collect values
for i = 1:Nsystems
    sys = results.(system_names{i});
    kappa_vals(i)   = sys.kappa;
    overlap_vals(i) = sys.overlap2;
    FQ_vals(i)      = sys.FQ_NHSE;
    if isfield(sys, 'regime')
        regimes(i) = sys.regime;
    else
        regimes(i) = "unknown";
    end
end

% Create table
T = table(system_names, kappa_vals, overlap_vals, FQ_vals, regimes, ...
    'VariableNames', {'System', 'Kappa', 'Overlap2', 'FQ_NHSE', 'Regime'});
disp('NHSE Analysis Summary Table:');
disp(T);

% ------------------- Plotting -------------------
figure;

% Bar plots for kappa and QFI
subplot(1,3,1);
bar(kappa_vals);
set(gca, 'XTickLabel', system_names, 'XTickLabelRotation', 45);
ylabel('\kappa (localization)');
title('NHSE Localization');

subplot(1,3,2);
bar(FQ_vals);
set(gca, 'XTickLabel', system_names, 'XTickLabelRotation', 45);
ylabel('F_Q^{NHSE}');
title('QFI Suppression');

% Scatter plot: kappa vs QFI, colored by regime
subplot(1,3,3);
hold on;

colors = containers.Map({'weak','moderate','strong','very strong'}, ...
                        {[0 0.7 0],[0.2 0.2 1],[1 0.3 0],[0.8 0 0.8]});

unique_regimes = unique(regimes);
for r = 1:numel(unique_regimes)
    idx = regimes == unique_regimes(r);
    if colors.isKey(unique_regimes(r))
        c = colors(unique_regimes(r));
    else
        c = [0.5 0.5 0.5]; % default gray for unknown
    end
    scatter(kappa_vals(idx), FQ_vals(idx), 100, c, 'filled', 'DisplayName', unique_regimes(r));
end

text(kappa_vals + 0.02, FQ_vals, system_names, 'FontSize', 10);
xlabel('\kappa');
ylabel('F_Q^{NHSE}');
title('QFI Suppression vs Localization');
legend('Location','best');
grid on;
hold off;

sgtitle('NHSE Analysis Summary');

end
