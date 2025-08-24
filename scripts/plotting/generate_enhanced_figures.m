function generate_enhanced_figures(results, output_dir)
% Generate publication-quality figures for NHSE vs PT enhancement AND multiparameter QFI
    try
       %% FIGURE 1: Non-Hermitian Quantum Metrology Dichotomy
        %% FIGURE 1: Non-Hermitian Quantum Metrology Dichotomy
       %% Parameters
        fig1 = figure('Position', [100, 100, 1200, 500]); 
        t = 1;                    % Hopping amplitude
        gamma = 0.001 * t;        % Non-Hermitian parameter
        
        % Left panel - NHSE suppression parameters (CORRECTED)
        kappa_values = [0.44, 0.62, 0.96];  % CORRECTED: Match PDF values
        N_NHSE = 1:1:120;                    % CORRECTED: Extended range 0-120
        
        % Right panel - PT enhancement parameters (UNCHANGED)
        delta_values = [1e-4, 5e-4, 1e-3];  % Detuning parameters
        N_PT = 10:1:100;                     % System sizes for PT enhancement
        
        %% Create figure with two subplots
        figure('Position', [100, 100, 1200, 500]);
        
        %% Left Panel: NHSE Suppression (CORRECTED TO MATCH PDF)
        subplot(1, 2, 1);
        hold on;
        
        % Colors for different kappa values (CORRECTED for 3 values)
        colors_NHSE = [0.8, 0.2, 0.2; 0.6, 0.4, 0.8; 0.2, 0.6, 0.8];
        line_styles = {'-', '--', ':'};
        
        % Calculate and plot NHSE suppression for each kappa
        for i = 1:length(kappa_values)
            kappa = kappa_values(i);
            
            % NHSE QFI formula: F_Q = 4*N^3*exp(-2*kappa*N)/(3*t^2*sinh(kappa)^2)
            F_NHSE = 4 * N_NHSE.^3 .* exp(-2*kappa*N_NHSE) ./ (3 * t^2 * sinh(kappa)^2);
            
            % Plot analytical curve
            semilogy(N_NHSE, F_NHSE, 'LineWidth', 2.5, 'Color', colors_NHSE(i,:), ...
                'LineStyle', line_styles{i}, 'DisplayName', ...
                sprintf('κ=%.2f', kappa));
        end
        
        % Add SQL baseline
        F_SQL = N_NHSE;
        semilogy(N_NHSE, F_SQL, 'k--', 'LineWidth', 2, 'DisplayName', 'SQL');
        
        % Formatting for left panel (CORRECTED)
        set(gca, 'YScale', 'log');  % Y-axis log, X-axis linear (CORRECTED)
        xlabel('System Size N', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Quantum Fisher Information F_Q', 'FontSize', 12, 'FontWeight', 'bold');
        title('NHSE Suppression (Eq. 17)', 'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        xlim([0, 120]);          % CORRECTED: Match PDF range
        ylim([1e-16, 1e-6]);     % CORRECTED: Match PDF Y-range
        legend('Location', 'northeast', 'FontSize', 10);
        set(gca, 'FontSize', 11);
        
        %% Right Panel: PT-Symmetric Enhancement (CORRECTED RANGES)
        subplot(1, 2, 2);
        hold on;
        
        % Colors for different delta values
        colors_PT = [0.2, 0.8, 0.2; 0.8, 0.6, 0.2; 0.8, 0.2, 0.8];
        
        % Calculate and plot PT enhancement for each delta
        for i = 1:length(delta_values)
            delta = delta_values(i);
            
            % PT QFI formula: F_Q = t*N^2/(6*delta)
            F_PT = t * N_PT.^2 / (6 * delta);
            
            % Plot analytical curve
            loglog(N_PT, F_PT, 'LineWidth', 2.5, 'Color', colors_PT(i,:), ...
                'DisplayName', sprintf('δ=%.0e', delta));
        end
        
        % Add reference lines
        F_SQL_PT = N_PT;
        F_Heisenberg = N_PT.^2;
        
        loglog(N_PT, F_SQL_PT, 'k--', 'LineWidth', 2, 'DisplayName', 'SQL');
        loglog(N_PT, F_Heisenberg, 'k:', 'LineWidth', 2, 'DisplayName', 'Heisenberg N²');
        
        % Formatting for right panel (CORRECTED RANGES)
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('System Size N', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Quantum Fisher Information F_Q', 'FontSize', 12, 'FontWeight', 'bold');
        title('PT Enhancement (Eq. 46)', 'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        xlim([10, 100]);         % CORRECTED: Match PDF range
        ylim([1e2, 1e8]);        % CORRECTED: Match PDF Y-range
        legend('Location', 'northwest', 'FontSize', 10);
        set(gca, 'FontSize', 11);
        
        %% Overall figure formatting
       % sgtitle('Core Dichotomy in Non-Hermitian Quantum Metrology', ...
        %    'FontSize', 16);
        
        % Adjust subplot spacing
        subplot(1,2,1);
        pos1 = get(gca, 'Position');
        pos1(3) = 0.35;  % Adjust width
        set(gca, 'Position', pos1);
        
        subplot(1,2,2);
        pos2(2) = pos2(2) - 0.05;  
        pos2 = get(gca, 'Position');
        pos2(1) = 0.58;  % Adjust left position
        pos2(3) = 0.35;  % Adjust width
        set(gca, 'Position', pos2);
                
        %% FIGURE 2: Multiparameter QFI Matrix Analysis
        fig2 = figure('Position', [150, 150, 1400, 1000]);
        
        % Extended system size range for multiparameter analysis
        N_vals = 10:50;
        
        % Parameters from paper
        Delta = 0.01;  % pairing potential
        t = 1;         % hopping amplitude
        
        % QFI matrix diagonal elements (Eqs. 11 from paper)
        F_mumu = N_vals.^2 ./ (4*Delta^2);
        F_phiphi = 3*N_vals.^2 * t^4 / 2;
        F_gg = N_vals.^2 * Delta^2 / (4*t^2);
        
        % Panel 1: QFI Matrix Diagonal Elements
        subplot(2,2,1);
        loglog(N_vals, F_mumu, 'o-', 'LineWidth', 2, 'DisplayName', 'F_{\mu\mu}');
        hold on
        loglog(N_vals, F_phiphi, 's-', 'LineWidth', 2, 'DisplayName', 'F_{\phi\phi}');
        loglog(N_vals, F_gg, 'd-', 'LineWidth', 2, 'DisplayName', 'F_{gg}');
        loglog(N_vals, N_vals.^2, 'k:', 'LineWidth', 1, 'DisplayName', 'N^2 scaling');
        
        xlabel('System Size N', 'FontSize', 12);
        ylabel('QFI Matrix Diagonal Elements', 'FontSize', 12);
        title('Heisenberg Scaling Verification', 'FontSize', 14);
        legend('Location', 'best', 'FontSize', 10);
        grid on
        
        % Panel 2: Minimum Sensitivities (Cramér-Rao bound)
        subplot(2,2,2);
        % From corrected Eqs. (13-15) in paper
        delta_mu_min = (2*Delta*sqrt(6)) ./ (N_vals * sqrt(6 - Delta^4));
        delta_phi_min = 2 ./ (N_vals * t^2 * sqrt(6 - Delta^4));
        delta_g_min = (2*t) ./ (N_vals * Delta);
        
        loglog(N_vals, delta_mu_min, 'o-', 'LineWidth', 2, 'DisplayName', '\delta\mu (Cramér-Rao)');
        hold on
        loglog(N_vals, delta_phi_min, 's-', 'LineWidth', 2, 'DisplayName', '\delta\phi (Cramér-Rao)');
        loglog(N_vals, delta_g_min, 'd-', 'LineWidth', 2, 'DisplayName', '\deltag (Cramér-Rao)');
        loglog(N_vals, 1./N_vals, 'k:', 'LineWidth', 1, 'DisplayName', '1/N scaling');
        
        xlabel('System Size N', 'FontSize', 12);
        ylabel('Minimum Sensitivity', 'FontSize', 12);
        title('Heisenberg-Limited Sensitivities', 'FontSize', 14);
        legend('Location', 'best', 'FontSize', 10);
        grid on
        
        % Panel 3: Experimental Enhancement Factors
        subplot(2,2,3);
        % Enhancement factors from paper
        eta_mu = 20 * sqrt(N_vals);                    % chemical potential
        eta_phi = t^2 * sqrt(3*N_vals/2);             % Peierls phase  
        eta_g = Delta * sqrt(N_vals) / (2*t);         % gain/loss
        
        plot(N_vals, eta_mu, 'o-', 'LineWidth', 2, 'DisplayName', '\eta_\mu');
        hold on
        plot(N_vals, eta_phi, 's-', 'LineWidth', 2, 'DisplayName', '\eta_\phi');
        plot(N_vals, eta_g, 'd-', 'LineWidth', 2, 'DisplayName', '\eta_g');
        
        xlabel('System Size N', 'FontSize', 12);
        ylabel('Enhancement Factor', 'FontSize', 12);
        title('Experimental Enhancement Factors', 'FontSize', 14);
        legend('Location', 'best', 'FontSize', 10);
        grid on
        
        % Panel 4: QFI Matrix Heatmap (for N=50)
        subplot(2,2,4);
        N_example = 50;
        % Construct QFI matrix from Eq. (11)
        F_matrix = [N_example^2/(4*Delta^2), -N_example^2*Delta*t^2/4, 0;
                   -N_example^2*Delta*t^2/4, 3*N_example^2*t^4/2, 0;
                   0, 0, N_example^2*Delta^2/(4*t^2)];
        
        imagesc(log10(abs(F_matrix)));
        colorbar;
        title(sprintf('QFI Matrix (log_{10}|F|), N=%d', N_example), 'FontSize', 14);
        xlabel('\mu, \phi, g', 'FontSize', 12);
        ylabel('\mu, \phi, g', 'FontSize', 12);
        set(gca, 'XTick', 1:3, 'XTickLabel', {'\mu', '\phi', 'g'});
        set(gca, 'YTick', 1:3, 'YTickLabel', {'\mu', '\phi', 'g'});
        
        sgtitle('Multiparameter Quantum Fisher Information Matrix', 'FontSize', 16);
        
        %% Save both figures
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        
        % Save Figure 1
        exportgraphics(fig1, fullfile(output_dir, 'figure1_dichotomy.png'), 'Resolution', 300);
        savefig(fig1, fullfile(output_dir, 'figure1_dichotomy.fig'));
        
        % Save Figure 2  
        exportgraphics(fig2, fullfile(output_dir, 'figure2_qfi_matrix.png'), 'Resolution', 300);
        savefig(fig2, fullfile(output_dir, 'figure2_qfi_matrix.fig'));
        
        fprintf(' Enhanced figures (Fig 1 and Fig 2) generated successfully\n');
        fprintf('   - Figure 1: NHSE vs PT dichotomy\n');
        fprintf('   - Figure 2: Multiparameter QFI matrix analysis\n');
        
    catch ME
        warning('Enhanced figure generation failed: %s', char(ME.message));
        rethrow(ME);
    end
end
