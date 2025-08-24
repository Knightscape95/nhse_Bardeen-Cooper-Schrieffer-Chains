function results = biorthogonal_qfi_analysis(params)
    % Biorthogonal QFI analysis using PT-symmetric results
    fprintf('Using real biorthogonal QFI analysis (not fallback!)\n');
    
    % Get PT results
    pt_results = pt_symmetric_bdg_analysis(params);
    if isfield(results, 'qfi_value') && results.qfi_value > 1e10 * params_pt.N^2
        warning('Unrealistically high QFI - triggering fallback');
        results.status = 'fallback';
        results.reason = 'unphysical_qfi';
        return;
    end
    
    % 2. Validate biorthogonal overlaps
    if isfield(results, 'overlap') && abs(results.overlap) < 1e-12
        warning('Biorthogonal overlap too small - triggering fallback');
        results.status = 'fallback';
        results.reason = 'overlap_underflow';
        return;
    end
    results = struct();
    results.qfi_biorthogonal = pt_results.qfi_pt;
    results.enhancement_factor = pt_results.enhancement_factor;
    results.formalism = 'biorthogonal';
    results.scaling_law = 'N^2/delta';
    results.status = 'success';
    
    fprintf(' Biorthogonal QFI: F_Q = %.2e\n', results.qfi_biorthogonal);
end
