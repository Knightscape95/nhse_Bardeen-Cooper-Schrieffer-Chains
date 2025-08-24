function test_exceptional_points()
%% TEST_EXCEPTIONAL_POINTS Comprehensive test suite for exceptional_points function
%   Tests all functionality including edge cases and error handling

    fprintf('=== Testing exceptional_points function ===\n\n');
    
    % Test counter
    test_count = 0;
    passed = 0;
    
    %% Test 1: Default behavior (no arguments)
    test_count = test_count + 1;
    fprintf('Test %d: Default behavior (no arguments)\n', test_count);
    try
        ep_results = exceptional_points();
        
        % Check output structure
        assert(isstruct(ep_results), 'Output should be a struct');
        assert(isfield(ep_results, 'g_c_exact'), 'Missing g_c_exact field');
        assert(isfield(ep_results, 'g_c_asympt'), 'Missing g_c_asympt field');
        assert(isfield(ep_results, 'alpha'), 'Missing alpha field');
        assert(isfield(ep_results, 'details'), 'Missing details field');
        
        % Check values are reasonable
        assert(ep_results.g_c_exact > 0, 'g_c_exact should be positive');
        assert(ep_results.g_c_asympt > 0, 'g_c_asympt should be positive');
        assert(ep_results.alpha > 0, 'alpha should be positive');
        
        fprintf('  ✓ PASSED - Default case works correctly\n');
        passed = passed + 1;
    catch ME
        fprintf('  ✗ FAILED - %s\n', ME.message);
    end
    
    %% Test 2: Single system case
    test_count = test_count + 1;
    fprintf('\nTest %d: Single system case\n', test_count);
    try
        ep_results = exceptional_points('bdg_002_small_N8');
        
        % Validate against expected values for N=8, t=1
        N = 8;
        t = 1;
        k1_expected = pi/(N+1);
        g_c_exact_expected = 2 * t * cos(k1_expected);
        g_c_asympt_expected = 2 * t * (1 - (k1_expected^2)/2);
        alpha_expected = sqrt((2 * t^2 * sin(k1_expected)^2) / (N+1));
        
        % Check calculations
        tolerance = 1e-12;
        assert(abs(ep_results.g_c_exact - g_c_exact_expected) < tolerance, ...
               'g_c_exact calculation incorrect');
        assert(abs(ep_results.g_c_asympt - g_c_asympt_expected) < tolerance, ...
               'g_c_asympt calculation incorrect');
        assert(abs(ep_results.alpha - alpha_expected) < tolerance, ...
               'alpha calculation incorrect');
        
        % Check details structure
        assert(ep_results.details.N == N, 'N value incorrect in details');
        assert(strcmp(ep_results.details.system_key, 'bdg_002_small_N8'), ...
               'System key incorrect in details');
        
        fprintf('  ✓ PASSED - Single system calculation correct\n');
        passed = passed + 1;
    catch ME
        fprintf('  ✗ FAILED - %s\n', ME.message);
    end
    
    %% Test 3: All systems case without table
    test_count = test_count + 1;
    fprintf('\nTest %d: All systems case (no table)\n', test_count);
    try
        ep_results = exceptional_points('all');
        
        % Check output structure
        assert(isstruct(ep_results), 'Output should be a struct');
        
        % Check that multiple systems are present
        systems = fieldnames(ep_results);
        assert(length(systems) >= 4, 'Should have at least 4 BdG systems');
        
        % Validate each system
        for i = 1:length(systems)
            sys_key = systems{i};
            sys_result = ep_results.(sys_key);
            
            assert(isstruct(sys_result), 'Each system result should be a struct');
            assert(isfield(sys_result, 'g_c_exact'), 'Missing g_c_exact');
            assert(isfield(sys_result, 'g_c_asympt'), 'Missing g_c_asympt');
            assert(isfield(sys_result, 'alpha'), 'Missing alpha');
            
            % Values should be positive and reasonable
            assert(sys_result.g_c_exact > 0 && sys_result.g_c_exact < 2, ...
                   'g_c_exact out of expected range');
            assert(sys_result.alpha > 0 && sys_result.alpha < 1, ...
                   'alpha out of expected range');
        end
        
        fprintf('  ✓ PASSED - All systems processed correctly\n');
        passed = passed + 1;
    catch ME
        fprintf('  ✗ FAILED - %s\n', ME.message);
    end
    
    %% Test 4: All systems case with table
    test_count = test_count + 1;
    fprintf('\nTest %d: All systems case (with table)\n', test_count);
    try
        [ep_results, T] = exceptional_points('all', true);
        
        % Check table structure
        assert(istable(T), 'Second output should be a table');
        assert(ismember('System', T.Properties.VariableNames), 'Missing System column');
        assert(ismember('N', T.Properties.VariableNames), 'Missing N column');
        assert(ismember('g_c_exact', T.Properties.VariableNames), 'Missing g_c_exact column');
        assert(ismember('g_c_asympt', T.Properties.VariableNames), 'Missing g_c_asympt column');
        assert(ismember('alpha', T.Properties.VariableNames), 'Missing alpha column');
        
        % Check table has data
        assert(height(T) >= 4, 'Table should have at least 4 rows');
        
        % Verify ascending trend in g_c values with increasing N
        g_c_vals = T.g_c_exact;
        N_vals = T.N;
        [~, sort_idx] = sort(N_vals);
        sorted_gc = g_c_vals(sort_idx);
        assert(all(diff(sorted_gc) > 0), 'g_c should increase with N');
        
        fprintf('  ✓ PASSED - Table generation and trends correct\n');
        passed = passed + 1;
    catch ME
        fprintf('  ✗ FAILED - %s\n', ME.message);
    end
    
    %% Test 5: Mathematical consistency checks
    test_count = test_count + 1;
    fprintf('\nTest %d: Mathematical consistency\n', test_count);
    try
        % Test specific known values
        ep_N6 = exceptional_points('bdg_001_minimal_N6');
        ep_N16 = exceptional_points('bdg_004_large_N16');
        
        % For larger N, g_c should approach 2*t = 2
        assert(ep_N16.g_c_exact > ep_N6.g_c_exact, ...
               'g_c should increase with N');
        assert(ep_N16.g_c_exact < 2.0, ...
               'g_c should be less than 2*t for finite N');
        
        % Asymptotic should be close to exact for large N
        rel_error_N16 = abs(ep_N16.g_c_exact - ep_N16.g_c_asympt) / ep_N16.g_c_exact;
        rel_error_N6 = abs(ep_N6.g_c_exact - ep_N6.g_c_asympt) / ep_N6.g_c_exact;
        assert(rel_error_N16 < rel_error_N6, ...
               'Asymptotic approximation should improve for larger N');
        
        % Alpha should decrease with N
        assert(ep_N16.alpha < ep_N6.alpha, ...
               'Alpha should decrease with increasing N');
        
        fprintf('  ✓ PASSED - Mathematical relationships verified\n');
        passed = passed + 1;
    catch ME
        fprintf('  ✗ FAILED - %s\n', ME.message);
    end
    
    %% Test 6: Error handling
    test_count = test_count + 1;
    fprintf('\nTest %d: Error handling\n', test_count);
    try
        % Test invalid system key
        try
            exceptional_points('invalid_system');
            assert(false, 'Should have thrown error for invalid system');
        catch ME
            assert(contains(ME.identifier, 'invalidKey'), ...
                   'Should throw invalidKey error');
        end
        
        fprintf('  ✓ PASSED - Error handling works correctly\n');
        passed = passed + 1;
    catch ME
        fprintf('  ✗ FAILED - %s\n', ME.message);
    end
    
    %% Test 7: Numerical accuracy against paper formulas
    test_count = test_count + 1;
    fprintf('\nTest %d: Accuracy against theoretical formulas\n', test_count);
    try
        [all_results, T] = exceptional_points('all', true);
        
        % Check against paper formulas from supplementary Eq. (63) and (71)
        for i = 1:height(T)
            N = T.N(i);
            t = 1.0; % assumed from code
            
            % Paper formula: gc = 2t*cos(π/(N+1))
            k1 = pi/(N+1);
            expected_gc = 2 * t * cos(k1);
            actual_gc = T.g_c_exact(i);
            
            rel_error = abs(actual_gc - expected_gc) / expected_gc;
            assert(rel_error < 1e-14, ...
                   sprintf('Relative error too large for N=%d: %e', N, rel_error));
            
            % Paper formula: α = sqrt(2t²sin²(π/(N+1))/(N+1))
            expected_alpha = sqrt(2 * t^2 * sin(k1)^2 / (N+1));
            actual_alpha = T.alpha(i);
            
            rel_error_alpha = abs(actual_alpha - expected_alpha) / expected_alpha;
            assert(rel_error_alpha < 1e-14, ...
                   sprintf('Alpha relative error too large for N=%d: %e', N, rel_error_alpha));
        end
        
        fprintf('  ✓ PASSED - Numerical accuracy verified\n');
        passed = passed + 1;
    catch ME
        fprintf('  ✗ FAILED - %s\n', ME.message);
    end
    
    %% Summary
    fprintf('\n=== Test Summary ===\n');
    fprintf('Tests run: %d\n', test_count);
    fprintf('Tests passed: %d\n', passed);
    fprintf('Tests failed: %d\n', test_count - passed);
    
    if passed == test_count
        fprintf('🎉 ALL TESTS PASSED! Function is working correctly.\n');
    else
        fprintf('❌ Some tests failed. Please review the function.\n');
    end
    
    %% Display sample output for verification
    fprintf('\n=== Sample Output ===\n');
    [sample_results, sample_table] = exceptional_points('all', true);
    fprintf('Sample table (first few rows):\n');
    if height(sample_table) >= 2
        disp(sample_table(1:min(2, height(sample_table)), :));
    end
    
end

%% Helper function to run the test
% To run this test, simply call: test_exceptional_points()
