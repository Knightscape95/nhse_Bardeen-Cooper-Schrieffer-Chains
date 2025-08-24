classdef test_multiparameter_qfi < matlab.unittest.TestCase
    % TEST_MULTIPARAMETER_QFI
    % Unit tests for multiparameter_qfi.m
    %
    % Runs in 'numerical' mode by default (stable).
    % To check analytical formulas: set `useAnalytical = true` below.

    properties (Constant)
        % Switch between 'numerical' (default) and 'analytical'
        useAnalytical = false;
    end

    methods (Test)
        function testDefaultDemo(testCase)
            fprintf('Running multiparameter_qfi with manuscript defaults...\n');
            [F, Finv, comp] = multiparameter_qfi([], ...
                ternary(testCase.useAnalytical,'analytical','numerical'));
            
            % Dimension check
            testCase.verifySize(F,[3 3]);
            testCase.verifySize(Finv,[3 3]);

            % Inverse consistency check
            testCase.verifyLessThanOrEqual(comp.err_rel,0.2);
        end

        function testConsistencyIdentity(testCase)
            params = struct('N',6,'t',1,'Delta',0.05,'mu',0.02,'g',0.01);
            [F, Finv, comp] = multiparameter_qfi(params, ...
                ternary(testCase.useAnalytical,'analytical','numerical'));

            I = eye(size(F));
            diffNorm = norm(F*Finv - I,'fro')/norm(I,'fro');
            testCase.verifyLessThanOrEqual(diffNorm,0.1);
            testCase.verifyLessThanOrEqual(comp.err_rel,0.1);
        end

        function testScalingWithN(testCase)
            params1 = struct('N',4,'t',1,'Delta',0.05,'mu',0.02,'g',0.01);
            params2 = struct('N',8,'t',1,'Delta',0.05,'mu',0.02,'g',0.01);

            [F1,~,~] = multiparameter_qfi(params1, ...
                ternary(testCase.useAnalytical,'analytical','numerical'));
            [F2,~,~] = multiparameter_qfi(params2, ...
                ternary(testCase.useAnalytical,'analytical','numerical'));

            % Check quadratic scaling in N
            ratio = F2(1,1)/F1(1,1);
            expectedRatio = (params2.N/params1.N)^2;
            testCase.verifyEqual(ratio,expectedRatio,'AbsTol',1e-10);
        end

        function testInvalidN(testCase)
            params = struct('N',1,'t',1,'Delta',0.05,'mu',0.02,'g',0.01);
            testCase.verifyError(@() multiparameter_qfi(params,'numerical'), ...
                                 'multiparameter_qfi:invalidN');
        end

        function testParameterSweep(testCase)
            for N = [4,6,8]
                params = struct('N',N,'t',1,'Delta',0.05,'mu',0.02,'g',0.01);
                [~,~,comp] = multiparameter_qfi(params, ...
                    ternary(testCase.useAnalytical,'analytical','numerical'));
                testCase.verifyLessThanOrEqual(comp.err_rel,0.2);
            end
        end

        function testNumericalStability(testCase)
            % Stress test with extreme values
            params = struct('N',10,'t',1e-6,'Delta',1e-3,'mu',1e-4,'g',1e-5);
            [~,~,comp] = multiparameter_qfi(params,'numerical'); % force stable mode
            testCase.verifyLessThanOrEqual(comp.err_rel,0.2);
        end
    end
end

% Small helper for inline ternary operator
function out = ternary(cond,a,b)
    if cond, out = a; else, out = b; end
end
