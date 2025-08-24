classdef test_biorthogonal_qfi < matlab.unittest.TestCase
    % TEST_BIORTHOGONAL_QFI Unit tests for biorthogonal_qfi function
    %
    % Verifies single-state and all-state QFI calculations for various systems.

    methods (Test)
        function testSingleStateQFI(testCase)
            % --------------------------
            % Single-state QFI positive check
            % --------------------------
            sys = example_system();       % helper: returns struct with psi_R, psi_L, eigenvals, H_params
            param = 'gamma';
            opts.state = 1;

            F = biorthogonal_qfi(sys, param, opts);

            tolQFI = -1e-10;  % allow tiny negative due to numerical noise
            testCase.verifyGreaterThanOrEqual(F, tolQFI);
        end

        function testAllStatesQFI_Positive(testCase)
            % --------------------------
            % All states QFI positive check
            % --------------------------
            sys = example_system();
            param = 'gamma';
            N = numel(sys.eigenvals);

            tolQFI = -1e-10;

            for n = 1:N
                opts.state = n;
                F = biorthogonal_qfi(sys, param, opts);
                testCase.verifyGreaterThanOrEqual(F, tolQFI);
            end
        end
    end
end

%% ------------------------ Helper function ------------------------
function sys = example_system()
    % EXAMPLE_SYSTEM Returns a minimal 2-site PT-symmetric system for testing
    N = 2;
    t = 1.0;
    gamma = 0.5;
    Delta = 0.1;

    H_params.model_type = 'pt_symmetric';
    H_params.N = N;
    H_params.t = t;
    H_params.gamma = gamma;
    H_params.Delta = Delta;

    H = construct_hamiltonian(H_params);

    [V,D] = eig(full(H));
    psi_R = V;
    psi_L = inv(V)';       % biorthogonal left eigenvectors
    eigenvals = diag(D);

    sys.psi_R = psi_R;
    sys.psi_L = psi_L;
    sys.eigenvals = eigenvals;
    sys.H_params = H_params;
end
