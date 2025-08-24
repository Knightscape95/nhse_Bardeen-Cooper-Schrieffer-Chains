function H = construct_hamiltonian(params)
    % ===============================================
    % construct_hamiltonian.m
    %
    % Build sparse non-Hermitian BdG Hamiltonian
    % Supports:
    %   - NHSE model (nonreciprocal hopping)
    %   - PT-symmetric model (gain/loss)
    %
    % INPUT:
    %   params.model_type = 'nhse' | 'pt_symmetric'
    %   params.N          = number of sites
    %   params.t          = hopping amplitude
    %   params.gamma      = asymmetry or gain/loss strength
    %   params.Delta      = pairing amplitude
    %
    % OUTPUT:
    %   H  : sparse BdG Hamiltonian (2N x 2N)
    % ===============================================

    arguments
        params struct
    end

    N     = params.N;
    t     = params.t;
    gamma = params.gamma;
    Delta = params.Delta;

    % Hamiltonian dimension (BdG doubling)
    dim = 2 * N;
    H = sparse(dim, dim);

    % ------------------------
    % Hopping terms
    % ------------------------
    for i = 1:N-1
        switch params.model_type
            case 'nhse'
                % Nonreciprocal hopping: t ± γ/2
                H(i, i+1)   = t + gamma/2;   % right hopping
                H(i+1, i)   = t - gamma/2;   % left hopping

            case 'pt_symmetric'
                % Symmetric hopping
                H(i, i+1)   = t;
                H(i+1, i)   = t;

            otherwise
                error('Unknown model_type: %s', params.model_type);
        end
    end

    % ------------------------
    % On-site gain/loss
    % ------------------------
    if strcmp(params.model_type, 'pt_symmetric')
        for i = 1:N
            if mod(i,2)==0
                H(i,i) = 1i*gamma;   % gain
            else
                H(i,i) = -1i*gamma;  % loss
            end
        end
    end

    % ------------------------
    % Pairing terms (BdG structure)
    % ------------------------
    for i = 1:N
        H(i, N+i)     = Delta;   % particle-hole coupling
        H(N+i, i)     = conj(Delta);
    end
end
