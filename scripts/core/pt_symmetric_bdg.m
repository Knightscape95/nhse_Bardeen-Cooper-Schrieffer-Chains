function [H_pt, eigvals, psi_R, psi_L] = pt_symmetric_bdg(params)
% Constructs PT-symmetric BdG Hamiltonian and computes biorthogonal eigenpairs

N     = params.N;
t     = params.t;
Delta = params.Delta;
mu    = params.mu;
g     = params.g;

% Single-particle Hamiltonian h (N×N)
stag = diag(((-1).^[1:N])');       % +1,-1,+1,...
h = -mu*eye(N) + 1i*g*stag;

% Nearest-neighbor hopping
for j = 1:N-1
    h(j,j+1) = h(j,j+1) + t;
    h(j+1,j) = h(j+1,j) + t;
end

% Pairing matrix Delta_mat
Delta_mat = zeros(N);
for j = 1:N-1
    Delta_mat(j,j+1) = Delta;
    Delta_mat(j+1,j) = -Delta;  % antisymmetric for BdG
end

% Full BdG Hamiltonian
H_pt = [ h,          Delta_mat;
        -Delta_mat', -h.' ];

% Right eigenvectors
[psi_R, D] = eig(H_pt);
eigvals    = diag(D);

% Left eigenvectors: eigenvectors of H'
[psi_L_raw, D2] = eig(H_pt');

% Match left eigenvectors to right eigenvalues
[~, idx] = ismembertol(conj(eigvals), diag(D2), 1e-12);
psi_L = psi_L_raw(:, idx);

% ----------------- Biorthonormalization -----------------
% Rescale columns so that diag(psi_L' * psi_R) = 1
overlaps = diag(psi_L' * psi_R);
if any(abs(overlaps) < 1e-12)
    error('Biorthogonal overlap too small');
end
S = diag(1 ./ sqrt(overlaps));
psi_R = psi_R * S;
psi_L = psi_L * conj(S);

% Check
err = norm(psi_L' * psi_R - eye(2*N), 'fro');
if err > 1e-10
    warning('Biorthonormalization error = %.2e', err);
end
end
