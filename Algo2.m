function [L, K, E, m] = Algo2(sys, m)
%ALGO2 Observer gain design for LPV systems via pole placement
%
%   [L, K, E, m] = Algo2(sys, m)
%
%   Computes observer gains L and K using eigenvalue assignment and checks
%   the invertibility of matrix E over a delay horizon m.
%
% Inputs:
%   sys : structure with system matrices (A, B, C) and MM operator
%   m   : maximum delay horizon to check
%
% Outputs:
%   L, K : observer gains (cell arrays)
%   E    : matrix used for set-membership update
%   m    : chosen delay satisfying invertibility condition

n  = length(sys.A);         % Number of LPV vertices
na = size(sys.A{1}, 1);     % State dimension

% Step 1: Initialize eigenvalue bounds
delta_L    = 0.2;
delta_L0K  = 0.3;
delta_K    = 0.9;

delta_star = delta_K - 2 * delta_L0K;

lambda_Lmax = delta_L - 0.1 * delta_L;  % Initial max eigenvalue for L
lambda_Kmin = delta_L0K;                % Initial min eigenvalue for K

Notdone = true;

% Step 2: Loop over eigenvalue ranges and test invertibility
while Notdone && lambda_Lmax > 0 && lambda_Kmin < 1
    Lambda_L = cell(1, n);
    Lambda_K = cell(1, n);
    
    % Step 3: Assign eigenvalues and compute gains using `place`
    for i = 1:n
        Lambda_L{i} = lambda_Lmax - 0.01 * (0:na-1);  % Stable L poles
        Lambda_K{i} = lambda_Kmin + 0.01 * (0:na-1);  % Stable K poles

        % Observer gain matrices via pole placement
        L{i} = place(sys.A{i}', sys.C', Lambda_L{i})';
        K{i} = place(sys.A{i}', sys.C', Lambda_K{i})';
    end

    % Step 4: Check invertibility of E for 1 ≤ l ≤ m
    for l = 1:m
        E1 = eye(na);
        E2 = eye(na);
        
        for i = 1:l
            idx = l + 2 - i;
            E1 = E1 * (sys.MM(sys.A, idx) - sys.MM(L, idx) * sys.C);
            E2 = E2 * (sys.MM(sys.A, idx) - sys.MM(K, idx) * sys.C);
        end
        
        E = inv(E1) - inv(E2);
        
        % Check E's condition number and determinant
        if cond(E) < 10^na && abs(det(E)) > 1e-6
            Notdone = false;
            m = l;
            fprintf('DONE — Gains found with m = %d\n', m);
            break;
        end
    end

    % Slightly adjust eigenvalues for next iteration
    lambda_Lmax = lambda_Lmax - 0.0001;
    lambda_Kmin = lambda_Kmin + 0.0001;
end
end
