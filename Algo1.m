function [L, K, E, m] = Algo1(sys, m)
%ALGO1 Observer gain design for LPV systems via LMI optimization
%
%   [L, K, E, m] = Algo1(sys, m)
%
%   Designs observer gains L and K using LMIs for an LPV system and checks
%   invertibility of a key matrix E for a given delay horizon m.
%
% Inputs:
%   sys : structure with system matrices (A, B, C) and MM operator
%   m   : maximum horizon to test for invertibility
%
% Outputs:
%   L, K : observer gains (cell arrays of matrices)
%   E    : matrix used in set-membership update
%   m    : selected delay horizon that ensures invertibility

% Initialize alpha and beta (scalars affecting LMI feasibility)
alpha = 0.01;
beta = 0.012;

n  = length(sys.A);         % Number of vertices (polytope)
na = size(sys.A{1}, 1);     % State dimension
nc = size(sys.C, 1);        % Output dimension

% SDP options
options = sdpsettings('solver', 'SDPT3', 'verbose', 0);

Notdone = true;

while Notdone && (alpha < 1) && (beta < 1)
    %% Step 2: Solve LMI for observer gain L (inequality 18)
    const = [];
    for i = 1:n
        P{i} = sdpvar(na, na);
        Z{i} = sdpvar(na, nc);
        X{i} = sdpvar(na, na);
        const = [const; P{i} >= 0];
    end

    for i = 1:n
        for j = 1:n
            Aij = sys.A{i};
            F = [P{j} - X{i} - X{i}', (1/alpha) * (X{i} * Aij - Z{i} * sys.C);
                 (1/alpha) * (X{i} * Aij - Z{i} * sys.C)', -P{i}];
            const = [const; F <= 0];
        end
    end

    diagnostic = optimize(const, [], options);
    if diagnostic.problem ~= 0
        warning('LMI for gain L is infeasible.');
        alpha = alpha + 0.001;
        beta = beta + 0.001;
        continue;
    end

    for i = 1:n
        L{i} = value(X{i}) \ value(Z{i});
    end

    %% Step 3: Solve LMI for observer gain K (inequalities 19 and 24/25)
    const = [];
    const1 = [];
    const2 = [];

    for i = 1:n
        P{i} = sdpvar(na, na);
        S{i} = sdpvar(na, nc);
        Y{i} = sdpvar(na, na);
        const = [const; P{i} >= 0];
    end

    for i = 1:n
        for j = 1:n
            Aij = sys.A{i};
            F = [P{j} - Y{i} - Y{i}', Y{i} * Aij - S{i} * sys.C;
                 (Y{i} * Aij - S{i} * sys.C)', -P{i}];
            const = [const; F <= 0];
        end

        const1 = [const1; ...
            2 * beta * Y{i} - (Y{i} * Aij - S{i} * sys.C) - (Y{i} * Aij - S{i} * sys.C)' <= 0];
        const2 = [const2; ...
           -2 * beta * Y{i} + (Y{i} * Aij - S{i} * sys.C) + (Y{i} * Aij - S{i} * sys.C)' <= 0];
    end

    diagnostic = optimize([const; const1], [], options);
    if diagnostic.problem ~= 0
        diagnostic = optimize([const; const2], [], options);
    end

    if diagnostic.problem ~= 0
        warning('LMI for gain K is infeasible.');
        alpha = alpha + 0.001;
        beta = beta + 0.001;
        continue;
    end

    for i = 1:n
        K{i} = value(Y{i}) \ value(S{i});
    end

    %% Step 4: Check invertibility of E over [1, m]
    for l = 1:m
        E1 = eye(na);
        E2 = eye(na);
        for i = 1:l
            idx = l + 2 - i;
            E1 = E1 * (sys.MM(sys.A, idx) - sys.MM(L, idx) * sys.C);
            E2 = E2 * (sys.MM(sys.A, idx) - sys.MM(K, idx) * sys.C);
        end

        E = inv(E1) - inv(E2);

        % Check if E is numerically invertible and well-conditioned
        if cond(E) < 10^(na + 1) && abs(det(E)) > 1e-6
            Notdone = false;
            m = l;  % Store the feasible m value
            fprintf('Observer gains successfully found (m = %d)\n', m);
            break;
        end
    end

    % Increment parameters if solution is not feasible
    alpha = alpha + 0.001;
    beta = beta + 0.001;
end

end
