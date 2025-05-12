function [E1, E2, E] = invproduc(sys, k, m)
%INVPRODUC Compute product of observer matrices and composite inverse
%
%   [E1, E2, E] = invproduc(sys, k, m)
%
%   Computes the matrix products E1 and E2 for observer trajectories of length m,
%   and the composite matrix E used for set-membership estimation.
%
% Inputs:
%   sys : structure containing system matrices and observers
%   k   : current time step
%   m   : horizon length (number of past steps)
%
% Outputs:
%   E1  : product of dynamics with observer gain L
%   E2  : product of dynamics with observer gain K
%   E   : inverse of the difference between inv(E1) and inv(E2)

% Initialize products as identity (scalar 1 used for compatibility)
E1 = 1;
E2 = 1;

% Build the matrix products over the horizon m
for j = 1:m
    A_j  = sys.MM(sys.A, k - j);      % Time-varying system matrix
    L_j  = sys.MM(sys.L, k - j);      % Time-varying observer gain L
    K_j  = sys.MM(sys.K, k - j);      % Time-varying observer gain K
    C    = sys.C;

    E1 = E1 * (A_j - L_j * C);
    E2 = E2 * (A_j - K_j * C);
end

% Final composite matrix for correction
E = inv(inv(E1) - inv(E2));

end
