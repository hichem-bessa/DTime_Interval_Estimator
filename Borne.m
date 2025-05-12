function [B_bar, B_] = Borne(sys, k, invE, E1, E2, omega_borne)
%BORNE Computes bounding vectors for LPV set observer
%
%   [B_bar, B_] = Borne(sys, k, invE, E1, E2, omega_borne)
%
%   This function computes the bounds that define the size of the interval
%   set in which the system state lies at time step k, using matrices from
%   previous observer estimates and disturbance bounds.
%
% Inputs:
%   sys          : structure containing system model and observer settings
%   k            : current time step
%   invE         : inverse of (E1 - E2)
%   E1, E2       : observer matrices from invproduc
%   omega_borne  : matrix [omega_bar, omega_], i.e., upper and lower bounds of disturbance
%
% Outputs:
%   B_bar        : upper bound vector for state uncertainty
%   B_           : lower bound vector for state uncertainty

n = length(sys.A{1});     % System state dimension
b = inv(E1) - inv(E2);    % Initial difference in inverse matrices

% Accumulate matrix differences over m-1 past steps
for i = 2:sys.m
    [E11, E22] = invproduc(sys, k, i - 1);
    b = b + (E1 \ E11 - E2 \ E22);
end

% Apply inverse of E and identity transformation
b = invE * b * eye(n);

% Compute positive and negative parts
A_plus = max(b, 0);
A_minus = A_plus - b;

% Extract omega bounds
omega_bar = omega_borne(:, 1);  % Upper bound of disturbance
omega_    = omega_borne(:, 2);  % Lower bound of disturbance

% Final bounding terms
B_bar = A_plus * omega_bar - A_minus * omega_;
B_    = A_plus * omega_    - A_minus * omega_bar;

end
