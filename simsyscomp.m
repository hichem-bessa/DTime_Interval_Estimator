function [xk, b] = simsyscomp(sys, t, u)
%SIMSYSCOMP Simulates the LPV observer-based system with finite-time bounding.
%
%   [xk, b] = simsyscomp(sys, t, u)
%
% Inputs:
%   sys : system structure containing dynamics, observers, and parameters
%   t   : time vector
%   u   : control input vector over time
%
% Outputs:
%   xk : state evolution matrix [x; xi; nu; xi_bounded]
%   b  : bounding set over time (upper and lower bounds)

na = length(sys.A{1});      % System order
T = length(t);              % Number of time steps

% Initialize state vectors
xk = zeros(4 * na, T);
xk(1:na, 1) = sys.x0;                         % Real state
xk(na+1:2*na, 1) = sys.xi0;                   % Observer 1
xk(2*na+1:3*na, 1) = sys.nu0;                 % Observer 2
xk(3*na+1:4*na, 1) = sys.xi0;                 % Bounded state

% Initialize bounding box history
b = zeros(2, na, T);
b(:, :, 1) = sys.b0;

for k = 1:T-1
    % Real system evolution
    xk(1:na, k+1) = sys.AA(k-1) * xk(1:na, k) + sys.BB(k-1) * u(k);

    % Observer with gain L
    xk(na+1:2*na, k+1) = ...
        sys.MM(sys.A, k-1) * xk(na+1:2*na, k) + ...
        sys.BB(k-1) * u(k) + ...
        sys.MM(sys.L, k-1) * sys.C * (xk(1:na, k) - xk(na+1:2*na, k));

    % Observer with gain K
    xk(2*na+1:3*na, k+1) = ...
        sys.MM(sys.A, k-1) * xk(2*na+1:3*na, k) + ...
        sys.BB(k-1) * u(k) + ...
        sys.MM(sys.K, k-1) * sys.C * (xk(1:na, k) - xk(2*na+1:3*na, k));

    % Bounding interval computation
    if k < sys.m
        % Before enough steps have passed, use initial bounding
        xk(3*na+1:4*na, k+1) = sys.b0(:, 2);
        b(1, :, k+1) = sys.b0(:, 1);
        b(2, :, k+1) = sys.b0(:, 2);
    else
        % After `m` steps, compute bounding with inverse matrices
        [E1, E2, invE] = invproduc(sys, k, sys.m);
        invE1 = inv(E1);
        invE2 = inv(E2);

        % Finite-time set observer state update
        xk(3*na+1:4*na, k+1) = ...
            invE * (invE1 * xk(na+1:2*na, k+1) ...
                  - invE2 * xk(2*na+1:3*na, k+1) ...
                  - xk(na+1:2*na, k - sys.m + 1) ...
                  + xk(2*na+1:3*na, k - sys.m + 1));

        % Compute bounding box
        [B_bar, B_] = Borne(sys, k, invE, E1, E2, sys.omega_borne);
        b(1, :, k+1) = xk(3*na+1:4*na, k+1) + B_bar;
        b(2, :, k+1) = xk(3*na+1:4*na, k+1) + B_;
    end
end
end
