function dx = obser(t, x, obs, u)
%OBSER Implements an interval observer with H-infinity design for LPV systems.
%
% dx = obser(t, x, obs, u)
%
% Inputs:
%   t   - Current time
%   x   - State vector [x_r; x_hat^+; x_hat^-] where:
%         x_r        : True system state (2x1)
%         x_hatplus  : Upper interval observer state estimate (2x1)
%         x_hatminus : Lower interval observer state estimate (2x1)
%   obs - Structure containing system dynamics, observer parameters
%   u   - Input function handle u(t)
%
% Output:
%   dx  - Derivative of the full observer state [dx_r; dx_hatplus; dx_hatminus]

% Extract components from state vector
x_r        = x(1:2);
x_hatplus  = x(3:4);
x_hatminus = x(5:6);

% Bounding system matrices and inputs
[A_plus, A_barplus]         = bounding(obs.A_plus(t));
[A_minus, A_barminus]       = bounding(obs.A_minus(t));
[x_plus, x_barplus]         = bounding(x_hatplus);
[x_minus, x_barminus]       = bounding(x_hatminus);
[B_plus, B_barplus]         = bounding(obs.B_plus(t));
[B_minus, B_barminus]       = bounding(obs.B_minus(t));
[u_plus, u_barplus]         = bounding(u(t));
[u_minus, u_barminus]       = bounding(u(t));

% True system dynamics
dx_r = obs.AA(t) * x_r + obs.B0(t) * u(t);

% Interval observer dynamics (upper and lower bounds)
dx_hatplus = ...
    (obs.A0(t) - obs.L(t) * obs.C) * x_hatplus ...
    + (A_barplus * x_barplus - A_plus * x_barminus ...
       - A_barminus * x_plus + A_minus * x_minus) ...
    + (B_barplus * u_barplus - B_plus * u_barminus ...
       - B_barminus * u_plus + B_minus * u_minus) ...
    + obs.L(t) * obs.C * x_r;

dx_hatminus = ...
    (obs.A0(t) - obs.L(t) * obs.C) * x_hatminus ...
    + (A_plus * x_plus - A_barplus * x_minus ...
       - A_minus * x_barplus + A_barminus * x_barminus) ...
    + (B_plus * u_plus - B_barplus * u_minus ...
       - B_minus * u_barplus + B_barminus * u_barminus) ...
    + obs.L(t) * obs.C * x_r;

% Output full state derivative
dx = [dx_r; dx_hatplus; dx_hatminus];

end
