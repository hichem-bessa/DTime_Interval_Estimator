% main_simulation.m
% Simulation of Discrete-Time Interval Estimator for Vehicle Side-Slip Angle Estimation
% Author: Hichem BESSAFA
% Date: 2022
% Reference: https://doi.org/10.1016/j.ifacol.2022.11.294

% Clear workspace and command window
close all;
clear;
clc;

%% System Parameters
% Vehicle physical parameters
sys = initializeSystemParameters();

%% System Dynamics Calculation
sys = calculateSystemDynamics(sys);

%% Polytopic Model Generation
sys = generatePolytopicModel(sys);

%% State Estimation and Gain Calculation
% Apply Algorithm 1 for estimation
[sys.L, sys.K, E, sys.m] = Algo1(sys, 30);

%% Simulation Initial Conditions
sys = setInitialConditions(sys);

%% Simulation Parameters
[t, u] = setupSimulationParameters(sys);

%% Run Simulation
[x, borne] = simsyscomp(sys, t, u);

% Extract simulation results
[x_r, xi, nu, x_hat] = extractSimulationResults(x);

%% Plotting Results
plotResults(t, x_r, x_hat, borne);

% Save results
%saveResults(t, x_r, borne);

%% Helper Functions
function sys = initializeSystemParameters()
    % Vehicle physical parameters
    M = 1529.98;      % Vehicle mass (kg)
    Iz = 4607.47;     % Moment of inertia (kg·m²)
    lf = 1.13906;     % Distance from CG to front axle (m)
    lr = 2.77622 - lf;% Distance from CG to rear axle (m)

    % Cornering stiffness with uncertainty
    Cf = 1.024661677993701e+05;  % Front cornering stiffness
    deltaCf = Cf * 0.05;         % Uncertainty in front cornering stiffness
    Cr = 1.024661677993701e+05;  % Rear cornering stiffness 
    deltaCr = Cr * 0.05;         % Uncertainty in rear cornering stiffness
    
    sys = struct('M', M, 'Iz', Iz, 'lf', lf, 'lr', lr, ...
                 'Cf', Cf, 'deltaCf', deltaCf, 'Cr', Cr, 'deltaCr', deltaCr);
end

function sys = calculateSystemDynamics(sys)
    % Calculate system matrices and parameters
    sys.a_ = (sys.deltaCf + sys.deltaCr) / sys.M;
    a = (sys.Cf + sys.Cr) / sys.M;

    sys.b_ = -sys.deltaCf * sys.lf + sys.deltaCr * sys.lr;
    b = -sys.Cf * sys.lf + sys.Cr * sys.lr;

    sys.c_ = -(sys.deltaCf * sys.lf^2 + sys.deltaCr * sys.lr^2) / sys.Iz;
    c = -(sys.Cf * sys.lf^2 + sys.Cr * sys.lr^2) / sys.Iz;

    % Sampling parameters
    sys.Ts = 0.1;  % Sampling time (s)

    % Uncertainty function
    sys.rho{1} = @(k) 0.06 - 1/20 * abs(sin(0.2 * k));
    sys.rho{2} = @(k) sys.rho{1}(k)^2;

    % Uncertainty bounds
    sys.borne_rho = [0.06 0.01; 0.06^2 0.01^2];
end

function sys = generatePolytopicModel(sys)
    % Compute system matrices with uncertainty
    sys.omega = @(k) [sys.a_ / 1 * abs(cos(k)); 
                      sys.b_ / 1 * abs(sin(k)); 
                      sys.c_ / 1 * abs(sin(k))];

    % Define matrix generation functions
    sys.deltaA = @(k) [...
        -sys.omega(k)' * [1; 0; 0] * sys.rho{1}(k), sys.omega(k)' * [0; 1; 0] / sys.M * sys.rho{2}(k);
        sys.omega(k)' * [0; 1; 0] / sys.Iz, sys.omega(k)' * [0; 0; 1] * sys.rho{1}(k)];

    a = (sys.Cf + sys.Cr) / sys.M;
    b = -sys.Cf * sys.lf + sys.Cr * sys.lr;
    c = -(sys.Cf * sys.lf^2 + sys.Cr * sys.lr^2) / sys.Iz;
    
    sys.AA = @(k) round((...
        [-a * sys.rho{1}(k), b / sys.M * sys.rho{2}(k) - 1; 
         b / sys.Iz, c * sys.rho{1}(k)] + ...
        sys.deltaA(k)) * sys.Ts + eye(2), 14);

    % Initial system matrices
    sys.A{1} = [0, -0.5; b / (2 * sys.Iz), 0] * sys.Ts + 0.5 * eye(2);
    sys.B{1} = [0; sys.Cf * sys.lf / (2 * sys.Iz)] * sys.Ts;

    % Auxiliary matrix generation functions
    A_rho = @(rho) [-a * rho(1), b / sys.M * rho(2) - 1; 
                    b / sys.Iz, c * rho(1)] * sys.Ts + eye(2) - sys.A{1};
    B_rho = @(rho) [sys.Cf / sys.M * rho(1); sys.Cf * sys.lf / sys.Iz] * sys.Ts - sys.B{1};

    % Generate polytopic model
    [mu, A, B] = PolytopicModel(sys, sys.borne_rho, A_rho, B_rho);

    % Update system matrices
    sys.MM = @(M, k) MM_0(M, mu, k);
    sys.A = [sys.A, A];
    sys.B = [sys.B, B];

    % Input matrix
    sys.BB = @(k) [sys.Cf / sys.M * sys.rho{1}(k); sys.Cf * sys.lf / sys.Iz] * sys.Ts;
    sys.C = [0, 1];
end

function sys = setInitialConditions(sys)
    sys.x0 = [0.1; 0.6];      % Initial state
    sys.xi0 = [0.015; 0.2];   % Initial auxiliary state
    sys.nu0 = [0.02; 0.01];   % Initial noise state
    sys.b0 = zeros(2, 2);     % Initial bias

    % Compute uncertainty bounds
    omega1_bar = (sys.a_ * 0.05 * 0.06 + sys.b_ / sys.Iz * 0.45 * 0.06^2) * sys.Ts;
    omega2_bar = (sys.b_ / sys.M * 0.05 - sys.c_ * 0.45 * 0.06) * sys.Ts;
    sys.omega_borne = [omega1_bar, -omega1_bar; 
                       omega2_bar, -omega2_bar];
end

function [t, u] = setupSimulationParameters(sys)
    tf = 6;  % Simulation time (s)
    t = 0 : sys.Ts : tf;
    u = 0.1 * sin(10 * t);  % Input signal
end

function [x_r, xi, nu, x_hat] = extractSimulationResults(x)
    x_r = x(1:2, :);      % Real state
    xi = x(3:4, :);       % Auxiliary state
    nu = x(5:6, :);       % Noise state
    x_hat = x(7:8, :);    % Estimated state
end

function plotResults(t, x_r, x_hat, borne)
    % Plot 1: State Comparison
    figure;
    subplot(2, 1, 1);
    plot(t, x_r(1, :), 'b-^', t, x_hat(1, :), 'r--', 'LineWidth', 2);
    legend('$x_1$', '$\hat{x}_1$', 'Interpreter', 'Latex', 'FontSize', 17, 'Location', 'Best');
    ylabel('Side-Slip angle ($\beta$) (rad)', 'Interpreter', 'Latex', 'FontSize', 17, 'FontWeight', 'bold');

    subplot(2, 1, 2);
    plot(t, x_r(2, :), 'b-^', t, x_hat(2, :), 'r--', 'LineWidth', 2);
    legend('$x_2$', '$\hat{x}_2$', 'Interpreter', 'Latex', 'FontSize', 17, 'Location', 'Best');
    xlabel('Time (s)', 'Interpreter', 'Latex', 'FontSize', 17, 'FontWeight', 'bold');
    ylabel('Yaw rate ($r$) (rad/s)', 'Interpreter', 'Latex', 'FontSize', 17, 'FontWeight', 'bold');

    % Plot 2: Estimation Error
    figure;
    plot(t, abs(x_r(2, :) - x_hat(2, :)), 'LineWidth', 2);
    legend('$|x_2 -\hat{x}_2|$', 'Interpreter', 'Latex');
    title('Estimation Error');

    % Plot 3: Detailed State and Bounds Analysis
    bmax = squeeze(borne(1, :, :));
    bmin = squeeze(borne(2, :, :));

    % Calculate mean absolute errors
    e_11 = abs(x_r(1, :) - bmax(1, :) + x_r(1, :) - bmin(1, :)) / 2;
    e_21 = abs(x_r(2, :) - bmax(2, :) + x_r(2, :) - bmin(2, :)) / 2;

    figure;
    subplot(3, 1, 1);
    plot(t, x_r(1, :), 'k-^', t, x_hat(1, :), t, bmin(1, :), 'r--', t, bmax(1, :), 'b--', 'LineWidth', 2);
    ylabel('Side-Slip angle ($\beta$) (rad)', 'Interpreter', 'Latex', 'FontSize', 17, 'FontWeight', 'bold');
    legend({'$x_1$', '$\hat{x}_1$', '$x^-_1$', '$x^+_1$'}, 'Interpreter', 'Latex', 'FontSize', 12, 'Orientation', 'horizontal');

    subplot(3, 1, 2);
    plot(t, x_r(2, :), 'k-^', t, x_hat(2, :), t, bmin(2, :), 'r--', t, bmax(2, :), 'b--', 'LineWidth', 2);
    ylabel('Yaw rate ($r$) (rad/s)', 'Interpreter', 'Latex', 'FontSize', 17, 'FontWeight', 'bold');
    legend({'$x_2$', '$\hat{x}_2$', '$x^-_2$', '$x^+_2$'}, 'Interpreter', 'Latex', 'FontSize', 12, 'Orientation', 'horizontal');

    subplot(3, 1, 3);
    plot(t, e_11, t, e_21, 'LineWidth', 2);
    ylabel('Mean absolute error', 'Interpreter', 'Latex', 'FontSize', 17, 'FontWeight', 'bold');
    xlabel('Time (s)', 'Interpreter', 'Latex', 'FontSize', 17, 'FontWeight', 'bold');
    legend({'$e_{\beta} $', '$e_r$'}, 'Interpreter', 'Latex', 'FontSize', 12, 'Orientation', 'horizontal');
end

function saveResults(t, x_r, borne)
    bmax = squeeze(borne(1, :, :));
    bmin = squeeze(borne(2, :, :));
    e_11 = abs(x_r(1, :) - bmax(1, :) + x_r(1, :) - bmin(1, :)) / 2;
    e_21 = abs(x_r(2, :) - bmax(2, :) + x_r(2, :) - bmin(2, :)) / 2;
    save('FLPV', 't', 'x_r', 'bmin', 'bmax', 'e_11', 'e_21');
end