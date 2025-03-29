%% Example 1 of paper "Minimally Conservative Controlled-Invariant Set Synthesis Using Control Barrier Certificates"
%% Trajectory simulations
clear; close all; clc

% Define symbolic variables
syms x1 x2

% Define unsafe regions
sx1 = (-x1^2 + 4);              % Unsafe 1: (-x1^2 + 4) < 0
sx2 = (-x2^2 + 4);              % Unsafe 2: (-x2^2 + 4) < 0
sx3 = (x1 - 1)^2 + (x2 - 1)^2 - 0.04; % Unsafe 3: (x1 - 1)^2 + (x2 - 1)^2 - 0.04 < 0
sx4 = (x1 + 1)^2 + (x2 + 1)^2 - 0.04; % Unsafe 4: (x1 + 1)^2 + (x2 + 1)^2 - 0.04 < 0
sx5 = (x1 + 1)^2 + (x2 - 1)^2 - 0.04; % Unsafe 5: (x1 + 1)^2 + (x2 - 1)^2 - 0.04 < 0

% Define Control Barrier Function (CBF)
hx = -95.7091*x1^4 - 105.2701*x1^3*x2 - 653.8866*x1^2*x2^2 + ...
    175.8905*x1*x2^3 - 229.5341*x2^4 - 20.1123*x1^3 - 59.5746*x1^2*x2 - ...
    52.3664*x1*x2^2 - 59.8847*x2^3 + 258.8242*x1^2 - 48.5611*x1*x2 - ...
    127.3491*x2^2 + 68.2461*x1 + 71.5848*x2 + 419.7535;

% Calculate h_dot
[A_cbc, b_cbc] = h_dot_calculation(hx);
uM = 1; % Maximum control input

% Initial and goal states
X0 = [-0.992, -0.75];
color = 'r';
Xg = [0, 0];

% Controller type
controller = "PD_CBF_QP";

% Simulation parameters
tf = 2;                        % Final time [s]
tf_PD = 2.6;                   % Final time for PD controller [s]
fR = 50;                       % Frame rate [fps]
dt = 1/fR;                     % Time resolution [s]
time = linspace(0, tf, tf*fR); % Time vector [s]
size = tf/dt;                  % Number of time steps
size_PD = tf_PD/dt;            % Number of time steps for PD controller

% Initialize variables for saving data
state_hist = zeros(size+1, 2);
state_hist(1, :) = X0;
u_hist = zeros(size, 1);
CBC_hist = zeros(size, 1);
u_pd_hist = zeros(size, 1);
B_hist = zeros(size, 1);

% Only PD Controller
u_pd_hist_only_PD = zeros(size_PD, 1);
state_hist_only_PD = zeros(size_PD+1, 2);
state_hist_only_PD(1, :) = X0;
[u_pd_hist_only_PD, state_hist_only_PD] = ...
    only_PD(X0, Xg, size_PD, dt, u_pd_hist_only_PD, state_hist_only_PD, uM);

% Initialize state and alpha
X = X0; alpha = 1;

% Simulation loop
for i = 1:size
    ccc = i/fR;
    if mod(i, fR) == 0
        disp([num2str(ccc), ' second simulated!'])
    end

    % Calculate and save Barrier Function value
    x1 = X(1); x2 = X(2);
    BB = double(subs(subs(hx, x1), x2));
    B_hist(i, :) = BB;

    % Nominal Controller --> PD Controller
    u_pd = control_pd(X, Xg, uM);
    u_pd_hist(i, :) = u_pd;
    u_pd = 0;

    if controller == "PD"
        X = sim(X, dt, u_pd);
        state_hist(i+1, :) = X;
    elseif controller == "PD_CBF_QP"
        [u, BB] = control2(X, u_pd, A_cbc, b_cbc, uM, BB, alpha);
        X = sim(X, dt, u);
        u_hist(i, :) = u;
        state_hist(i+1, :) = X;
        CBC_hist(i, :) = BB;
    end
end

%% Plot results

% Plot Trajectories
figure(1)
ii = plot(X0(1), X0(2), 'o', 'Color', color, 'LineWidth', 1.5);
ii.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
ss = plot(state_hist(:, 1), state_hist(:, 2), 'LineWidth', 1.5, 'LineStyle', '--', 'Color', color);
ss.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

% Plot unsafe regions and CBF
fimplicit(sx1*sx2, 'r', 'LineWidth', 1);
hold on
fimplicit(sx3, 'r', 'LineWidth', 1);
hold on
fimplicit(sx4, 'r', 'LineWidth', 1);
hold on
fimplicit(sx5, 'r', 'LineWidth', 1);
hold on
fimplicit(hx, 'b', 'LineWidth', 1);

axis(2.5*[-1 1 -1 1])

return

%% Function to calculate h_dot
function [A_cbc, b_cbc] = h_dot_calculation(hx)
    syms x1 x2
    mu = 1;
    fx = [x2;            % Van der Pol oscillator dynamics
          mu*(1-x1^2)*x2 - x1];
    gx = [0;
          1];

    diffh = [diff(hx, x1), diff(hx, x2)];
    A_cbc = -diffh * gx;
    b_cbc = diffh * fx;
end

%% Function for PD control only
function [u_pd_hist_only_PD, state_hist_only_PD] = ...
    only_PD(Y, Yg, size, dt, u_pd_hist_only_PD, state_hist_only_PD, uM)
    for i = 1:size
        u_pd = control_pd(Y, Yg, uM);
        u_pd_hist_only_PD(i, :) = u_pd;
        Y = sim(Y, dt, u_pd);
        state_hist_only_PD(i+1, :) = Y;
    end
end

%% PD control function
function u_pd = control_pd(state, Xg, uM)
    K = [0.4142, 2.6818];  % LQR gain for Van der Pol Oscillator System
    u_pd = -K * (state - Xg)';
end

%% Simulation function
function Y = sim(Y, dt, u)
    t = [0, dt];
    [~, yy] = ode45(@(t, Y) VanDerPol(t, Y, u), t, Y);
    Y = [yy(end, 1), yy(end, 2)];
end

%% Van der Pol dynamics
function dy = VanDerPol(~, y, u)
    y1 = y(1); y2 = y(2);
    dy1 = y2;
    dy2 = (1 - y1^2) * y2 - y1 + u;
    dy = [dy1; dy2];
end

%% CBF-QP control function
function [u, BB] = control2(state, u_pd, A_cbc, b_cbc, uM, BB, alpha)
    x1 = state(1); x2 = state(2);
    H = 1;
    f = -u_pd;

    A = double(subs(subs(A_cbc, x1), x2));
    A = [A; 1; -1];

    b = double(subs(subs(b_cbc, x1), x2));
    b = b + alpha * BB;
    b = [b; uM; uM];

    options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', 'Display', 'off');
    [u, ~, ~, ~, ~] = quadprog(H, f, A, b, [], [], [], [], [], options);
end