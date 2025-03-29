%% Example 2 of paper "Minimally Conservative Controlled-Invariant Set Synthesis Using Control Barrier Certificates"
%% Trajectory simulations
clear; clc
close all;

% Define symbolic variables
syms x1 x2 x3

% Define system dynamics
fx = [x2 - x3^3;
      x3 - x1^2;
      -x1 - 2*x2 - x3 + x2^3]; 
gx = [0, 0;
      1, 0;
      0, 1]; 
uM = 1; % Maximum control input

% Define Control Barrier Function (CBF)
B_sol = -3.4787*x1^4 - 1.4794*x1^3*x2 + 1.1806*x1^3*x3 - 8.5279*x1^2*x2^2 ...
  - 3.2141*x1^2*x2*x3 - 4.7107*x1^2*x3^2 - 6.1398*x1*x2^3 - ...
  0.37263*x1*x2^2*x3 - 5.8863*x1*x2*x3^2 + 1.0879*x1*x3^3 - 3.1362*x2^4 - ...
  0.98845*x2^3*x3 - 2.5898*x2^2*x3^2 + 1.8654*x2*x3^3 - 1.9447*x3^4 - ...
  4.6316*x1^3 - 4.1814*x1^2*x2 + 3.091*x1^2*x3 - 6.5449*x1*x2^2 + ...
  1.3564*x1*x2*x3 - 4.7198*x1*x3^2 - 0.47284*x2^3 - 1.5177*x2^2*x3 - ...
  1.0353*x2*x3^2 - 0.9676*x3^3 + 1.1363*x1^2 - 2.7483*x1*x2 + ...
  3.7891*x1*x3 - 3.0139*x2^2 - 1.0478*x2*x3 - 1.0987*x3^2 + 25.9565*x1 + ...
  9.7687*x2 + 3.2276*x3 + 23.3425;

% Calculate h_dot
[A_cbc, b_cbc, B] = b_dot_calculation(B_sol, fx, gx);

% Initial and goal states
X0 = [1.35, 0.82, -0.38];
color = 'm';
Xg = [0, 0, 0];

% Simulation parameters
tf = 2;                        % Final time [s]
tf_PD = 2.6;                   % Final time for PD controller [s]
fR = 50;                       % Frame rate [fps]
dt = 1/fR;                     % Time resolution [s]
time = linspace(0, tf, tf*fR); % Time vector [s]
size = tf/dt;                  % Number of time steps
size_PD = tf_PD/dt;            % Number of time steps for PD controller

% Controller type
controller = "PD_CBF_QP";

% Initialize variables for saving data
state_hist = zeros(size+1, 3);
state_hist(1, :) = X0;
u_hist = zeros(size, 2);
CBC_hist = zeros(size, 1);
u_pd_hist = zeros(size, 2);
B_hist = zeros(size, 1);
V_hist = zeros(size, 1);
delta_hist = zeros(size, 1);

% Initialize state and parameters
X = X0; alpha = 0; gamma = 1;

% Simulation loop
for i = 1:size
    ccc = i/fR;
    if mod(i, fR) == 0
        disp([num2str(ccc), ' second simulated!'])
    end

    if controller == "PD_CBF_QP"
        % Calculate and save Barrier Function value
        x1 = X(1); x2 = X(2); x3 = X(3);
        BB = double(subs(subs(subs(B, x1), x2), x3));
        B_hist(i, :) = BB;

        % Nominal Controller --> PD Controller
        u_pd = control_pd(X', Xg');
        u_pd_hist(i, :) = u_pd;

        % CBF-QP Controller
        [u, BB] = control2(X, u_pd, A_cbc, b_cbc, uM, BB, alpha);
        X = sim(X, dt, u);
        u_hist(i, :) = u;
        state_hist(i+1, :) = X;
        CBC_hist(i, :) = BB;
    end
end

%% Plot results

% Plot Trajectories
figure
plot3(state_hist(:, 1), state_hist(:, 2), state_hist(:, 3), 'Color', color, 'LineWidth', 1.5)
hold on
plot3(X0(1), X0(2), X0(3), 'o', 'Color', color, 'LineWidth', 1.5)
xlabel("$x_1$", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'latex')
ylabel("$x_2$", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'latex')
zlabel("$x_3$", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'latex')
hold on

% Plot CBF
ss = fimplicit3(B_sol, 'b', 'FaceAlpha', 0.01, 'LineStyle', ':', 'EdgeColor', 'interp', 'MeshDensity', 20);
ss.Annotation.LegendInformation.IconDisplayStyle = 'off';
grid on

ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
ax.ZAxis.LineWidth = 1;

return

%% Function to calculate h_dot
function [A_cbc, b_cbc, B] = b_dot_calculation(B_sol, fx, gx)
    syms x1 x2 x3

    diffB = [diff(B_sol, x1), diff(B_sol, x2), diff(B_sol, x3)];
    A_cbc = -diffB * gx;
    b_cbc = diffB * fx;
    B = B_sol;
end

%% Function to calculate V_dot (not used in this example)
function [A_clf, b_clf, V] = V_dot_calculation(V_sol, fx, gx)
    syms x1 x2 x3

    diffV = [diff(V_sol, x1), diff(V_sol, x2), diff(V_sol, x3)];
    A_clf = diffV * gx;
    b_clf = -diffV * fx;
    V = V_sol;
end

%% PD control function
function u_pd = control_pd(state, Xg)
    uM = evalin('base', 'uM');
    Kx = [3.0645, 4.0364, -0.0413;
          0.2685, -0.0413, 2.3039];
    u_pd = -Kx * state;
end

%% Simulation function
function Y = sim(Y, dt, u)
    t = [0, dt];
    [~, yy] = ode45(@(t, Y) Permisive(t, Y, u), t, Y);
    Y = [yy(end, 1), yy(end, 2), yy(end, 3)];
end

%% System dynamics function
function dy = Permisive(~, y, u)
    y1 = y(1); y2 = y(2); y3 = y(3);
    dy1 = y2 - y3^3;
    dy2 = y3 - y1^2 + u(1);
    dy3 = -y1 - 2*y2 - y3 + y2^3 + u(2);
    dy = [dy1; dy2; dy3];
end

%% CBF-QP control function
function [u, BB] = control2(state, u_pd, A_cbc, b_cbc, uM, BB, alpha)
    x1 = state(1); x2 = state(2); x3 = state(3);

    H = eye(2);
    f = -u_pd;

    A = double(subs(subs(subs(A_cbc, x1), x2), x3));
    b = double(subs(subs(subs(b_cbc, x1), x2), x3)) + alpha * BB;

    options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', 'Display', 'off');
    lb = [-uM; -uM];
    ub = [uM; uM];

    [u, ~, ~, ~, ~] = quadprog(H, f, A, b, [], [], lb, ub, [], options);
end