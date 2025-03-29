%% Example 2 of paper "Minimally Conservative Controlled-Invariant Set Synthesis Using Control Barrier Certificates"
clear; close all; clc

%% Adding SeDuMi, Mosek, SOSTOOLS, and Yalmip:
% addpath(genpath('PATH TO SOSTOOLS'))
% addpath(genpath('PATH TO Mosek'))

options.solver='mosek'; 

%% Hyper-Parameters
% --------- degree of multipliers and polynomials -------------
it_num = 100;      % Number of iterations for enlarging the invariant set 
deg_poly = 4;      % Degree of the polynomial multipliers 
deg_sos = deg_poly/2; % Degree of the SOS multipliers 
deg_control = 5;    % Degree of the control polynomial 
Deg_B = 4;        % Degree of the barrier function polynomial 
uM = 1;           % Upper bound on the control input (input constraint) 
eps1 = 0.001;      % Small positive constant used in SOSP constraints 
eps2 = 0.00;      % Small positive constant used in SOSP constraints 
eps3 = 0.00;       % Small positive constant used in SOSP constraints 
eps_clean = 1e-10;  % Tolerance for cleaning polynomial coefficients 
gamma_treshold = 1e-10; % Threshold for stopping the iterative enlargement 

% saving data
gam_hist = zeros(it_num,1); % History of gamma values 
B_hist = [];        % History of barrier function polynomials 

% --------- Constructing the vector field (dynamical system) -------------
pvar x1 x2 x3;     % Define polynomial variables x1, x2, x3 
dpvar gam          % Define decision variable gamma 
x = [x1; x2; x3];   % State vector 

% Define the dynamics of the system (Example 2 from ) 
fx = [x2-x3^3;
      x3-x1^2;
      -x1-2*x2-x3+x2^3]; 
gx = [0, 0;
      1, 0;
      0, 1]; 

% --------- defining safe sets and initial set -------------
sx = 81- transpose(x)*x;  % Define a safe set (not used in this example) 
sx1 = (x1-2)^2+(x2-1)^2+(x3-2)^2-1; % Define unsafe set 1 
sx2 = (x1+1)^2+(x2+2)^2+(x3+1)^2-1; % Define unsafe set 2 
sx3 = (x1-0)^2+(x2-0)^2+(x3-6)^2-9; % Define unsafe set 3 
sx4 = (x1-0)^2+(x2-0)^2+(x3+5)^2-9; % Define unsafe set 4 

% Initial barrier function (defines the initial safe set) 
S =[2.1246    1.7522    0.3407;
    1.7522    3.3188    0.8334;
    0.3407    0.8334    1.2238];
B_sol = 0.01 - transpose(x)*S*x;

% Plot the unsafe sets 
fimplicit3(matlabFunction(p2s(sx1)),'r','EdgeColor','w','LineStyle',':','FaceAlpha',0.3,'MeshDensity',30, 'DisplayName', 'Unsafe Zone');
hold on
aa = fimplicit3(matlabFunction(p2s(sx2)),'r','EdgeColor','w','LineStyle',':','FaceAlpha',0.3,'MeshDensity',30);
aa.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
bb = fimplicit3(matlabFunction(p2s(sx3)),'r','EdgeColor','w','LineStyle',':','FaceAlpha',0.3);
bb.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
cc = fimplicit3(matlabFunction(p2s(sx4)),'r','EdgeColor','w','LineStyle',':','FaceAlpha',0.3);
cc.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

%% Initialize the SOSP:
% Initialization: Given b0(x), find u(x) and λ1(x), λ2(x) by solving (6)
% ==================================================

prog = sosprogram(x); % Initialize a sum-of-squares program 

% Define polynomial multipliers lambda1 to lambda5 
[prog,lamb1] = sospolyvar(prog,monomials(x,0:deg_poly));
[prog,lamb2] = sospolyvar(prog,monomials(x,0:deg_poly));
[prog,lamb3] = sospolyvar(prog,monomials(x,0:deg_poly));
[prog,lamb4] = sospolyvar(prog,monomials(x,0:deg_poly));
[prog,lamb5] = sospolyvar(prog,monomials(x,0:deg_poly));

% Define SOS multipliers sigma0 to sigma4 
[prog,sigma0] = sossosvar(prog,monomials(x,0:deg_sos));
[prog,sigma1] = sossosvar(prog,monomials(x,0:deg_sos));
[prog,sigma2] = sossosvar(prog,monomials(x,0:deg_sos));
[prog,sigma3] = sossosvar(prog,monomials(x,0:deg_sos));
[prog,sigma4] = sossosvar(prog,monomials(x,0:deg_sos));

% Define control input polynomials ux1 and ux2 
[prog,ux1] = sospolyvar(prog,monomials(x,0:deg_control));
[prog,ux2] = sospolyvar(prog,monomials(x,0:deg_control));

% Calculate the time derivative of the barrier function 
diffB = jacobian(B_sol); 
Bdot_sol = diffB*(fx + gx*[ux1;ux2]);

% =============================
% Inequality constraints
% =============================

% Constraint 1: Bdot(x) > 0 on boundary of B 
constr1 = (Bdot_sol - lamb1*B_sol -eps2);
prog = sosineq(prog, constr1); 

% Constraint 2: feasible u on boundary of B (input constraints) 
a_u1=[+1; -1]; % Coefficients for input constraint 1
c_u1 = uM*[+1; +1]; % Constant terms for input constraint 1
prog = sosineq(prog, a_u1(1)*ux1+c_u1(1) - lamb2*B_sol-eps3);
prog = sosineq(prog, a_u1(2)*ux1+c_u1(2) - lamb3*B_sol-eps3);

a_u2=[+1; -1]; % Coefficients for input constraint 2
c_u2 = uM*[+1; +1]; % Constant terms for input constraint 2
prog = sosineq(prog, a_u2(1)*ux2 +c_u2(1) - lamb4*B_sol-eps3);
prog = sosineq(prog, a_u2(2)*ux2 +c_u2(2) - lamb5*B_sol-eps3);

% Constraint 3: B<0 in S<0 (B<0 in unsafe set) 
prog = sosineq(prog, -B_sol + sigma0*sx -eps1);
prog = sosineq(prog, -B_sol + sigma0*sx1 -eps1);
prog = sosineq(prog, -B_sol + sigma0*sx2 -eps1);
prog = sosineq(prog, -B_sol + sigma0*sx3 -eps1);
prog = sosineq(prog, -B_sol + sigma0*sx4 -eps1);

% Solve the SOSP 
prog = sossolve(prog, options);

% Get the solutions 
ux1_sol = sosgetsol(prog,ux1); 
ux2_sol = sosgetsol(prog,ux2);
lamb1_sol = sosgetsol(prog,lamb1);
lamb2_sol = sosgetsol(prog,lamb2);
lamb3_sol = sosgetsol(prog,lamb3);
lamb4_sol = sosgetsol(prog,lamb4);
lamb5_sol = sosgetsol(prog,lamb5);
sigma0_sol = sosgetsol(prog,sigma0);
sigma1_sol = sosgetsol(prog,sigma1);
sigma2_sol = sosgetsol(prog,sigma2);
sigma3_sol = sosgetsol(prog,sigma3);
sigma4_sol = sosgetsol(prog,sigma4);

clear prog % Clear the SOSP variable 

%% While loop in the Algorithm 1 of the paper:
% Initialization: Given u(x), λ1(x), λ2(x) , find b(x) and γ by solving (8)
% ==================================================

for k=1:it_num % Iterate to enlarge the invariant set 

    disp("========================================================")
    fprintf('Iteration = %d \n', k) 
    disp("========================================================")

    % ==================================================
    % Fix the poly multipliers and control input, solve for B and SOS multiplyers
    % ==================================================

    prog_B = sosprogram(x); % Initialize a new SOSP for enlarging B 

    % Define the barrier function as a polynomial variable 
    [prog_B,B] = sospolyvar(prog_B,monomials(x,0:Deg_B)); 

    % Calculate the time derivative of B 
    diffB = jacobian(B); 
    Bdot = diffB*(fx + gx*[ux1_sol;ux2_sol]);

    % Define SOS multipliers 
    [prog_B,sigma00] = sossosvar(prog_B,monomials(x,0:deg_sos)); 
    [prog_B,sigma0] = sossosvar(prog_B,monomials(x,0:deg_sos)); 
    [prog_B,sigma1] = sossosvar(prog_B,monomials(x,0:deg_sos));
    [prog_B,sigma2] = sossosvar(prog_B,monomials(x,0:deg_sos));
    [prog_B,sigma3] = sossosvar(prog_B,monomials(x,0:deg_sos)); 
    [prog_B,sigma4] = sossosvar(prog_B,monomials(x,0:deg_sos)); 

    % =============================
    % Inequality constraints
    % =============================

    % Constraint 1: Bdot(x) > 0 on boundary of B 
    prog_B = sosineq(prog_B, Bdot - lamb1_sol*B -eps2); 

    % Constraint 2: feasible u on boundary of B 
    prog_B = sosineq(prog_B, a_u1(1)*ux1_sol+c_u1(1) - lamb2_sol*B-eps3); 
    prog_B = sosineq(prog_B, a_u1(2)*ux1_sol+c_u1(2) - lamb3_sol*B-eps3); 
    prog_B = sosineq(prog_B, a_u2(1)*ux2_sol +c_u2(1) - lamb4_sol*B-eps3);
    prog_B = sosineq(prog_B, a_u2(2)*ux2_sol +c_u2(2) - lamb5_sol*B-eps3);

    % Constraint 3: B in B_prev [14]
    prog_B = sosineq(prog_B, B - sigma00*B_sol );

    % Constraint 4: making B larger [14]
    [prog_B,lamb4] = sospolyvar(prog_B,monomials(x,0:deg_poly)); 
    prog_B = sosdecvar(prog_B,gam); % Define gamma as a decision variable 
    prog_B = sosineq(prog_B, B - lamb4*B_sol -gam); 
    prog_B = sosineq(prog_B, gam);
    prog_B= sossetobj(prog_B,-gam); % Maximize gamma to enlarge B [14]

    % Constraint 5: B<0 in S<0 (B<0 in unsafe set) 
    prog_B = sosineq(prog_B, -B + sigma0*sx -eps1);
    prog_B = sosineq(prog_B, -B + sigma1*sx1 -eps1);
    prog_B = sosineq(prog_B, -B + sigma2*sx2 -eps1);
    prog_B = sosineq(prog_B, -B + sigma3*sx3 -eps1);
    prog_B = sosineq(prog_B, -B + sigma4*sx4 -eps1);

    % Solve the problem 
    prog_B = sossolve(prog_B, options);

    % Get the value of gamma 
    gamma = double(sosgetsol(prog_B,gam)); 

    % Check if gamma is below the threshold 
    if gamma <=gamma_treshold 
        fun = p2s(B_sol_prev);
        myfun = matlabFunction(fun);
        ss = fimplicit3(myfun,'b', 'FaceAlpha',0.15,'LineStyle',':','EdgeColor','interp','MeshDensity',30, 'DisplayName', "Our CBC"); 
        disp("gamma <= gamma_treshold")
        break 
    end

    gam_hist(k,1)=gamma; % Store the gamma value 

    % Get the solution for the barrier function 
    B_sol = sosgetsol(prog_B,B);

    % Check if a valid solution for B was found [15]
    if p2s(B_sol) == 0 
        temp = "K = ";
        disp("B_sol can not found!")
        display(temp + k)
        B_sol = B_sol_prev; % Use the previous solution if no valid solution is found
        break 
    end

    B_sol_prev = B_sol; % Store the current solution for the next iteration
    B_hist = [B_hist; p2s(B_sol)]; % Store the solution in the history
    clear prog_B; % Clear the SOSP variable 

    % ==================================================
    % Fix B, solve for the control input, and poly multipliers
    % ==================================================
    prog_u = sosprogram(x); % Initialize a new SOSP for refining u and lambdas 

    % Define polynomial multipliers 
    [prog_u,lamb1] = sospolyvar(prog_u,monomials(x,0:deg_poly)); 
    [prog_u,lamb2] = sospolyvar(prog_u,monomials(x,0:deg_poly)); 
    [prog_u,lamb3] = sospolyvar(prog_u,monomials(x,0:deg_poly));
    [prog_u,lamb4] = sospolyvar(prog_u,monomials(x,0:deg_poly)); 
    [prog_u,lamb5] = sospolyvar(prog_u,monomials(x,0:deg_poly));

    % Define control input polynomials 
    [prog_u,ux1] = sospolyvar(prog_u,monomials(x,0:deg_control)); 
    [prog_u,ux2] = sospolyvar(prog_u,monomials(x,0:deg_control)); 

    % Calculate the time derivative of B 
    diffB = jacobian(B_sol);
    Bdot_sol = diffB*(fx + gx*[ux1; ux2]);

    % =============================
    % Inequality constraints
    % =============================

    % Constraint 1: Bdot(x) > 0 on boundary of B 
    prog_u = sosineq(prog_u, Bdot_sol - lamb1*B_sol -eps2); 

    % Constraint 2: feasible u on boundary of B 
    prog_u = sosineq(prog_u, a_u1(1)*ux1+c_u1(1) - lamb2*B_sol-eps3);
    prog_u = sosineq(prog_u, a_u1(2)*ux1+c_u1(2) - lamb3*B_sol-eps3);
    prog_u = sosineq(prog_u, a_u2(1)*ux2+c_u2(1) - lamb4*B_sol-eps3);
    prog_u = sosineq(prog_u, a_u2(2)*ux2+c_u2(2) - lamb5*B_sol-eps3);

    % Solve the SOSP 
    prog_u = sossolve(prog_u, options);

    % Get the solutions 
    ux1_sol = sosgetsol(prog_u,ux1);
    ux2_sol = sosgetsol(prog_u,ux2);
    lamb1_sol = sosgetsol(prog_u,lamb1);
    lamb2_sol = sosgetsol(prog_u,lamb2);
    lamb3_sol = sosgetsol(prog_u,lamb3);
    lamb4_sol = sosgetsol(prog_u,lamb4);
    lamb5_sol = sosgetsol(prog_u,lamb5);

    clear prog_u; % Clear the SOSP variable

    % Plot the final barrier function 
    fun = p2s(B_sol); 
    myfun = matlabFunction(fun);
    if k == it_num
        ss = fimplicit3(myfun,'b', 'FaceAlpha',0.15,'LineStyle',':','EdgeColor','interp','MeshDensity',30, 'DisplayName', "Our CBC");
    end
    hold on;
end

% Add labels to the plot [16]
xlabel("$x_1$",'FontWeight', 'bold', 'FontSize', 14,'Interpreter','latex')
ylabel("$x_2$",'FontWeight', 'bold', 'FontSize', 14, 'Interpreter','latex')
zlabel("$x_3$",'FontWeight', 'bold', 'FontSize', 14, 'Interpreter','latex')















