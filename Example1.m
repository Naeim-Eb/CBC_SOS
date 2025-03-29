%% Example 1 of paper "Minimally Conservative Controlled-Invariant Set Synthesis Using Control Barrier Certificates"
clear; close all; clc  
%% Adding SeDuMi, Mosek, SOSTOOLS, and Yalmip:
% addpath(genpath('PATH TO SOSTOOLS'))
% addpath(genpath('PATH TO Mosek'))

options.solver='mosek'; 

tic
%% Hyper-Parameters
obstacles = 5;     % <-- either 1, 2 or 5
it_num = 50; 

deg_poly = 4; deg_sos = deg_poly/2;  deg_control = 3; Deg_B = 4; uM = 1;
eps1 = 0.001;   eps2 = 0.00;   eps3 = 0.00; eps_clean = 1e-6;  gamma_treshold = 1e-6;
% saving data
gam_hist = zeros(it_num,1);
B_hist = [];

pvar x1 x2;
dpvar gam %gam2 gam3
x = [x1;x2];

%% Define Van Der Pole Occilator Dynamics:
mu = 1;
fx = [x2;            % Van der Pol oscillator
      mu*(1-x1^2)*x2 - x1];
gx = [0;
      1];

% We use LQR cost-to-go as initial guess (calculation can be found in lqr_calc.m)
S = lqr_calc();


b0 = 0.01 - transpose(x)*S*x;

% b0 = 0.01 - transpose(x)*x; % --> this also gives us an answer but not
% always


% --------- defining safe and unsafe sets
sx = (-x1^2-x2^2+4);          % unsafe 1: (-x1^2-x2^2+4) < 0
sx1 = (-x1^2+4);              % unsafe 1: (-x1^2+4) < 0
sx2 = (-x2^2+4);              % unsafe 2: (-x2^2+4) < 0
sx3 = (x1-1)^2+(x2-1)^2-0.04; % unsafe 3: (x1-1)^2+(x2-1)^2-0.04 < 0
sx4 = (x1+1)^2+(x2+1)^2-0.04; % unsafe 4: (x1-1)^2+(x2+1)^2-0.04 < 0
sx5 = (x1+1)^2+(x2-1)^2-0.04; % unsafe 5: (x1+1)^2+(x2-1)^2-0.04 < 0


%% Initialize the SOSP:
% Initialization: Given b0(x), find u(x) and λ1(x), λ2(x) by solving (6)
% ==================================================
prog_Initialization = sosprogram(x);

% The poly multipliers
[prog_Initialization,lamb1] = sospolyvar(prog_Initialization,monomials(x,0:deg_poly));
[prog_Initialization,lamb2] = sospolyvar(prog_Initialization,monomials(x,0:deg_poly));
[prog_Initialization,lamb3] = sospolyvar(prog_Initialization,monomials(x,0:deg_poly));

% The SOS multipliers
[prog_Initialization,sigma1] = sossosvar(prog_Initialization,monomials(x,0:deg_sos));

% control poly
[prog_Initialization,ux] = sospolyvar(prog_Initialization,monomials(x,0:deg_control));

% Bdot
diffB = [diff(b0,x1) diff(b0,x2)];
Bdot_sol = diffB*(fx + gx*ux);

% =============================
% Inequality constraints
% Constraint 1 ( (6d) in Paper ): Bdot(x) > 0 on boundry of B
prog_Initialization = sosineq(prog_Initialization, Bdot_sol - lamb1*b0 -eps2);  % <-- CBC

% Constraint 2 ( (6c) in Paper ): feasible u on boundry of B
% -uM<u<uM
a_u=[+1; -1];
c_u = uM*[+1; +1];
prog_Initialization = sosineq(prog_Initialization, a_u(1)*ux+c_u(1) - lamb2*b0-eps3); 
prog_Initialization = sosineq(prog_Initialization, a_u(2)*ux+c_u(2) - lamb3*b0-eps3);

% Constraint 3 ( (6b) in Paper ): B in B
% prog_Initialization = sosineq(prog_Initialization, b0 - sigma1*(b0-0.01));
prog_Initialization = sosineq(prog_Initialization, b0 - sigma1*b0);

% Constraint 4 ( (6a) in Paper ) : B<0 in S<0 (B<0 in unsafe set)
if obstacles == 1
    [prog_Initialization,sigma2] = sossosvar(prog_Initialization,monomials(x,0:deg_sos));
    prog_Initialization = sosineq(prog_Initialization, -b0 + sigma2*sx -eps1);  
end

if obstacles == 2
    [prog_Initialization,sigma2] = sossosvar(prog_Initialization,monomials(x,0:deg_sos));
    [prog_Initialization,sigma3] = sossosvar(prog_Initialization,monomials(x,0:deg_sos));

    prog_Initialization = sosineq(prog_Initialization, -b0 + sigma2*sx1 -eps1); 
    prog_Initialization = sosineq(prog_Initialization, -b0 + sigma3*sx2 -eps1); 
end

if obstacles == 5
    [prog_Initialization,sigma2] = sossosvar(prog_Initialization,monomials(x,0:deg_sos));
    [prog_Initialization,sigma3] = sossosvar(prog_Initialization,monomials(x,0:deg_sos));
    [prog_Initialization,sigma4] = sossosvar(prog_Initialization,monomials(x,0:deg_sos));     
    [prog_Initialization,sigma5] = sossosvar(prog_Initialization,monomials(x,0:deg_sos));     
    [prog_Initialization,sigma6] = sossosvar(prog_Initialization,monomials(x,0:deg_sos));     

    prog_Initialization = sosineq(prog_Initialization, -b0 + sigma2*sx1 -eps1); 
    prog_Initialization = sosineq(prog_Initialization, -b0 + sigma3*sx2 -eps1); 
    prog_Initialization = sosineq(prog_Initialization, -b0 + sigma4*sx3-eps1); 
    prog_Initialization = sosineq(prog_Initialization, -b0 + sigma5*sx4-eps1);
    prog_Initialization = sosineq(prog_Initialization, -b0 + sigma6*sx5-eps1);
end

% Solve the problem
prog_Initialization = sossolve(prog_Initialization, options);
% Get the resulting invariant set
ux_sol = sosgetsol(prog_Initialization,ux);
lamb1_sol = sosgetsol(prog_Initialization,lamb1);
lamb2_sol = sosgetsol(prog_Initialization,lamb2);
lamb3_sol = sosgetsol(prog_Initialization,lamb3);

B_sol = b0;

clear prog_Initialization

%% While loop in the Algorithm 1 of the paper:
% Initialization: Given u(x), λ1(x), λ2(x) , find b(x) and γ by solving (8)
% ==================================================
for k=1:it_num
    disp("========================================================")
    fprintf('Iteration = %d \n', k)
    disp("========================================================")

    % ==================================================
    %% Fix u(x), λ1(x), λ2(x), solve for B and SOS multiplyers
    % Initialize the SOS
    prog_B = sosprogram(x);

    % The barrier function
    % we can choose polyvar for B, as follows:
    [prog_B,B] = sospolyvar(prog_B,monomials(x,0:Deg_B));
    diffB = [diff(B,x1) diff(B,x2)];
    Bdot = diffB*(fx + gx*ux_sol);

    % The SOS multipliers
    [prog_B,sigma1] = sossosvar(prog_B,monomials(x,0:deg_sos));
    
    % =============================
    % Inequality constraints
    % Constraint 1: Bdot(x) > 0 on boundry of B
    prog_B = sosineq(prog_B, Bdot - lamb1_sol*B -eps2);  % <-- CBC

    % Constraint 2: feasible u on boundry of B
    % -uM<u<uM
    prog_B = sosineq(prog_B, a_u(1)*ux_sol+c_u(1) - lamb2_sol*B-eps3); 
    prog_B = sosineq(prog_B, a_u(2)*ux_sol+c_u(2) - lamb3_sol*B-eps3); 

    % Constraint 3: B in B_
    prog_B = sosineq(prog_B, B - sigma1*B_sol ); 
    
    % Constraint 4: making B larger
    [prog_B,lamb4] = sospolyvar(prog_B,monomials(x,0:deg_poly));
    prog_B = sosdecvar(prog_B,gam);
    prog_B = sosineq(prog_B, B - lamb4*B_sol -gam);
    prog_B = sosineq(prog_B, gam);
    prog_B= sossetobj(prog_B,-gam);


    % Constraint 5: B<0 in S<0 (B<0 in unsafe set)
    if obstacles == 1
        [prog_B,sigma2] = sossosvar(prog_B,monomials(x,0:deg_sos));
        prog_B = sosineq(prog_B, -B + sigma2*sx -eps1);  
    end

    if obstacles == 2
        [prog_B,sigma2] = sossosvar(prog_B,monomials(x,0:deg_sos));
        [prog_B,sigma3] = sossosvar(prog_B,monomials(x,0:deg_sos));
        prog_B = sosineq(prog_B, -B + sigma2*sx1 -eps1); 
        prog_B = sosineq(prog_B, -B + sigma3*sx2 -eps1);  
    end

    if obstacles == 5
        [prog_B,sigma2] = sossosvar(prog_B,monomials(x,0:deg_sos));
        [prog_B,sigma3] = sossosvar(prog_B,monomials(x,0:deg_sos));
        [prog_B,sigma4] = sossosvar(prog_B,monomials(x,0:deg_sos));
        [prog_B,sigma5] = sossosvar(prog_B,monomials(x,0:deg_sos));
        [prog_B,sigma6] = sossosvar(prog_B,monomials(x,0:deg_sos));
        prog_B = sosineq(prog_B, -B + sigma2*sx1 -eps1); 
        prog_B = sosineq(prog_B, -B + sigma3*sx2 -eps1); 
        prog_B = sosineq(prog_B, -B + sigma4*sx3-eps1); 
        prog_B = sosineq(prog_B, -B + sigma5*sx4-eps1);
        prog_B = sosineq(prog_B, -B + sigma6*sx5-eps1);
    end
    
    
    % Solve the problem
    prog_B = sossolve(prog_B, options);

    % Save the temporary variable
    gamma = double(sosgetsol(prog_B,gam));
    if gamma <= gamma_treshold
        fun = p2s(B_sol_prev);
        myfun = matlabFunction(fun);
        ss = fimplicit(myfun, 'b', 'LineWidth', 1, 'DisplayName', "Our CBC");
        ss.Annotation.LegendInformation.IconDisplayStyle = 'off';

        disp("gamma <= gamma_treshold")
        break
    end
    gam_hist(k,1)=gamma;

    % Get the solutions
    B_sol = sosgetsol(prog_B,B);
    cleanpoly(B_sol,eps_clean);

    if p2s(B_sol) == 0
        temp = "K = ";
        disp("B_sol can not found!")
        display(temp + k)
        B_sol = B_sol_prev;
        break
    end
    B_sol_prev = B_sol;
    B_hist = [B_hist; p2s(B_sol)];
    clear prog_B;


    % ==================================================
    %% Fix B, solve for the control input, and poly multipliers (for the next step)
    prog_u = sosprogram(x);

    [prog_u,lamb1] = sospolyvar(prog_u,monomials(x,0:deg_poly));
    [prog_u,lamb2] = sospolyvar(prog_u,monomials(x,0:deg_poly));
    [prog_u,lamb3] = sospolyvar(prog_u,monomials(x,0:deg_poly));

    [prog_u,ux] = sospolyvar(prog_u,monomials(x,0:deg_control));
    
    diffB = [diff(B_sol,x1) diff(B_sol,x2)];
    Bdot_sol = diffB*(fx + gx*ux);

    % =============================
    % Inequality constraints
    % Constraint 1: Bdot(x) > 0 on boundry og B
    prog_u = sosineq(prog_u, Bdot_sol - lamb1*B_sol -eps2);

    % Constraint 2: feasible u on boundry of B
    % -uM<u<uM
    prog_u = sosineq(prog_u, a_u(1)*ux+c_u(1) - lamb2*B_sol-eps3); 
    prog_u = sosineq(prog_u, a_u(2)*ux+c_u(2) - lamb3*B_sol-eps3); 

    % Solve the problem
    prog_u = sossolve(prog_u, options);
    
    % Get the solutions
    ux_sol = sosgetsol(prog_u,ux);
    cleanpoly(ux_sol,eps_clean);
    lamb1_sol = sosgetsol(prog_u,lamb1);
    lamb2_sol = sosgetsol(prog_u,lamb2);
    lamb3_sol = sosgetsol(prog_u,lamb3);

    if p2s(lamb2_sol) == 0
        temp = "K = ";
        disp("u(x) can not found!")
        display(temp + k)
        break
    end
    
    % Delete the temporary variable
    clear prog_u;

    fun = p2s(B_sol);
    myfun = matlabFunction(fun);
    if k == it_num
        ss = fimplicit(myfun, 'b', 'LineWidth', 1, 'DisplayName', "Our CBC");

    else
        ss = fimplicit(myfun, 'k', 'LineStyle', '--');
        ss.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    hold on


end
toc

B_sol
cleanpoly(B_sol,eps_clean)

fun2 = p2s(b0);
myfun2 = matlabFunction(fun2);
ss = fimplicit(myfun2, 'g', 'LineWidth', 1, 'DisplayName', "b_0(x)");
ss.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
% hold on
if obstacles == 1
    sx = evalin('base', 'sx');
    fun = p2s(sx);
    myfun = matlabFunction(fun);
    
    plot_S = fimplicit(myfun, 'r', 'LineWidth', 1, 'DisplayName','State Constraints');
    hold on

end

if obstacles == 2
    sx1 = evalin('base', 'sx1');
    sx2 = evalin('base', 'sx2');
    fun = p2s(sx1*sx2);
    myfun = matlabFunction(fun);
    
    plot_S = fimplicit(myfun, 'r', 'LineWidth', 1, 'DisplayName','State Constraints');
    hold on
    % plot_I = fimplicit(myfun,'g', 'LineWidth', 1, 'DisplayName','Initial invariant set');

    text( 0-1.5 , 3 , "Unsafe", "color", "r", FontSize=15)
    text( 0-0.3 , -3 , "Unsafe", "color", "r", FontSize=15)
    text( 2.5 , 0-0.3 , "Unsafe", "color", "r", "Rotation", 90, FontSize=15)
    text( -2.5 , 0-0.3 , "Unsafe", "color", "r", "Rotation", 90, FontSize=15)
end

if obstacles == 5
    sx1 = evalin('base', 'sx1');
    sx2 = evalin('base', 'sx2');
    sx3 = evalin('base', 'sx3');
    sx4 = evalin('base', 'sx4');
    sx5 = evalin('base', 'sx5');
    fun = p2s(sx1*sx2);
    myfun = matlabFunction(fun);
    
    plot_S = fimplicit(myfun, 'r', 'LineWidth', 1, 'DisplayName','State Constraints');
    hold on

    fun = p2s(sx3);
    myfun = matlabFunction(fun);
    ss1 = fimplicit(myfun,'r', 'LineWidth', 1);
    ss1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on

    fun = p2s(sx4);
    myfun = matlabFunction(fun);
    ss2 = fimplicit(myfun,'r', 'LineWidth', 1);
    ss2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on

    fun = p2s(sx5);
    myfun = matlabFunction(fun);
    ss3 = fimplicit(myfun,'r', 'LineWidth', 1);
    ss3.Annotation.LegendInformation.IconDisplayStyle = 'off';

    text( 0-0.3 , 2.2 , "Unsafe", "color", "r", FontSize=10)
    text( 0-0.3 , -2.2 , "Unsafe", "color", "r", FontSize=10)
    text( 2.2 , 0-0.3 , "Unsafe", "color", "r", "Rotation", 90, FontSize=10)
    text( -2.2 , 0-0.3 , "Unsafe", "color", "r", "Rotation", 90, FontSize=10)
    text( -1-0.18 , 1 , "Unsafe", "color", "r", FontSize=5)
    text( 1-0.18 , 1 , "Unsafe", "color", "r", FontSize=5)
    text( -1-0.18 , -1 , "Unsafe", "color", "r", FontSize=5)

end

axis(2.4*[-1 1 -1 1])
xlabel("$x_1$",'FontWeight', 'bold', 'FontSize', 14,'Interpreter','latex')
ylabel("$x_2$",'FontWeight', 'bold', 'FontSize', 14, 'Interpreter','latex')
legend




%% Calculation of LQR cost-to-go
function S = lqr_calc()

    n = 2;
    syms x [n 1]
    syms u

    x_e = zeros(n,1);    
    u_e = 0;

    mu = 1;
    fx = [x2;            % Van der Pol oscillator
          mu*(1-x1^2)*x2 - x1];
    gx = [0;
          1];

    f = fx+gx*u;

    A = vpa(jacobian(f,x),3);
     
    A = subs(A, x, x_e);
    A = vpa(subs(A, u, u_e), 3);
    
    B = jacobian(f,u);
    B = subs(B, x, x_e);
    B = vpa(subs(B, u, u_e),3);
    
        
    Co = ctrb(A, B);
    rank_Co = rank(Co);

    if rank_Co < length(A)
        error('The system is not controllable!');
    end

    Q = eye(n);     
    R =1;
    [K,S,P] = lqr(double(A),double(B),Q,R);

end


















