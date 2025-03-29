%% Example 1 of paper "Minimally Conservative Controlled-Invariant Set Synthesis Using Control Barrier Certificates"
% solved by the method proposed in Reference [9]
clear; close all; clc  
%% Adding SeDuMi, Mosek, SOSTOOLS, and Yalmip:
% addpath(genpath('PATH TO SOSTOOLS'))
% addpath(genpath('PATH TO Mosek'))

options.solver='mosek'; 


tic
%% Hyper-Parameters
obstacles = 5;     % <-- either 1, 2 or 5
it_num = 30; 

deg_poly = 20; deg_sos = deg_poly/2;  Deg_h = 4; uM = 1;
eps_clean = 1e-5;  d_tol = 1e-10;
% saving data
gam_hist = zeros(it_num,1);
B_hist = [];

pvar x1 x2; x = [x1;x2];
dpvar d t

%% Define Van Der Pole Occilator Dynamics:
mu = 1;
fx = [x2;            % Van der Pol oscillator
      mu*(1-x1^2)*x2 - x1];
gx = [0;
      1];

% We use LQR cost-to-go as initial guess (calculation can be found in lqr_calc.m)
S = lqr_calc();
% S =[3.3784    0.4142;
%     0.4142    2.6818];

h0 = 0.01 - transpose(x)*S*x;

% --------- defining safe and unsafe sets
sx = (-x1^2-x2^2+4);          % unsafe 1: (-x1^2-x2^2+4) < 0
sx1 = (-x1^2+4);              % unsafe 1: (-x1^2+4) < 0
sx2 = (-x2^2+4);              % unsafe 2: (-x2^2+4) < 0
sx3 = (x1-1)^2+(x2-1)^2-0.04; % unsafe 3: (x1-1)^2+(x2-1)^2-0.04 < 0
sx4 = (x1+1)^2+(x2+1)^2-0.04; % unsafe 4: (x1-1)^2+(x2+1)^2-0.04 < 0
sx5 = (x1+1)^2+(x2-1)^2-0.04; % unsafe 5: (x1+1)^2+(x2-1)^2-0.04 < 0

x_eps = zeros(2,1);  S_eps = eye(2);
h_sol = h0; d_prev = 0;
beta_minus = -0.1; %beta_plus = 1;
kappa = 0.1;



%% While loop in the Algorithm 3 of the paper for producing CBF
% Fix h(i), solve SOS program (29) to find d(i)
% ==================================================
for k=1:it_num
    disp("========================================================")
    fprintf('Iteration = %d \n', k)
    disp("========================================================")
    
    
    % Constraint (29b) - ensure the ellipsoid is contained within the 
    % safe region defined by `h_sol` 
    prog_29 = sosprogram(x);
    prog_29 = sosdecvar(prog_29,d);
    prog_29 = sossetobj(prog_29,-d);
    [prog_29,psi] = sossosvar(prog_29,monomials(x,0:deg_sos));
    consr1 = transpose(x-x_eps)*S_eps*(x-x_eps) - d + psi*h_sol;
    prog_29 = sosineq(prog_29, consr1);

    % solve program 29
    prog_29 = sossolve(prog_29, options);
    d_sol = sosgetsol(prog_29, d)

    clear prog_29
    
    % if d(i) − d(i−1) < tol then converged = True, else
    % Fix h(i)(x), solve SOS program (30) to find μ(i)(x), φ(i)(x)
    % Fix μ(i)(x), φ(i)(x), solve SOS program (31) to find h(i+1), i = i + 1

     if double(d_sol) - d_prev < d_tol   
        ss = fimplicit(matlabFunction(p2s(h_sol)), 'b', 'LineWidth', 1, 'DisplayName',"Method [9]");
        
        disp("========================================================")
        disp("The program is Converged")
        disp("========================================================")

        break;
     end

    % plot inner ellipsoid
    sss = +transpose(x-x_eps)*S_eps*(x-x_eps) - double(d_sol);
    ss = fimplicit(matlabFunction(p2s(sss)), 'm', 'LineWidth', 1);
    ss.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on

    d_prev = double(d_sol);

    % =================================================================
    % Fix h(i)(x), solve SOS program (30) to find μ(i)(x), φ(i)(x)
    prog_30 = sosprogram(x);

    % Calculate Lie derivatives of `h_sol` with respect to the system dynamics
    diffh_sol = [diff(h_sol,x1) diff(h_sol,x2)];
    hdot_sol_u = diffh_sol*(fx + gx*uM);
    hdot_sol_l = diffh_sol*(fx + gx*(-uM));

    % Constraint (27a),(27b) - enforce the CBF condition using the S-procedure
    [prog_30,mu0] = sossosvar(prog_30,monomials(x,0:deg_sos));
    [prog_30,mu1] = sossosvar(prog_30,monomials(x,0:deg_sos));
    [prog_30,mu2] = sossosvar(prog_30,monomials(x,0:deg_sos));

    constr1 = (1+mu0)*(beta_minus-h_sol) + mu1*(hdot_sol_u + kappa*h_sol) + ...
    mu2*(hdot_sol_l + kappa*h_sol);
    prog_30 = sosineq(prog_30, constr1);

    % Constraint (28a),(28b)- ensure the unsafe regions are
    % contained within the sublevel set {x | h(x) <= 0}
    [prog_30,phi10] = sossosvar(prog_30,monomials(x,0:deg_sos));
    [prog_30,phi11] = sossosvar(prog_30,monomials(x,0:deg_sos));
    constr2 = -(1+phi10)*h_sol + phi11*sx1; 
    prog_30 = sosineq(prog_30, constr2);

    % Repeat the process of defining sum-of-squares variables and adding
    % constraints for the remaining unsafe regions sx2,..., sx5.
    [prog_30,phi20] = sossosvar(prog_30,monomials(x,0:deg_sos));
    [prog_30,phi21] = sossosvar(prog_30,monomials(x,0:deg_sos));
    constr3 = -(1+phi20)*h_sol + phi21*sx2;
    prog_30 = sosineq(prog_30, constr3);

    [prog_30,phi30] = sossosvar(prog_30,monomials(x,0:deg_sos));
    [prog_30,phi31] = sossosvar(prog_30,monomials(x,0:deg_sos));
    constr4 = -(1+phi30)*h_sol + phi31*sx3;
    prog_30 = sosineq(prog_30, constr4);
    
    [prog_30,phi40] = sossosvar(prog_30,monomials(x,0:deg_sos));
    [prog_30,phi41] = sossosvar(prog_30,monomials(x,0:deg_sos));
    constr5 = -(1+phi40)*h_sol + phi41*sx4;
    prog_30 = sosineq(prog_30, constr5);
    
    [prog_30,phi50] = sossosvar(prog_30,monomials(x,0:deg_sos));
    [prog_30,phi51] = sossosvar(prog_30,monomials(x,0:deg_sos));
    constr6 = -(1+phi50)*h_sol + phi51*sx5;
    prog_30 = sosineq(prog_30, constr6);

     % solve program 30 [7]
    prog_30 = sossolve(prog_30, options);

    % Extract the solutions for the polynomial variables. [7]
    mu0_sol = sosgetsol(prog_30, mu0);
    mu1_sol = sosgetsol(prog_30, mu1);
    mu2_sol = sosgetsol(prog_30, mu2);
    phi10_sol = sosgetsol(prog_30, phi10);
    phi11_sol = sosgetsol(prog_30, phi11);
    phi20_sol = sosgetsol(prog_30, phi20);
    phi21_sol = sosgetsol(prog_30, phi21);
    phi30_sol = sosgetsol(prog_30, phi30);
    phi31_sol = sosgetsol(prog_30, phi31);
    phi40_sol = sosgetsol(prog_30, phi40);
    phi41_sol = sosgetsol(prog_30, phi41);
    phi50_sol = sosgetsol(prog_30, phi50);
    phi51_sol = sosgetsol(prog_30, phi51);

    clear prog_30

    % =================================================================
    % Fix μ(i)(x), φ(i)(x), solve SOS program (31) to find h(i+1), i = i + 1

    % Define a new SOS program [10]
    prog_31 = sosprogram(x);
    prog_31 = sosdecvar(prog_31,t);
    prog_31 = sossetobj(prog_31,-t);

    [prog_31,h] = sospolyvar(prog_31,monomials(x,0:Deg_h));

    [prog_31,nu] = sossosvar(prog_31,monomials(x,0:deg_sos));

    diffh = [diff(h,x1) diff(h,x2)];
    hdot_u = diffh*(fx + gx*uM);
    hdot_l = diffh*(fx + gx*(-uM));

    % Constraint (31b)
    constr1 = h -t - nu*(double(d_sol) - transpose(x-x_eps)*S_eps*(x-x_eps));
    prog_31 = sosineq(prog_31, constr1);

    % Constraint (27a) 
    constr2 = (1+mu0_sol)*(beta_minus-h) + mu1_sol*(hdot_u + kappa*h) + ...
    mu2_sol*(hdot_l + kappa*h);
    prog_31 = sosineq(prog_31, constr2);


    % Constraints (28a)
    constr3 = -(1+phi10_sol)*h + phi11_sol*sx1;
    prog_31 = sosineq(prog_31, constr3);
    
    constr4 = -(1+phi20_sol)*h + phi21_sol*sx2;
    prog_31 = sosineq(prog_31, constr4);
    
    constr5 = -(1+phi30_sol)*h + phi31_sol*sx3;
    prog_31 = sosineq(prog_31, constr5);
    
    constr6 = -(1+phi40_sol)*h + phi41_sol*sx4;
    prog_31 = sosineq(prog_31, constr6);
    
    constr7 = -(1+phi50_sol)*h + phi51_sol*sx5;
    prog_31 = sosineq(prog_31, constr7);

    % Constraint (31e)
    h_anchor = subs(h,x1,0);
    h_anchor = subs(h_anchor,x2,0);
    constr9 = 1-h_anchor;
    prog_31 = sosineq(prog_31, constr9);

    % Solve program 31
    prog_31 = sossolve(prog_31, options);

    % Extract the solution for h(x) [10]
    h_sol = sosgetsol(prog_31, h);
    h_sol = cleanpoly(h_sol,eps_clean);

    % Plot the updated CBF 
    if k == it_num
        ss = fimplicit(matlabFunction(p2s(h_sol)), 'b', 'LineWidth', 1, 'DisplayName',"Method [9]");
        % ss.Annotation.LegendInformation.IconDisplayStyle = 'off';
    else
        ss = fimplicit(matlabFunction(p2s(h_sol)), 'k', 'LineStyle', '--');
        ss.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end

    hold on
    clear prog_31

end

toc

h_sol
cleanpoly(h_sol,eps_clean)

%% Plotting
fun2 = p2s(h0);
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
       
    % fixed points
    x_e = zeros(n,1);    u_e = 0;

    % System Dynamics
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



