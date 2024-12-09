x0 = [1 0.5 ]'
mu0 = [1 1]'
epsi = 1e-5
[xs, fs, k] = sqp_ie2(['han1', 'g_han1','y_han1', x0 ,mu0, epsi])
%%
% Program: f_ex15_2.m
% Description: Evaluates the objective and constraint functions 
% for Example 15.2.  It is used to test functions sqp_ie,  
% sqp_ie_c, sqp_ie_h, and sqp_ie_p. 
%=======================================
function z = han1(x)
    x1 = x(1);
    x2 = x(2);

    ck = [x1^2 - x2; x1^2 + x2^2 -1];
    fk = -x1 - x2;
    z = [fk; ck];
end
%%
% Program: y_ex15_2.m
% Description: Evaluates the Hessian of the objective and 
% constraint functions for Example 15.2.  It is used to 
% test function sqp_ie. 
%=======================================
function z = y_han1(xmuk)
    mu = xmuk(5:6);
    mu1 = mu(1);
    mu2 = mu(2);

    z = [2*mu1 + 2*mu2, 0; 0 , 2*mu2];
end
%%
% Program: g_ex15_2.m
% Description: Evaluates the gradient of the  objective and 
% constraint functions for Example 15.2. It is used to test 
% function sqp_ie sqp_ie, sqp_ie_c, sqp_ie_h, and sqp_ie_p. 
%=======================================
function z = g_han1(x)
    x1 = x(1);
    x2 = x(2);

    Ak = [2*x1, -1; 2*x1, 2*x2]; 
    gk = [-1, -1]';
    z = [gk Ak'];
end
%%
% Program: sqp_ie.m
% Title: SQP algorithm for nonlinear problems with 
% inequality constraints.
% Description: Implements the SQP algorithm for nonlinear problems 
% with inequality constraints (Algorithm 15.2). 
% Theory: See Practical Optimization Sec. 15.2.2.
% Input:
%   fcname: function that evaluates the objective and constraint 
%           functions
%    gname: function that evaluates the gradient of the objective 
%           and constraint functions
%    yname: function that evaluates the Hessian of the objective 
%           and constraint functions
%       x0: initial point
%      mu0: initial Lagrange multiplier
%     epsi: termination tolerance
% Output:
%   xs: solution point
%   fs: objective function evaluated at xs.
%   k: number of iterations at convergence
% Example:
% Apply Algorithm 15.2 to solve the minimization problem 
% in Example 15.2.
% Solution:
% Execute the following commands:
% x0 = [1 0.5 2 3]'
% mu0 = [1 1]'
% epsi = 1e-5
% [xs,fs,k] = sqp_ie('f_ex15_2','g_ex15_2','y_ex15_2',x0,mu0,epsi)
% ========================================================
function [xs,fs,k] = sqp_ie2(fcname,gname,yname,x0,mu0,epsi)
disp(' ')
disp('Program sqp_ie.m')
xk = x0(:);
muk = mu0(:);
n = length(x0);
q = length(muk);
q1 = q + 1;
k = 0;
In = eye(n);
d = 1;
while d >= epsi,
    fck = feval(fcname,xk);
    Gk = feval(gname,xk);
    ck = fck(2:q1);
    gk = Gk(:,1);
    Ak = Gk(:,2:q1)';
    xmuk = [xk; muk];
    Yk = feval(yname,xmuk);
    emin = min(eig(Yk));
    if emin <= 0,
        Yk = Yk + (1e-6 - 1.01*emin)*In;
    end
    % d_x = quadprog(Yk,gk,-Ak,ck);
    d_x = qp_path_ie(Yk,gk,Ak,-ck,zeros(n,1),epsi);
    xk = xk + d_x;
    muk = pinv(Ak')*(Yk*d_x + gk);
    k = k + 1;
    d = norm(d_x);
end
format long
disp('Solution point:')
xs = xk
disp('Objective function at the solution point:')
fck = feval(fcname,xs);
fs = fck(1)
format short
disp('Number of iterations at convergence:')
k
end