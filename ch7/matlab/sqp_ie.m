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
function [xs,fs,k] = sqp_ie(fcname,gname,yname,x0,mu0,epsi)
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