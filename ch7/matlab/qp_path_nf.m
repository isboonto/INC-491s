% Program: qp_path_nf.m
% Title: Nonfeasible-initialization primal-dual 
% path-following algorithm for convex QP problems
% Description: Implements Algorithm 13.3 for convex QP 
% problems.
% Theory: See Practical Optimization Sec. 13.4.3.
% Input:     
%             H -- positive semidefinite Hessian matrix
%             A -- full row-rank constraint matrix A
%          p, b -- input vectors
% (x0,lmd0,mu0) -- initial point with x0 > 0 and mu0 > 0
%           rho -- parameter involved in determining 
%                  the value of tau, rho is required 
%                  to be no less than the square root of n
%          epsi -- tolerance for duality gap
% Output:   
%            xs -- solution vector
%            fs -- value of objective function at xs
%             k -- number of iterations at convergence
% Example:
% Find a minimizer of the QP problem
% minimize 0.5*x'*H*x +x'*p
% subject to A*x = b, x >= 0
% where
% H = [4 0 0; 0 1 -1; 0 -1 1]
% p = [-8 -6 -6]'
% A = [1 1 1]
% b = 3
% using the initial value
% x0 = [1 2 2]'
% Solution:
% Execute the following commands:
% H = [4 0 0; 0 1 -1; 0 -1 1]
% p = [-8 -6 -6]'
% A = [1 1 1]
% b = 3
% x0 = [1 2 2]'
% lmd0 = -1
% mu0 = [0.2 0.2 0.2]'
% rho = 3
% epsi = 1e-6
% [xs,fs,k]= qp_path_nf(H,p,A,b,x0,lmd0,mu0,rho,epsi)
% ===================================================
function [xs,fs,k]= qp_path_nf(H,p,A,b,x0,lmd0,mu0,rho,epsi)
disp(' ')
disp('Program qp_path_nf.m')
% Data initialization
n = length(x0);
rho_n = n + rho;
a_max = 1-1e-6;
x = x0(:);
lmd = lmd0(:);
mu = mu0(:);
xm = x.*mu;
gap = sum(xm);
k = 0;
% Iteration begins.
while gap > epsi,
  % Compute Y0, yd, etc.
  tau = gap/rho_n;
  r_p = b - A*x;
  r_d = H*x + p - A'*lmd-mu;
  M = diag(mu);
  X = diag(x);
  Gam = inv(M+X*H);
  GaX = Gam*X*A';
  Y0 = inv(A*GaX);
  yd = Gam*(x.*(mu+r_d) - tau);
  % Calculate d_x,d_lmd, and d_mu.
  d_lmd = Y0*(A*yd + r_p);
  d_x = GaX*d_lmd - yd;
  d_mu = H*d_x - A'*d_lmd + r_d;
  % Calculate step size.
  ind = find(d_x < 0);
  a_p = min(x(ind)./(-d_x(ind)));
  ind = find(d_mu < 0);
  a_d = min(mu(ind)./(-d_mu(ind)));
  a_k = a_max*min([a_p a_d]);
  % Form new iterate. 
  x = x + a_k*d_x;
  lmd = lmd + a_k*d_lmd;
  mu = mu + a_k*d_mu;
  % Compute duality gap and check convergence.
  xm = x.*mu;
  gap = sum(xm);
  % Update iteration index, etc.
  k = k + 1;
  xsa(:,k) = x;
end
xs = xsa;
fs = 0.5*xs(:,end)'*(H*xs(:,end) + 2*p);