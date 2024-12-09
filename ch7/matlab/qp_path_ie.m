% Program: qp_path_ie.m
% Title: Nonfeasible-initialization primal-dual 
% path-following algorithm for the convex QP problem 
% in Eq. (13.12)
% Description: Implements a primal-dual path-following
% algorithm that applies the method in Sec. 13.4.3 to 
% the QP problem in Eq. (13.12).
% Input:     
%         H: positive semidefinite Hessian matrix
%         A: full row-rank constraint matrix A
%      p, b: input vectors
%        x0: initial point such that A*x0 - b > 0
%      epsi: tolerance for duality gap
% Output:   
%        xs: solution vector
%        fs: value of objective function at xs
%         k: number of iterations at convergence
% Example:
% Find a minimizer of the QP problem
%   minimize 0.5*x'*H*x +x'*p
%   subject to A*x >= b
% where
%   H = [2 -1; -1 1]; p = [-3 -1]'; 
%   A = [-1 -1; -2 1; 1 0; 0 1]; b = [-2 -2 0 0]';
% Solution:
% Execute the following commands:
% H = [2 -1; -1 1]; p = [-3 -1]'; 
% A = [-1 -1; -2 1; 1 0; 0 1]; b = [-2 -2 0 0]';
% x0 = [0.5 1]'; epsi = 1e-6;
% xs= qp_path_ie(H,p,A,b,x0,epsi)
% ======================================
function xs = qp_path_ie(H,p,A,b,x0,epsi)
% disp(' ')
% disp('Program qp_path_ie.m')
% Data initialization.
n = length(x0);
q = length(b);
rho_q = q + 1.5*sqrt(q);
a_max = 1-1e-6;
x = x0(:);
y = A*x - b;
mu = ones(q,1);
ym = y.*mu;
gap = sum(ym);
k = 0;
% Iteration begins.
while gap > epsi,
  tau = gap/rho_q;
  rd = -H*x - p + A'*mu;
  ksi = tau - y.*mu;
  ymi = mu./y;
  yksi = ksi./y;
  YMI = diag(ymi);
  G = inv(H + A'*YMI*A);
  % Calculate d_x,d_lmd, and d_mu.
  d_x = G*(A'*yksi + rd);
  d_y = A*d_x;
  d_mu = (ksi - (mu.*d_y))./y;
  % Calculate step size.
  ind = find(d_y < 0);
  a_p = min(y(ind)./(-d_y(ind)));
  ind = find(d_mu < 0);
  a_d = min(mu(ind)./(-d_mu(ind)));
  a_k = a_max*min([1 a_p a_d]);
  % Form new iterate. 
  x = x + a_k*d_x;
  mu = mu + a_k*d_mu;
  y = A*x - b;
  % Compute duality gap and check convergence.
  ym = y.*mu;
  gap = sum(ym);
  % Update iteration index, etc.
  k = k + 1;
end
xs = x;
% fs = 0.5*xs'*(H*xs + 2*p);