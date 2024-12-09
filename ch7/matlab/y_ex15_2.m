% Program: y_ex15_2.m
% Description: Evaluates the Hessian of the objective and 
% constraint functions for Example 15.2.  It is used to 
% test function sqp_ie. 
%=======================================
function z = y_ex15_2(xmuk)
mu = xmuk(5:6);
mu1 = mu(1);
mu2 = mu(2);
y11 = 1+mu1/2;
y22 = 1+2*mu1;
y33 = 1+5*mu2/4;
y34 = 3*mu2/4;
z = [y11 0 -1 0; 
0 y22 0 -1; 
-1 0 y33 y34; 
0 -1 y34 y33];