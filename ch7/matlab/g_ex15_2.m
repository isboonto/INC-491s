% Program: g_ex15_2.m
% Description: Evaluates the gradient of the  objective and 
% constraint functions for Example 15.2. It is used to test 
% function sqp_ie sqp_ie, sqp_ie_c, sqp_ie_h, and sqp_ie_p. 
%=======================================
function z = g_ex15_2(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
Ak = [-x1/2+.5 -2*x2 0 0;
0 0 -(5*x3+3*x4)/4+11/2 -(3*x3+5*x4)/4+13/2];
gk = [x1-x3, x2-x4, x3-x1, x4-x2]';
z = [gk Ak'];