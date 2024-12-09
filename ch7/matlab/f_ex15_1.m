% Program: f_ex15_1.m
% Description: Evaluates the objective and constraint functions 
% for Example 15.1. It is used to test function sqp_e. 
%=======================================
function z = f_ex15_1(x)
x1 = x(1);
x2 = x(2);
x3 = x(3); 
fk = -x1^4-2*x2^4-x3^4-x1^2*x2^2-x1^2*x3^2;
ak = [x1^4+x2^4+x3^4-25; 8*x1^2+14*x2^2+7*x3^2-56];
z = [fk; ak];