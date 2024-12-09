% Program: f_ex15_2.m
% Description: Evaluates the objective and constraint functions 
% for Example 15.2.  It is used to test functions sqp_ie,  
% sqp_ie_c, sqp_ie_h, and sqp_ie_p. 
%=======================================
function z = f_ex15_2(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
ck = [-(x1^2/4+x2^2)+x1/2+0.75; 
-(5*x3^2+6*x3*x4+5*x4^2)/8+(11*x3+13*x4)/2-35/2];
fk = 0.5*((x1-x3)^2+(x2-x4)^2);
z = [fk; ck];