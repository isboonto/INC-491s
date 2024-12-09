% Program: y_ex15_2.m
% Description: Evaluates the Hessian of the objective and 
% constraint functions for Example 15.2.  It is used to 
% test function sqp_ie. 
%=======================================
function z = y_han1(xmuk)
    mu = xmuk(3:4);
    mu1 = mu(1);
    mu2 = mu(2);

    z = [2*mu1 + 2*mu2, 0; 0 , 2*mu2];
end