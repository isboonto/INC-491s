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