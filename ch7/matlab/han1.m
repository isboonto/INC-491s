% Program: f_ex15_2.m
% Description: Evaluates the objective and constraint functions 
% for Example 15.2.  It is used to test functions sqp_ie,  
% sqp_ie_c, sqp_ie_h, and sqp_ie_p. 
%=======================================
function z = han1(x)
    x1 = x(1);
    x2 = x(2);

    ck = [x1^2-x2; 
        x1^2+x2^2-1];
    fk = x1 + x2;
    z = [fk; ck];
end
