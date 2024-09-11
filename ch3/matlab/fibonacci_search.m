% Fibonacci_search
function [a, b] = fibonacci_search(f, a, b, n, epsilon)
    switch nargin 
        case  4
        epsilon = 0.01;
    end
    golden = 1.61803;
    s = (1-sqrt(5))/(1+sqrt(5));
    rho = 1 / (golden * (1-s^(n+1))/(1-s^n));
    d = rho*b + (1-rho)*a;
    yd = f(d);

    for i = 1:n-1
        if i == n-1
            c = epsilon*a + (1-epsilon)*d;
        else
            c = rho*a + (1-rho)*b;
        end
        yc = f(c);
        if yc < yd
            b = d; d = c; yd = yc;
        else
            a = b; b = c;
        end
        rho = 1/ (golden*(1-s^(n-i+1))/(1-s^(n-i)));
    end

    if a > b
        [a, b] = swap(a,b); % [a, b] = deal(b,a)
    end
end