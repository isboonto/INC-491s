function [am, bm] = bisection_search(f, a, b, epsilon)
    if nargin < 4
        epsilon = 1e-8;
    end

    if a > b
        [a, b] = swap(a,b);     % ensure a < b
    end

    ya = f(a); yb = f(b);
    if ya == 0; b = a; end
    if yb == 0; a = b; end
    i = 0;
    while b - a > epsilon
        i = i + 1;
        x = (a + b) /2;
        y = f(x);
        if y == 0 
            a = x; b = x;
        elseif sign(y) == sign(ya)
            a = x;
        else
            b = x;
        end
        am(i,:) = a;
        bm(i,:) = b;
    end
end
