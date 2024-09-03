function [a, b] = golden_section_search(f, a, b, n)
    golden = 1.61803;
    rho = golden - 1;
    d = rho * b + (1 - rho)*a;
    yd = f(d);

    for i = 1:n-1
        c = rho*a + (1-rho)*b;
        yc = f(c);
        if yc < yd
            b = d; d = c; yd = yc;
        else
            a = b; b = c;
        end
    end
    if a > b
        [a, b] = swap(a,b);
    end
end
