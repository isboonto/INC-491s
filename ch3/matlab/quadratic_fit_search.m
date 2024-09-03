function [a, b, c] = quadratic_fit_search(f, a, b, c, n)
    ya = f(a); yb = f(b); yc = f(c);
    for i = 1:n-3
        x = 0.5*(ya*(b^2 - c^2) + yb*(c^2 - a^2) + yc*(a^2 - b^2)) ...
            / (ya*(b - c) + yb*(c - a) + yc*(a - b));
        yx = f(x);
        if x > b
            if yx > yb
                c = x; yc = yx;
            else
                a = b; ya = yb; b = x; yb = yx;
            end
        elseif x < b
            if yx > yb
                a = x; ya = yx;
            else
                c = b; yc = yb; b = x; yb = yx;
            end
        end
    end
end
