function alpha = line_search(f, x, d)
    p1 = x(1); p2 = x(2); p3 = d(1); p4 = d(2);
    objective = @(alpha) f(p1 + alpha*p3, p2 + alpha*p4);
    [a, b] = bracket_minimum(objective);

    alpha = fminbnd(objective, a, b);
end