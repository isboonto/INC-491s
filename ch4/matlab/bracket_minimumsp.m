function [a, c] = bracket_minimumsp(f, x, s, k)
% x: starting point 
% s: step size 
% k: gamma (weight multiplier) 

  switch nargin
      case 4
          k = 2.0;
      case 3
          s = 1e-2; k = 2.0;
      case 2
          x = 0; s = 1e-2; k = 2.0;
  end

  a = x; ya = f(x); b = a + s; yb = f(a + s);

  % change the direction 
  if yb > ya
      [a, b] = deal(b,a);
      [ya, yb] = deal(yb, ya);
      s = -s;
  end
  % yb < ya 
  flag = 1;
  while flag
    s = s * k;
    c = b + s; yc = f(b + s);
    if yc < yb
        a = b; ya = yb; b = c; yb = yc;
    else 
        if c < a
            [a, c] = deal(c,a);
        end
        flag = 0;
    end
  end
end