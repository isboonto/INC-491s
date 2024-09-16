function [a, c] = bracket_minimum(f, x, s, k)
% x: starting point 
% s: step size 
% k: gamma (weight multiplier) 
  switch nargin
      case 3
          k = 2.0;
      case 2
          s = 1e-2; k = 2.0;
      case 1
          x = 0; s = 1e-2; k = 2.0;
  end
  
   
  a = x; ya = f(x); b = a + s; yb = f(a + s);

  % change the direction 
  if yb > ya
      [a, b] = swap(a,b);
      [ya, yb] = swap(ya, yb);
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