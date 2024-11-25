% input:
% [x0,lam0]: initial point
% tau: weight in merit function for constraint violation
% mu: weight for penalty term
% rou: bound on \delta
% epsi: tolerance for convergence
% output:
% xs: local minimizer obtained.
% fs: objective function at xs.
% k: number of iterations performed.
% X: trajectory of the iterates.
% Example: [xs,fs,k,X] = prob15_4([3 2 3]',[1 1]',20,1,0.5,1e-6);
%%
function [xs, fx, k, X] = prob15_4(x0, lam0, eta, tau, rou, epsi)
xk = x0(:);
lam1k = lam0(1); lam2k = lam0(2); k = 0; X = xk;
f = xk(1)^4 + 2*xk(2)^4 + xk(3)^4 - (xk(1)^2)*(xk(2)^2) - (xk(1)^2)*(xk(3)^2);
err = 1;
while err >= epsi
    xk1 = xk(1); xk2 = xk(2); xk3 = xk(3);
    gf1 = 4*xk1^3 - 2*xk1*(xk2^2) - 2*xk1*(xk3^2);
    gf2 = 8*xk2^2 - 2*(xk1^2)*xk2;
    gf3 = 4*xk3^3 - 2*(xk1^2)*xk3;
    gf = [gf1, gf2, gf3]'; 
    ak1 = xk1^4 + xk2^4 + xk3^4 - 25; ak2 = 8*xk1^2 + 14*xk2^2 + 7*xk3^2 - 56;
    ga11 = 4*xk1^3; ga12 = 4*xk2^3; ga13 = 4*xk3^3; ga21 = 16*xk1; ga22 = 28*xk2; ga23 = 14*xk3;
    Ae1 = [ga11, ga12, ga13]; Ae2 = [ga21, ga22, ga23];
    cvx_begin quiet
        variable dt(3,1)
        variable u1(1,1)
        variable v1(1,1)
        variable u2(1,1)
        variable v2(1,1)
        variable r(1,1)
        dual variable y1
        dual variable y2

        minimize(0.5*dt'*Zk*dt + gf'*dt + eta*(u1+v1+u2+v2+r));
        subject to
            y1:ak1 + Ae1*dt == u1 - v1;
            y2:ak2 + Ae2*dt == u2 - v2;
            norm(dt) <= rou + r;
            u1 >= 0; v1 >= 0; u2 >= 0; r >= 0;
    cvx_end
    nuvw = norm([u1, v1, u2, v2, rou]);
    if nuvw <= 1e-6
        dtk = dt;
        lam1_qp = y1;
        lam2_qp = y2;
    else 
        ak1 = ak1 - u1 + v1;
        ak2 = ak2 - u2 + v2;
        rou_t = rou + r;
        cvx_begin quite
            variable dt(3,1)
            dual variable y1
            dual variable y2

            minimize(0.5*dt'*Zk*dt + gf'*dt)';
            subject to
                y1: Ae1*dt == -ak1;
                y2: Ae2*dt == -ak2;
                norm(dt) <= rou_t;
          cvx_end
          dtk = dt;
          lam1_qp = y1;
          lam2_qp = y2;
    end
    dlam1 = lam1_qp - lam1k;
    dlam2 = lam2_qp - lam1k;
    Kt = 40;
    alf = 0:1/(Kt-1):1;
    psi = zeros(Kt, 1);
    for i = 1:Kt
        alfi = alf(i);
        xw = xk + alfi*dtk;
        x1 = xw(1); x2 = xw(2); x3 = xw(3);
        t1 = x1^4 + 2*x2^4 + x3^4 - (x1^2)*(x2^2) - (x1^2)*(x3^2);
        t2 = abs(x1^4 + x2^2 + x3^4 - 25);
        t3 = abs(8*x1^2 + 14*x2^2 + 7*x3^2 - 56);
        psi(i) = t1 + tau*(t2+t3);
    end
    [~, ind] = min(psi); alfk = alf(ind); 
    xnew = xk + alfk*dtk;
    lam1k1 = lam1k + alfk*dlam1;
    lam1k1 = lam2k + alfk*dlam2;
    x1 = xnew(1); x2 = xnew(2); x3 = xnew(3);
    gf1 = 4*x1^3 - 2*x1*(x2^2) - 2*x1*(x3^2);
    gf2 = 8*x2^3 - 2*(x1^2)*x2;
    gf3 = 4*x3^2 - 2*(x1^2)*x3;
    gfnew = [gf1, gf2, gf3]';
    ga11 = 4*x1^3; ga12 = 4*x2^3; ga13 = 4*x3^3; ga21 = 16*x1; ga22 = 28*x2; ga23 = 14*x3;
    Ae1new = [ga11, ga12, ga13]; Ae2new = [ga21, ga22, ga23];
    gamk = (gfnew - gf) + (Ae1new - Ae1)'*lam1k1 + (Ae2new - Ae2)'*lam2k1;
    ct1 = dtk'*gamk; ct2 = dtk'*Zk*dtk;
    if ct1 >= 0.2*ct2
        thet = 1;
    else
        thet = (0.8*ct2)/(ct2 - ct1);
    end
    ht = Zk*dtk;
    ek = thet*gamk + (1-thet)*ht;
    Zk = Zk + (ek*ek')/(dtk'*ek) - (ht*ht')/ct2;
    err = norm([(xnew-xk); (lam1k1-lam1k); (lam2k1-lam2k)]);
    xk = xnew; lam1k = lam1k1; lam2k = lam2k1;
    fk = xk(1)^4 + 2*xk(2)^4 + xk(3)^4 - (xk(1)^2)*(xk(2)^2) - (xk(1)^2)*(xk(3)^2);
    f = [f fk]; X = [X xk];
    k = k + 1;
end
xs = x; fs = fk; 
figure(1)
set(gca, 'fontsize', 13, 'fontname', 'times');
plot(0:1:k, f, 'k-', 'linewidth', 1.4)
xlabel({'Iteration';'(b)'}, 'fontsize',14, 'fontname','times');
ylable('{\itf}{\rm(}{\itx}{\rm)}','fontsize', 14, 'fontname','times');
grid
axis square
% asix([0, k, -1, 3])
disp('satisfaction of 1st equality constraint:');
xs(1)^4 + xs(2)^4 + xs(3)^4 - 25;
disp('satisfaction of 2nd equality constraint:');
8*xs(1)^2 + 14*xs(2)^2 + 7*xs(3)^2 - 56;