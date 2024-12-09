% Program: qp_ie0.m
% Title: Primal active-set algorithm for convex
% QP problems with inequality constraints.
% Description: Implements Algorithm 13.1 for
% convex QP problems with inequality constraints.
% Theory: See Practical Optimization Sec. 13.3.1.
% Input:
% H -- positive semidefinite Hessian matrix
% A -- full row-rank constraint matrix A
% p, b -- input vectors
% x0 -- a feasible initial point
% Output:
% xs -- solution vector
% fs -- value of objective function at xs
% k -- number of iteration at convergence
% Example:
% Find a minimizer of the QP problem
% minimize 0.5*x'*H*x +x'*p
% subject to A*x >= b
% where
% H = [1 0 0; 0 1 0; 0 0 0.01]
% p = [2 1 -1]'
% A = [0 1 1]
% b = 1
% Use the initial point
% x0 = [1 1 1]'
% Solution:
% Execute the following commands:
% H = [1 0 0; 0 1 0; 0 0 0.01]
% p = [2 1 -1]'
% A = [0 1 1]
% b = 1
% x0 = [1 1 1]'
% [xs,fs,k] = qp_ie0(H,p,A,b,x0)
% =====================================
function [xs,fs,k] = qp_ie0(H,p,A,b,x0)
    disp(' ')
    disp('Program qp_ie0.m')
    % Data initialization.
    n = length(p);
    p1 = length(b);
    x = x0(:);
    k = 0;
    cm = 1;
    r = A*x-b;
    cha = ones(p1,1);
    zr = zeros(p1,1);
    Aa = [];
    epsilon = 1e-15;
    % Form active cconstraint matrix for x0.
    ind1 = find(abs(r) < epsilon);
    cn = length(ind1);
    Aa = [Aa; A(ind1,:)];
    cha(ind1) = 0;
    % Iteratiion begins.
    while cm > 0,
        flag = 0;
        % Check rank consistency condition (13.21).
        g = H*x+p;
        r_a = rank(Aa);
        r_ag = rank([Aa' g]);
        % If (13.21) holds, then compute Lagrange multiplier mu.
        if r_a == r_ag
            mu=pinv(Aa')*g;
            % Check if mu >= 0.
            ind2 = find(mu < -epsilon);
            cm = length(ind2);
            % If mu contains negative entries, remove one active constraint
            % from the active set.
        if cm > 0
            flag = 1;
            [y1,I1] = min(mu);
            wk = ind1(I1);
            cha(wk) = 1;
            ind1 = find(cha == 0);
            Aa = A(ind1,:);
            cn = cn-1;
        end
    end
    % Solve QP subproblem with equality constraints.
    if flag == 1 | r_ag - r_a > 0
        d = qp_e(H+1e-14*eye(n),g,Aa,zeros(cn,1));
        % Determine step length alpha_k.
        adi = A*d;
        i3 = find(adi < -epsilon);
        ic = zr;
        ic(i3) = 1;
        I3 = ic.*cha;
        ind3 = find(I3 == 1);
        if ~isempty(ind3)
            alp = (b(ind3) - A(ind3,:)*x)./adi(ind3);
        else alp = [];
        end
        [y2,I2] = min(alp);
        if y2 <= 1,
            is = ind3(I2);
            Aa = [Aa;A(is,:)];
            ind1 = [ind1; is];
            cha(is) = 0;
            cn = cn+1;
            alpha = y2;
        else
            alpha = 1;
        end
        % Set next iterate and update iteration number k.
        x = x + alpha*d;
    end
    k = k + 1;
end
xs = x;
fs = 0.5*xs'*(H*xs + 2*p);
end