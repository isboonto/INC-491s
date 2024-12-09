function xs = qp_e(H,p,A,b)
    if min(size(A))== 0
        xs = -pinv(H)*p;
    else
        n = length(p);
        p1 = size(A)*[1 0]';
        [Q,R] = qr(A');
        R = R(1:p1,:);
        Q1 = Q(:,1:p1);
        Q2 = Q(:,p1+1:n);
        xs = Q1*((R')\b);
        Hh = Q2'*H*Q2;
        ph = Q2'*(H*xs+p);
        phi_s = -inv(Hh)*ph;
        xs = Q2*phi_s+xs;
    end
end