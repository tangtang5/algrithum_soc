function [sV] = sV_fucn(A,B,s,sV,N,I,Q,n)
% Êµ¼Ê×´Ì¬
    for k = 1:N
        sV(:,k) = s;
        %s = f(s) + q*randn(n,1);
        s = A*s-B*I(k)+sqrt(Q)*randn(n,1);
    end



end