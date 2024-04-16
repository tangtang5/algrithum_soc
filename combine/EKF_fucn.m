function [xV,z,z1v] = EKF_fucn(n,Q,R,I,z,x,P,xV,zV,zIV,z1v,A,B,h,N)
    
    for k = 1:N
        %if I(k) ~=0
            %I(k) = I(k)+ e(k);
        %end
        
        zz=  z(k);
        f = @(x)(A*x-B*I(k));
        
        % ״̬����ֵ,��ѹֵUk,�����������
        zV(:,k) = zz;
        % ��������I
        zIV(:,k) = I(k);
        x1 = f(x);
        % ���̷���Ԥ�⣬+-
        P = A*P*A'+Q;
        % ����h���ſɱȾ��󣬣��ſɱ�����ʽ����
        [z1,C] = jaccsd(h,x1);
        z1v(:,k) = z1;
        %C =[-1 -1 -1 5*4.715*(x(4))^4-4*14.53*(x(4))^3+3*16.32*(x(4))^2-2*8.031*(x(4))+2.438];
        % ���������棬��Ӧline4
        % inv���������
        K = P*C'/(C*P*C'+R);   %4*1/1*1
        % ״̬EKF����ֵ����Ӧline5
        x = x1+K*(zz-z1);
        % EKF�����Ӧline6
        P = P-K*C*P;

        % save
        xV(:,k) = x;
        % update process ���Կ��ǼӲ��������������Ϊ�����������

    end
end