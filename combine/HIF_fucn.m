function [xV,zV,z1v] = HIF_fucn(n,Q,R,z,x,P,xV,zV,zIV,z1v,L,S,A,B,h,N,I)
    
 

    for k = 1:N
        %if I(k) ~=0
           %I(k) = I(k)+ e(k);
        %end
        zz=  z(k);
        f = @(x)(A*x-B*I(k));
        %zz=  z(k);
        % ״̬����ֵ,��ѹֵUk
        zV(:,k) =zz;
        % ��������I
        zIV(:,k) = I(k);
        x1 = f(x);
        % ���̷���Ԥ�⣬+-
        P = A*P*A'+Q;
        % ����h���ſɱȾ��󣬣��ſɱ�����ʽ����
        [z1,C] = jaccsd(h,x1);
        z1v(:,k) = z1;
        % ���������棬��Ӧline4
        % inv���������
        % K = P*C'/(C*P*C'+R);   %4*1/1*1
        %K = A*P/(eye(n)-L*S*P+C'/(R)*C*P)*C'/(R);
        K = A*P/(eye(n)-L*S*P+C'/R*C*P)*C'/(R);
        % ״̬����ֵ����Ӧline5
        x = x1+K.*(zz-z1);
        % EKF�����Ӧline6
        %P = P-K*C*P;
        %P = P*inv(eye(n)-L*S*P+C'*inv(R)*C*P);
        P = P*inv(eye(n)-L*S*P+C'/R*C*P);
        % save
        xV(:,k) = x;
end