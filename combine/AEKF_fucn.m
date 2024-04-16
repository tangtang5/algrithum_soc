function [xV,z,z1v] = AEKF_fucn(n,Q,R,q0,r0,I,z,x,P,xV,zV,zIV,z1v,b,T,A,B,h,N)

    for k = 1:N
        zz=  z(k);
        f = @(x)(A*x-B*I(k));

        % ״̬����ֵ,��ѹֵUk
        zV(:,k) = zz;
        % ��������I
        zIV(:,k) = I(k);


        % ���ﲻһ��

        %x1 = f(x);
        %x1 = f(x)+q0;    %��һ������
        x1 = f(x);                 %֮��ĸ���������
        % ���̷���Ԥ�⣬+-
        P = A*P*A'+Q;
        %P = A*P*A'+T*Q*T';
        % ����h���ſɱȾ��󣬣��ſɱ�����ʽ����
        [z1,C] = jaccsd(h,x1);
        z1v(:,k) = z1;
        % ���������棬��Ӧline4
        % inv���������
        K = P*C'/(C*P*C'+R);   %4*1/1*1

        % ����Ӧ����Э�������,����Mȡ2
        %Hk = ((z(k)-z1)*(z(k)-z1)'+(z(k)-z1)*(z(k)-z1)')/2;
        E = zz-z1;  %�в�


         % EKF�����Ӧline6

        % ״̬EKF����ֵ����Ӧline5
        %x = x1+K*(z(k)-z1);
        x = x1+K*E;
        P = P-K*C*P;
    %% �����q,r,Q,R���£��õ��Ĳ��ǲвע������
    % �ο�  ���ڸĽ���Sage-Husa����Ӧ��չ�������˲��ĳ���״̬���ƣ����
        %G = inv(T'*T)*T';
        d = (1-b)/(1-b^(k));
        %q0 = (1-d)*q0+d*(x-x0);
        %Q = (1-d)*Q+d*G*(K*E*E'*K'+P-A*P*A')*G';
        %Q = (1-d)*Q+d*(K*E*E'*K');
        r0 = (1-d)*r0+d*(zz-z1);
        %R = (1-d)*R+d*(E*E'-C*P*C');
        R = (1-d)*R+d*(E*E');

    %%   
        % save

        xV(:,k) = x;
        % update process 


    end
end