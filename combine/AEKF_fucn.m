function [xV,z,z1v] = AEKF_fucn(n,Q,R,q0,r0,I,z,x,P,xV,zV,zIV,z1v,b,T,A,B,h,N)

    for k = 1:N
        zz=  z(k);
        f = @(x)(A*x-B*I(k));

        % 状态测量值,电压值Uk
        zV(:,k) = zz;
        % 电流激励I
        zIV(:,k) = I(k);


        % 这里不一样

        %x1 = f(x);
        %x1 = f(x)+q0;    %下一步估计
        x1 = f(x);                 %之后的更新噪音用
        % 过程方差预测，+-
        P = A*P*A'+Q;
        %P = A*P*A'+T*Q*T';
        % 计算h的雅可比矩阵，，雅可比行列式规则
        [z1,C] = jaccsd(h,x1);
        z1v(:,k) = z1;
        % 卡尔曼增益，对应line4
        % inv返回逆矩阵
        K = P*C'/(C*P*C'+R);   %4*1/1*1

        % 自适应噪声协方差更新,这里M取2
        %Hk = ((z(k)-z1)*(z(k)-z1)'+(z(k)-z1)*(z(k)-z1)')/2;
        E = zz-z1;  %残差


         % EKF方差，对应line6

        % 状态EKF估计值，对应line5
        %x = x1+K*(z(k)-z1);
        x = x1+K*E;
        P = P-K*C*P;
    %% 这里的q,r,Q,R更新，用到的不是残差，注意区别
    % 参考  基于改进的Sage-Husa自适应扩展卡尔曼滤波的车辆状态估计，李刚
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