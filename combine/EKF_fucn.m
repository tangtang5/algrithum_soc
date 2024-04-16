function [xV,z,z1v] = EKF_fucn(n,Q,R,I,z,x,P,xV,zV,zIV,z1v,A,B,h,N)
    
    for k = 1:N
        %if I(k) ~=0
            %I(k) = I(k)+ e(k);
        %end
        
        zz=  z(k);
        f = @(x)(A*x-B*I(k));
        
        % 状态测量值,电压值Uk,加入测量噪音
        zV(:,k) = zz;
        % 电流激励I
        zIV(:,k) = I(k);
        x1 = f(x);
        % 过程方差预测，+-
        P = A*P*A'+Q;
        % 计算h的雅可比矩阵，，雅可比行列式规则
        [z1,C] = jaccsd(h,x1);
        z1v(:,k) = z1;
        %C =[-1 -1 -1 5*4.715*(x(4))^4-4*14.53*(x(4))^3+3*16.32*(x(4))^2-2*8.031*(x(4))+2.438];
        % 卡尔曼增益，对应line4
        % inv返回逆矩阵
        K = P*C'/(C*P*C'+R);   %4*1/1*1
        % 状态EKF估计值，对应line5
        x = x1+K*(zz-z1);
        % EKF方差，对应line6
        P = P-K*C*P;

        % save
        xV(:,k) = x;
        % update process 可以考虑加不加这个东西，人为的制造白噪音

    end
end