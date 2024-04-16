function [xV,zV,z1v] = HIF_fucn(n,Q,R,z,x,P,xV,zV,zIV,z1v,L,S,A,B,h,N,I)
    
 

    for k = 1:N
        %if I(k) ~=0
           %I(k) = I(k)+ e(k);
        %end
        zz=  z(k);
        f = @(x)(A*x-B*I(k));
        %zz=  z(k);
        % 状态测量值,电压值Uk
        zV(:,k) =zz;
        % 电流激励I
        zIV(:,k) = I(k);
        x1 = f(x);
        % 过程方差预测，+-
        P = A*P*A'+Q;
        % 计算h的雅可比矩阵，，雅可比行列式规则
        [z1,C] = jaccsd(h,x1);
        z1v(:,k) = z1;
        % 卡尔曼增益，对应line4
        % inv返回逆矩阵
        % K = P*C'/(C*P*C'+R);   %4*1/1*1
        %K = A*P/(eye(n)-L*S*P+C'/(R)*C*P)*C'/(R);
        K = A*P/(eye(n)-L*S*P+C'/R*C*P)*C'/(R);
        % 状态估计值，对应line5
        x = x1+K.*(zz-z1);
        % EKF方差，对应line6
        %P = P-K*C*P;
        %P = P*inv(eye(n)-L*S*P+C'*inv(R)*C*P);
        P = P*inv(eye(n)-L*S*P+C'/R*C*P);
        % save
        xV(:,k) = x;
end