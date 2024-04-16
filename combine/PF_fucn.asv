function [xhatPartArr,z1v] = PF_fucn(n,Q,R,W,t,I,z,x,P,xpart,xhatPartArr,xhatPart,xpartminus,z1v,N,NN,Nth,A,B,h)

   
    for k = 1 : N       %tf为时间长度，k可以理解为时间轴上的k时刻
        zz=  z(k);


 
        

        % Particle filter 生成100个粒子并根据预测和观测值差值计算各个粒子的权重
        w = ones(1,NN);
        for i = 1 : NN
            %xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
            xpartminus(:,i) = f(xpart(:,i))+sqrt(Q)*randn(n,1);
            %ypart = xpartminus(i)^2 / 20;
            zpart = h(xpartminus(:,i));
            vhat = zz - zpart; %观测和预测的差
            w(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R); %根据差值给出权重
        end



        % Normalize the likelihood of each a priori estimate.

        wsum = sum(w);
        for i = 1 : NN
            w(i) = w(i) / wsum;%归一化权重
        end



       %% Resample.
       %参考 https://blog.csdn.net/u011624019/article/details/80559397
       %确定Neff
        ww=zeros(1,NN);
        %w_new = zeros(1,N);
        for i = 1:NN
            ww(i)=w(i)^2;
        end
        Neff = 1/sum(ww);
        w_new =w;

        %开始判断重采样
        if Neff<Nth
            for i = 1 : NN
                u = rand; % uniform random number between 0 and 1
                wtempsum = 0;
                for j = 1 : NN
                    wtempsum = wtempsum + w(j);
                    if wtempsum >= u
                    %重采样对低权重进行剔除，同时保留高权重，防止退化的办法
                        xpart(:,i) = xpartminus(:,j);
                        %w_new(j) = w(j);
                        break;          
                    end          
                end
            end
        end


       % The particle filter estimate is the mean of the particles.

        for i = 1:n
            xhatPart(i) = mean(xpart(i,:));%经过粒子滤波处理后的均值
            %mm=zeros(1,N);
            %for j=1:N
               % mm(j)=xpart(i,j)*w_new(j);
           % end
          %  xhatPart(i) = sum(mm);
        end
        
    
        xhatPartArr(:,k) = xhatPart;
        x1 = f(xhatPart);
        [z1,Ccc] = jaccsd(h,x1);
        z1v(:,k) = z1;
        
    end
end

