function [ekpfxhatPartArr,zlv_ekpf] = EKPF_fucn(n,Q,R,I,z,x,ekpfp,ekpfxpart,ekpfxhatPartArr,ekpfxhatPart,ekpfxpartupdate,z1v,N,NN,Nth,A,B,h)

   
    for k = 1 : N       %tfΪʱ�䳤�ȣ�k�������Ϊʱ�����ϵ�kʱ��
        zz=  z(k);
        f = @(x)(A*x-B*I(k));
       
        x1 = f(x);
    %x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn; %״̬����(1)
    
    %y = x^2 / 20 + sqrt(R) * randn;%�۲ⷽ��(2)
    %P = A*P*A'+Q;
        [z1,Cc] = jaccsd(h,x1);
        zlv_ekpf(:,k) = z1;


        % Particle filter ����100�����Ӳ�����Ԥ��͹۲�ֵ��ֵ����������ӵ�Ȩ��
        w = ones(1,NN);
        for i = 1 : NN
            ekpfxpartupdate(:,i) = f(ekpfxpart(:,i))+sqrt(Q) * randn(n,1);
            ekpfp(:,:,i) = A * ekpfp(:,:,i) * A' + Q;% ����p
            [zpart,C] = jaccsd(h,ekpfxpartupdate(:,i));
            ekpfK = ekpfp(:,:,i)*C'/(C*ekpfp(:,:,i)*C'+R);
            ekpfxpartupdate(:,i) = ekpfxpartupdate(:,i) + ekpfK* (zz - zpart);
            ekpfp(:,:,i) = ekpfp(:,:,i) - ekpfK * C * ekpfp(:,:,i);
            zpart = h(ekpfxpartupdate(:,i))+sqrt(R)*randn;

            %����Ȩ��
            vhat = zz - zpart; %�۲��Ԥ��Ĳ�
            w(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R); %���ݲ�ֵ����Ȩ��
        end



        % Normalize the likelihood of each a priori estimate.

        wsum = sum(w);
        for i = 1 : NN
            w(i) = w(i) / wsum;%��һ��Ȩ��
        end



       %% Resample.
       %�ο� https://blog.csdn.net/u011624019/article/details/80559397
       %ȷ��Neff
        ww=zeros(1,NN);
        %w_new = zeros(1,N);
        for i = 1:NN
            ww(i)=w(i)^2;
        end
        Neff = 1/sum(ww);
        w_new =w;

        %��ʼ�ж��ز���
        if Neff<Nth
            for i = 1 : NN
                u = rand; % uniform random number between 0 and 1
                wtempsum = 0;
                for j = 1 : NN
                    wtempsum = wtempsum + w(j);
                    if wtempsum >= u
                    %�ز����Ե�Ȩ�ؽ����޳���ͬʱ������Ȩ�أ���ֹ�˻��İ취
                        ekpfxpart(:,i) = ekpfxpartupdate(:,j);
                        ekpfpupdate(:,:,j) = ekpfp(:,:,j);
                        %w_new(j) = w(j);
                        break;          
                    end          
                end
            end
        end


       % The particle filter estimate is the mean of the particles.

        for i = 1:n
            ekpfxhatPart(i) = mean(ekpfxpart(i,:));%���������˲������ľ�ֵ
            %mm=zeros(1,N);
            %for j=1:N
               % mm(j)=xpart(i,j)*w_new(j);
           % end
          %  xhatPart(i) = sum(mm);
        end
        
    
        ekpfxhatPartArr(:,k) = ekpfxhatPart;
    end
end

