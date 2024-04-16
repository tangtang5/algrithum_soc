% https://blog.csdn.net/zengxiantao1994/article/details/71170728
% function: simulating the process of EKF
clc;
clear;
clear all;

load('BDDST2.mat');
%load('BDDST2.mat');
%% ���г�ʼ��
n = 4;             % ״̬ά��
N =length(t);
Q = 0.000000000001;          % ���̷���  
%Q =  0.0000000001;  
W = sqrt(Q)*randn(n,N);          
R = 0.01;           % ��������
V = sqrt(R)*randn(1,N);
s = [0;0;0;1];   
%x = [0;0;0;0.4]; 
x = [0;0;0;0.6]; 
% 4��50�У�һ�д���һ������
% ���Ź���ֵ
xV = zeros(n,N);
% ��ʵֵ
sV = zeros(n,N);
% ״̬����ֵ
zV = zeros(n,N);
% ��������
zIV = zeros(n,N);

z1v = zeros(1,N);


%% ekf��ʼ��
% eye���ص�λ����
P_ekf = eye(n);

%% AEKF��ʼ��
%r_aekf = 1;   % ������׼�� ��ʼֵ,Ϊ�˺ø�ֵ����Ȼֱ��������
q0 = zeros(n,1);%���µ�״̬���� 
r0 = 1;    %���µĲ�������
b = 0.97;     % sage-husa������
T = 0.01*eye(n);
P_aekf = eye(n);

%% hif��ʼ�� ע��P������ʽ��֮ǰ��û��ע�⣬������ǻ������
L = 0.1;             %�����Ǳ߽�,����LԽ��KԽС��ʧȥ��ƫ����
%S = [0.2 0 0 0;0 0.2 0 0;0 0 0.2 0;0 0 0 0.5];                %δ֪ 
S=1;
P_hif = eye(n);
%% PF��ʼ��
%PP = 0.00000001;%��ʼ״̬����0��ֵ��sqrt(PP�������˹�ֲ�
NN = 100;
%s = [0;0;0;1];                        % ��ʼ״̬ U0 U1 U2 SoC
Nth =100;  %�ز�����ֵ
xhatPart = x;
%PP = 1;%��ʼ״̬����0��ֵ��sqrt(PP�������˹�ֲ�
P_pf = eye(n);
xpart = zeros(n,NN);
xhatPartArr = zeros(n,N); %�洢���Ƶ�SOCֵ
xpartminus = zeros(n,NN);
for i = 1 : NN
    xpart(:,i) = x + sqrt(Q)* randn(n,1);
end
%% EKPF��ʼ��
R_ekpf =2;
ekpfxhatPart = x;
ekpfxpart = zeros(n,NN);
ekpfp = zeros(n,n,NN);
ekpfxhatPartArr = zeros(n,N); %�洢���Ƶ�SOCֵ
z1v_ekpf = zeros(1,N);
for i = 1 : NN
    ekpfxpart(:,i) = x + sqrt(Q)* randn(n,1);
    ekpfp(:,:,i) = eye(n);
end
ekpfpupdate = ekpfp;
ekpfxpartupdate = ekpfxpart;

%% ��ʼ����
R0 = @(x)(-0.07495*(x(4))^4+0.2187*(x(4))^3-0.1729*(x(4))^2+0.01904*(x(4))+0.1973);
R1 = @(x)(0.07826*(x(4))^4-0.2208*(x(4))^3+0.217*(x(4))^2-0.08761*(x(4))+0.01664);
R2 = @(x)(0.1248*(x(4))^4-0.2545*(x(4))^3+0.1254*(x(4))^2-0.03868*(x(4))+0.05978);
C1 = @(x)(2431*(x(4))^4-4606*(x(4))^3+3084*(x(4))^2-589*(x(4))+209.8);
C2 = @(x)(681.1*(x(4))^4-3197*(x(4))^3+4595*(x(4))^2-3114*(x(4))+1444);
U_HPPC_OCV = @(x)(4.715*(x(4))^5-14.53*(x(4))^4+16.32*(x(4))^3-8.031*(x(4))^2+2.438*(x(4))+3.247);
U_OCV = @(x)(0.008488*(x(4))^5-0.03224*(x(4))^4+0.004118*(x(4))^3+0.06053*(x(4))^2+0.2167*(x(4))+3.683);
% ϵͳ��������
%f = Ax+bi+w        % ״̬���������Ե�ϵͳ
h = @(x)(4.715*(x(4))^5-14.53*(x(4))^4+16.32*(x(4))^3-8.031*(x(4))^2+2.438*(x(4))+3.247-x(3)-x(2)-x(1));                         % ��������

% ϵ������

% ��Ϊ�ŵ�Ϊ���ţ�������������ʱ��B�����һ�����������λ��A*s[����*��]
%A = [0 0 0 0;0 exp(-1/(R1*C1)) 0 0;0 0 exp(-1/(R2*C2)) 0;0 0 0 1];
A = [0 0 0 0;0 exp(-1/(R1(x)*C1(x))) 0 0;0 0 exp(-1/(R2(x)*C2(x))) 0;0 0 0 1];    
%A = [0 0 0 0;0 (R1*C1/(1+R1*C1)) 0 0;0 0 (R2*C2/(1+R2*C2)) 0;0 0 0 1];
%B = [R0;R1*(1-exp(-1/(R1*C1)));R2*(1-exp(-1/(R2*C2)));-1/3.18/3600];
%B = [R0;R1*(1-exp(-1/(R1*C1)));R2*(1-exp(-1/(R2*C2)));-1/3.18/3600];
B = [R0(x);R1(x)*(1-exp(-1/(R1(x)*C1(x))));R2(x)*(1-exp(-1/(R2(x)*C2(x))));-1/3.18/3600];
As = [0 0 0 0;0 exp(-1/(R1(s)*C1(s))) 0 0;0 0 exp(-1/(R2(s)*C2(s))) 0;0 0 0 1];
Bs = [R0(s);R1(s)*(1-exp(-1/(R1(s)*C1(s))));R2(s)*(1-exp(-1/(R2(s)*C2(s))));-1/3.18/3600];


%% ��ʵֵ����
sV = sV_fucn(As,Bs,s,sV,N,I,Q,n);

%for i =1:N
    %I(k) = I(k)+0.35;
%end

%% �����㷨
%xV�ǹ���ֵ��sV����ʵֵ
[xV,z,z1v] = EKF_fucn(n,Q,R,I,z,x,P_ekf,xV,zV,zIV,z1v,A,B,h,N);
[xV_aekf,z_aekf,z1v_aekf] = AEKF_fucn(n,Q,R,q0,r0,I,z,x,P_aekf,xV,zV,zIV,z1v,b,T,A,B,h,N);
[xV_hif,z_hif,z1v_hif] = HIF_fucn(n,Q,R,z,x,P_hif,xV,zV,zIV,z1v,L,S,A,B,h,N,I);
[xV_pf,z1v_pf] = PF_fucn(n,Q,R,W,t,I,z,x,P_pf,xpart,xhatPartArr,xhatPart,xpartminus,z1v,N,NN,Nth,A,B,h);
%[xV_ekpf,zlv_ekpf] = EKPF_fucn(n,Q,R_ekpf,I,z,x,ekpfp,ekpfxpart,ekpfxhatPartArr,ekpfxhatPart,ekpfxpartupdate,z1v,N,NN,Nth,A,B,h);




soc = zeros(1,N);
soc(1) = 1;
for k = 2:N
    soc(k) = soc(k-1)+I(k)/3.18/3600;
end

%%  ������
% RMSE
%RMSE_AH = sqrt(sum((soc-sV(4,:)).^2)/N);
RMSE_EKF = sqrt(sum((xV(4,:)-sV(4,:)).^2)/N);
RMSE_AEKF = sqrt(sum((xV_aekf(4,:)-sV(4,:)).^2)/N);
RMSE_HIF = sqrt(sum((xV_hif(4,:)-sV(4,:)).^2)/N);
RMSE_PF = sqrt(sum((xV_pf(4,:)-sV(4,:)).^2)/N);
RMSE = [RMSE_EKF,RMSE_AEKF,RMSE_HIF,RMSE_PF]';
%RMSE_EKPF = sqrt(sum((xV_ekpf(4,:)-sV(4,:)).^2)/N);
%%   MAE ƽ���������
%MAE_AH = sum(abs(soc-sV(4,:)))/N;
MAE_EKF = sum(abs(xV(4,:)-sV(4,:)))/N;
MAE_AEKF = sum(abs(xV_aekf(4,:)-sV(4,:)))/N;
MAE_HIF = sum(abs(xV_hif(4,:)-sV(4,:)))/N;
MAE_PF = sum(abs(xV_pf(4,:)-sV(4,:)))/N;
MAE = [MAE_EKF,MAE_AEKF,MAE_HIF,MAE_PF]';
%MAE_EKPF = sum(abs(xV_ekpf(4,:)-sV(4,:)))/N;
%% MAXE ������
%MAXE_AH = max(abs(soc-sV(4,:)));
MAXE_EKF = max(abs(xV(4,:)-sV(4,:)));
MAXE_AEKF = max(abs(xV_aekf(4,:)-sV(4,:)));
MAXE_HIF = max(abs(xV_hif(4,:)-sV(4,:)));
MAXE_PF = max(abs(xV_pf(4,:)-sV(4,:)));
MAXE = [MAXE_EKF,MAXE_AEKF,MAXE_HIF,MAXE_PF]';
%MAXE_EKPF = max(abs(xV_ekpf(4,:)-sV(4,:)));
%% ����ͼ��SOC

FontSize = 14;
LineWidth = 1;
figure(1);
% ������ʵֵ
plot(sV(4,:),'g');
hold on;
    
% �������Ź���ֵ
plot(xV(4,:),'LineWidth',LineWidth+1);
hold on;

plot(xV_aekf(4,:),'LineWidth',LineWidth);
hold on;

plot(xV_hif(4,:),'LineWidth',LineWidth);
hold on;
%plot(xV_pf(4,:),'LineWidth',LineWidth);
%hold on;

%plot(soc,'LineWidth',LineWidth);
%hold on;
%plot(xV_ekpf(4,:),'LineWidth',LineWidth);
%hold on;


%legend('��ʵ״̬', 'EKF���Ź��ƹ���ֵ','AEKF���Ź��ƹ���ֵ','hif���Ź��ƹ���ֵ','PF���Ź��ƹ���ֵ');
legend('��ʵ״̬', 'EKF���Ź��ƹ���ֵ','AEKF���Ź��ƹ���ֵ','hif���Ź��ƹ���ֵ');
xl = xlabel('t(s)');
xl = xlabel('t(s)');
% ����ֵת�����ַ����� ת�������ʹ��fprintf��disp�������������
yl = ylabel('SoC(%)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 0 1.1]);
axis([0 25000 0 1.1]);
%% ����EKF���ƶ˵�ѹ

figure(2);
% ����״̬����ֵ
plot(z1v(1,:),'k');
hold on;
plot(z,'r');
hold on;
legend('EKF�˵�ѹ����ֵ','��ʵ����ֵ');
xl = xlabel('t(s)');
% ����ֵת�����ַ����� ת�������ʹ��fprintf��disp�������������
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
axis([0 25000 2.5 4.5]);
%axis([0 12000 2.5 5]);
%% ����AEKF���ƶ˵�ѹ
figure(3);
% ����״̬����ֵ
plot(z1v_aekf(1,:),'k');
hold on;
plot(z,'r');
hold on;
legend('AEKF�˵�ѹ����ֵ','��ʵ����ֵ');
xl = xlabel('t(s)');
% ����ֵת�����ַ����� ת�������ʹ��fprintf��disp�������������
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 2.5 5]);
axis([0 25000 2.5 4.5]);
%% ����hif���ƶ˵�ѹ
figure(4);
% ����״̬����ֵ
plot(z1v_hif(1,:),'k');
hold on;
plot(z,'r');
hold on;
legend('HIF�˵�ѹ����ֵ','��ʵ����ֵ');
xl = xlabel('t(s)');
% ����ֵת�����ַ����� ת�������ʹ��fprintf��disp�������������
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 2.5 5]);
axis([0 25000 2.5 4.5]);
%% ����PF���ƶ˵�ѹ
figure(5);
% ����״̬����ֵ
plot(z1v_pf(1,:),'k');
hold on;
plot(z,'r');
hold on;
legend('PF�˵�ѹ����ֵ','��ʵ����ֵ');
xl = xlabel('t(s)');
% ����ֵת�����ַ����� ת�������ʹ��fprintf��disp�������������
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 2.5 5]);
axis([0 25000 2.5 4.5]);
%% SOC������
figure(6);
plot(abs(xV(4,:)-sV(4,:)),'LineWidth',LineWidth+1);
hold on;


plot(abs(xV_aekf(4,:)-sV(4,:)),'LineWidth',LineWidth);
hold on;

plot(abs(xV_hif(4,:)-sV(4,:)),'LineWidth',LineWidth);
hold on;
%plot(abs(xV_pf(4,:)-sV(4,:)),'LineWidth',LineWidth);
%hold on;

%plot(abs(soc-sV(4,:)),'LineWidth',LineWidth);
%hold on;

%plot(abs(xV_ekpf(4,:)-sV(4,:)),'LineWidth',LineWidth);
%hold on;

legend('EKF-SoC���','AEKF-SoC���','HIF-SoC���');
%legend('EKF-SoC���','AEKF-SoC���','HIF-SoC���','PF-SoC���');
xl = xlabel('t(s)');
% ����ֵת�����ַ����� ת�������ʹ��fprintf��disp�������������
yl = ylabel('SoC���');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 25000 -0.06 0.06]);
%axis([0 12000 -0.005 0.03]);
axis([0 25000 -0.005 0.03]);
%% ����EKF�˵�ѹ���
figure(7);
% ����״̬����ֵ
plot(abs(z1v'-z),'k');
hold on;
legend('EKF���ƶ˵�ѹ���');
xl = xlabel('t(s)');
% ����ֵת�����ַ����� ת�������ʹ��fprintf��disp�������������
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 0 0.1]);
axis([0 25000 0 0.18]);
%% ����AEKF�˵�ѹ���
figure(8);
% ����״̬����ֵ
plot(abs(z1v_aekf'-z),'k');
hold on;
legend('AEKF���ƶ˵�ѹ���');
xl = xlabel('t(s)');
% ����ֵת�����ַ����� ת�������ʹ��fprintf��disp�������������
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 0 0.1]);
axis([0 25000 0 0.18]);
%% ����hif�˵�ѹ���
figure(9);
% ����״̬����ֵ
plot(abs(z1v_hif'-z),'k');
hold on;
legend('HIF���ƶ˵�ѹ���');
xl = xlabel('t(s)');
% ����ֵת�����ַ����� ת�������ʹ��fprintf��disp�������������
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 0 0.1]);
axis([0 25000 0 0.18]);
%% ����Pf�˵�ѹ���
figure(10);
% ����״̬����ֵ
plot(abs(z1v_pf'-z),'k');
hold on;
legend('PF���ƶ˵�ѹ���');
xl = xlabel('t(s)');
% ����ֵת�����ַ����� ת�������ʹ��fprintf��disp�������������
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 0 0.1]);
axis([0 25000 0 0.18]);
%%
hold off;
set(gca,'FontSize',FontSize);
%% ���ֻ��˵��ϵĲ�������Ҫ��֤����
