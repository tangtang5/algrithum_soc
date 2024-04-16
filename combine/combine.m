% https://blog.csdn.net/zengxiantao1994/article/details/71170728
% function: simulating the process of EKF
clc;
clear;
clear all;

load('BDDST2.mat');
%load('BDDST2.mat');
%% 共有初始化
n = 4;             % 状态维度
N =length(t);
Q = 0.000000000001;          % 过程方差  
%Q =  0.0000000001;  
W = sqrt(Q)*randn(n,N);          
R = 0.01;           % 测量方差
V = sqrt(R)*randn(1,N);
s = [0;0;0;1];   
%x = [0;0;0;0.4]; 
x = [0;0;0;0.6]; 
% 4行50列，一列代表一个数据
% 最优估计值
xV = zeros(n,N);
% 真实值
sV = zeros(n,N);
% 状态测量值
zV = zeros(n,N);
% 电流激励
zIV = zeros(n,N);

z1v = zeros(1,N);


%% ekf初始化
% eye返回单位矩阵
P_ekf = eye(n);

%% AEKF初始化
%r_aekf = 1;   % 测量标准差 初始值,为了好改值，不然直接用数字
q0 = zeros(n,1);%更新的状态噪音 
r0 = 1;    %更新的测量噪音
b = 0.97;     % sage-husa更新用
T = 0.01*eye(n);
P_aekf = eye(n);

%% hif初始化 注意P矩阵形式，之前是没有注意，结果老是会出问题
L = 0.1;             %倒数是边界,但是L越大，K越小，失去纠偏能力
%S = [0.2 0 0 0;0 0.2 0 0;0 0 0.2 0;0 0 0 0.5];                %未知 
S=1;
P_hif = eye(n);
%% PF初始化
%PP = 0.00000001;%初始状态服从0均值，sqrt(PP）方差高斯分布
NN = 100;
%s = [0;0;0;1];                        % 初始状态 U0 U1 U2 SoC
Nth =100;  %重采样阈值
xhatPart = x;
%PP = 1;%初始状态服从0均值，sqrt(PP）方差高斯分布
P_pf = eye(n);
xpart = zeros(n,NN);
xhatPartArr = zeros(n,N); %存储估计的SOC值
xpartminus = zeros(n,NN);
for i = 1 : NN
    xpart(:,i) = x + sqrt(Q)* randn(n,1);
end
%% EKPF初始化
R_ekpf =2;
ekpfxhatPart = x;
ekpfxpart = zeros(n,NN);
ekpfp = zeros(n,n,NN);
ekpfxhatPartArr = zeros(n,N); %存储估计的SOC值
z1v_ekpf = zeros(1,N);
for i = 1 : NN
    ekpfxpart(:,i) = x + sqrt(Q)* randn(n,1);
    ekpfp(:,:,i) = eye(n);
end
ekpfpupdate = ekpfp;
ekpfxpartupdate = ekpfxpart;

%% 初始函数
R0 = @(x)(-0.07495*(x(4))^4+0.2187*(x(4))^3-0.1729*(x(4))^2+0.01904*(x(4))+0.1973);
R1 = @(x)(0.07826*(x(4))^4-0.2208*(x(4))^3+0.217*(x(4))^2-0.08761*(x(4))+0.01664);
R2 = @(x)(0.1248*(x(4))^4-0.2545*(x(4))^3+0.1254*(x(4))^2-0.03868*(x(4))+0.05978);
C1 = @(x)(2431*(x(4))^4-4606*(x(4))^3+3084*(x(4))^2-589*(x(4))+209.8);
C2 = @(x)(681.1*(x(4))^4-3197*(x(4))^3+4595*(x(4))^2-3114*(x(4))+1444);
U_HPPC_OCV = @(x)(4.715*(x(4))^5-14.53*(x(4))^4+16.32*(x(4))^3-8.031*(x(4))^2+2.438*(x(4))+3.247);
U_OCV = @(x)(0.008488*(x(4))^5-0.03224*(x(4))^4+0.004118*(x(4))^3+0.06053*(x(4))^2+0.2167*(x(4))+3.683);
% 系统测量方程
%f = Ax+bi+w        % 状态方程是线性的系统
h = @(x)(4.715*(x(4))^5-14.53*(x(4))^4+16.32*(x(4))^3-8.031*(x(4))^2+2.438*(x(4))+3.247-x(3)-x(2)-x(1));                         % 测量方程

% 系数矩阵

% 因为放电为负号，所以这里估算的时候B的最后一个变成正，单位是A*s[安培*秒]
%A = [0 0 0 0;0 exp(-1/(R1*C1)) 0 0;0 0 exp(-1/(R2*C2)) 0;0 0 0 1];
A = [0 0 0 0;0 exp(-1/(R1(x)*C1(x))) 0 0;0 0 exp(-1/(R2(x)*C2(x))) 0;0 0 0 1];    
%A = [0 0 0 0;0 (R1*C1/(1+R1*C1)) 0 0;0 0 (R2*C2/(1+R2*C2)) 0;0 0 0 1];
%B = [R0;R1*(1-exp(-1/(R1*C1)));R2*(1-exp(-1/(R2*C2)));-1/3.18/3600];
%B = [R0;R1*(1-exp(-1/(R1*C1)));R2*(1-exp(-1/(R2*C2)));-1/3.18/3600];
B = [R0(x);R1(x)*(1-exp(-1/(R1(x)*C1(x))));R2(x)*(1-exp(-1/(R2(x)*C2(x))));-1/3.18/3600];
As = [0 0 0 0;0 exp(-1/(R1(s)*C1(s))) 0 0;0 0 exp(-1/(R2(s)*C2(s))) 0;0 0 0 1];
Bs = [R0(s);R1(s)*(1-exp(-1/(R1(s)*C1(s))));R2(s)*(1-exp(-1/(R2(s)*C2(s))));-1/3.18/3600];


%% 真实值估计
sV = sV_fucn(As,Bs,s,sV,N,I,Q,n);

%for i =1:N
    %I(k) = I(k)+0.35;
%end

%% 加入算法
%xV是估计值，sV是真实值
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

%%  误差分析
% RMSE
%RMSE_AH = sqrt(sum((soc-sV(4,:)).^2)/N);
RMSE_EKF = sqrt(sum((xV(4,:)-sV(4,:)).^2)/N);
RMSE_AEKF = sqrt(sum((xV_aekf(4,:)-sV(4,:)).^2)/N);
RMSE_HIF = sqrt(sum((xV_hif(4,:)-sV(4,:)).^2)/N);
RMSE_PF = sqrt(sum((xV_pf(4,:)-sV(4,:)).^2)/N);
RMSE = [RMSE_EKF,RMSE_AEKF,RMSE_HIF,RMSE_PF]';
%RMSE_EKPF = sqrt(sum((xV_ekpf(4,:)-sV(4,:)).^2)/N);
%%   MAE 平均绝对误差
%MAE_AH = sum(abs(soc-sV(4,:)))/N;
MAE_EKF = sum(abs(xV(4,:)-sV(4,:)))/N;
MAE_AEKF = sum(abs(xV_aekf(4,:)-sV(4,:)))/N;
MAE_HIF = sum(abs(xV_hif(4,:)-sV(4,:)))/N;
MAE_PF = sum(abs(xV_pf(4,:)-sV(4,:)))/N;
MAE = [MAE_EKF,MAE_AEKF,MAE_HIF,MAE_PF]';
%MAE_EKPF = sum(abs(xV_ekpf(4,:)-sV(4,:)))/N;
%% MAXE 最大误差
%MAXE_AH = max(abs(soc-sV(4,:)));
MAXE_EKF = max(abs(xV(4,:)-sV(4,:)));
MAXE_AEKF = max(abs(xV_aekf(4,:)-sV(4,:)));
MAXE_HIF = max(abs(xV_hif(4,:)-sV(4,:)));
MAXE_PF = max(abs(xV_pf(4,:)-sV(4,:)));
MAXE = [MAXE_EKF,MAXE_AEKF,MAXE_HIF,MAXE_PF]';
%MAXE_EKPF = max(abs(xV_ekpf(4,:)-sV(4,:)));
%% 绘制图形SOC

FontSize = 14;
LineWidth = 1;
figure(1);
% 画出真实值
plot(sV(4,:),'g');
hold on;
    
% 画出最优估计值
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


%legend('真实状态', 'EKF最优估计估计值','AEKF最优估计估计值','hif最优估计估计值','PF最优估计估计值');
legend('真实状态', 'EKF最优估计估计值','AEKF最优估计估计值','hif最优估计估计值');
xl = xlabel('t(s)');
xl = xlabel('t(s)');
% 把数值转换成字符串， 转换后可以使用fprintf或disp函数进行输出。
yl = ylabel('SoC(%)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 0 1.1]);
axis([0 25000 0 1.1]);
%% 绘制EKF估计端电压

figure(2);
% 画出状态测量值
plot(z1v(1,:),'k');
hold on;
plot(z,'r');
hold on;
legend('EKF端电压估计值','真实测量值');
xl = xlabel('t(s)');
% 把数值转换成字符串， 转换后可以使用fprintf或disp函数进行输出。
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
axis([0 25000 2.5 4.5]);
%axis([0 12000 2.5 5]);
%% 绘制AEKF估计端电压
figure(3);
% 画出状态测量值
plot(z1v_aekf(1,:),'k');
hold on;
plot(z,'r');
hold on;
legend('AEKF端电压估计值','真实测量值');
xl = xlabel('t(s)');
% 把数值转换成字符串， 转换后可以使用fprintf或disp函数进行输出。
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 2.5 5]);
axis([0 25000 2.5 4.5]);
%% 绘制hif估计端电压
figure(4);
% 画出状态测量值
plot(z1v_hif(1,:),'k');
hold on;
plot(z,'r');
hold on;
legend('HIF端电压估计值','真实测量值');
xl = xlabel('t(s)');
% 把数值转换成字符串， 转换后可以使用fprintf或disp函数进行输出。
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 2.5 5]);
axis([0 25000 2.5 4.5]);
%% 绘制PF估计端电压
figure(5);
% 画出状态测量值
plot(z1v_pf(1,:),'k');
hold on;
plot(z,'r');
hold on;
legend('PF端电压估计值','真实测量值');
xl = xlabel('t(s)');
% 把数值转换成字符串， 转换后可以使用fprintf或disp函数进行输出。
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 2.5 5]);
axis([0 25000 2.5 4.5]);
%% SOC误差分析
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

legend('EKF-SoC误差','AEKF-SoC误差','HIF-SoC误差');
%legend('EKF-SoC误差','AEKF-SoC误差','HIF-SoC误差','PF-SoC误差');
xl = xlabel('t(s)');
% 把数值转换成字符串， 转换后可以使用fprintf或disp函数进行输出。
yl = ylabel('SoC误差');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 25000 -0.06 0.06]);
%axis([0 12000 -0.005 0.03]);
axis([0 25000 -0.005 0.03]);
%% 绘制EKF端电压误差
figure(7);
% 画出状态测量值
plot(abs(z1v'-z),'k');
hold on;
legend('EKF估计端电压误差');
xl = xlabel('t(s)');
% 把数值转换成字符串， 转换后可以使用fprintf或disp函数进行输出。
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 0 0.1]);
axis([0 25000 0 0.18]);
%% 绘制AEKF端电压误差
figure(8);
% 画出状态测量值
plot(abs(z1v_aekf'-z),'k');
hold on;
legend('AEKF估计端电压误差');
xl = xlabel('t(s)');
% 把数值转换成字符串， 转换后可以使用fprintf或disp函数进行输出。
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 0 0.1]);
axis([0 25000 0 0.18]);
%% 绘制hif端电压误差
figure(9);
% 画出状态测量值
plot(abs(z1v_hif'-z),'k');
hold on;
legend('HIF估计端电压误差');
xl = xlabel('t(s)');
% 把数值转换成字符串， 转换后可以使用fprintf或disp函数进行输出。
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 0 0.1]);
axis([0 25000 0 0.18]);
%% 绘制Pf端电压误差
figure(10);
% 画出状态测量值
plot(abs(z1v_pf'-z),'k');
hold on;
legend('PF估计端电压误差');
xl = xlabel('t(s)');
% 把数值转换成字符串， 转换后可以使用fprintf或disp函数进行输出。
yl = ylabel('U(V)');
set(xl,'fontsize',FontSize);
set(yl,'fontsize',FontSize);
%axis([0 12000 0 0.1]);
axis([0 25000 0 0.18]);
%%
hold off;
set(gca,'FontSize',FontSize);
%% 这个只是说拟合的不错，还需要验证工况
