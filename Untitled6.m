%   自适应多模型Kalman滤波的目标跟踪仿真

function Untitled6
clear;
T=2;%雷达扫描周期，即采样周期
M=20;%滤波次数
N=900/T;%总采样点数
N1=400/T;%第一转弯处采样起点
N2=600/T;%第一匀速处采样起点
N3=610/T;%第二转弯处采样起点
N4=660/T;%第二匀速处采样起点
Delta=30;%测量噪声标准差

%对真实的轨迹和观测轨迹数据的初始化，以向右向上为正
Rx=zeros(N,1);
Ry=zeros(N,1);
Zx=zeros(N,M);
Zy=zeros(N,M);
% 1-沿y轴匀速
t=2:T:400;
x0=2000+0*t';
y0=0+10*t';%初始速度Vy=10，位置（2000，0）
% 2-慢转弯
t=402:T:600;
x1=x0(N1)+0.075*((t'-400).^2)/2;%X轴正向加速度0.075
y1=y0(N1)+10*(t'-400)-0.075*((t'-400).^2)/2;%Y轴反向加速度0.075
% 3-沿X轴匀速
t=602:T:610;
vx=0.075*(600-400);
x2=x1(N2-N1)+vx*(t'-600);
y2=y1(N2-N1)+0*t';%此时Y轴速度正好为0，因此Y位置不变
% 4-快转弯
t=612:T:660;
x3=x2(N3-N2)+vx*(t'-610)-0.3*((t'-610).^2)/2;%X轴反向加速度0.3
y3=y2(N3-N2)-0.3*((t'-610).^2)/2;%Y轴反向加速度0.3
% 5-沿Y轴匀速
t=662:T:900;
vy=-0.3*(660-610);
x4=x3(N4-N3)+0*t';%此时X轴速度正好为0
y4=y3(N4-N3)+vy*(t'-660);
% 最终将所有轨迹整合成为一个列向量，即真实轨迹，Rx为Real-x,Ry为Real-y
Rx=[x0;x1;x2;x3;x4];
Ry=[y0;y1;y2;y3;y4];
% 对每次蒙特卡洛仿真的滤波估计位置的初始化
Mt_Est_Px=zeros(M,N);
Mt_Est_Py=zeros(M,N);
% 产生观测数据，要仿真M次，必须有M次
nx=randn(N,M)*Delta;%产生X轴和Y轴观测噪声
ny=randn(N,M)*Delta;
Zx=Rx*ones(1,M)+nx;%真实值的基础上叠加噪声，即构成仿真的观测值
Zy=Ry*ones(1,M)+ny;

for m=1:M
    %滤波初始化
    Mt_Est_Px(m,1)=Zx(1,m);%初始数据
    Mt_Est_Py(m,1)=Zx(2,m);
    xn(1)=Zx(1,m);%滤波初值
    xn(2)=Zx(2,m);
    yn(1)=Zy(1,m);
    yn(2)=Zy(2,m);
    %非机动模型参数
    phi=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];%状态转移矩阵，没有加速度
    h=[1,0,0,0;0,0,1,0];%观测矩阵
    g=[T/2,0;1,0;0,T/2;0,1];%噪声增益矩阵
    q=[Delta^2,0;0,Delta^2];%过程噪声协方差
    vx=(Zx(2)-Zx(1,m))/2;%非机动模型速度不变，所以每次滤波按位移比时间设定
    vy=(Zy(2)-Zy(1,m))/2;
    %初始状态估计
    x_est=[Zx(1,m);vx;Zy(1,m);vy];
    p_est=[Delta^2,Delta^2/T,0,0;Delta^2/T,2*Delta^2/(T^2),0,0;
        0,0,Delta^2,Delta^2/T;0,0,Delta^2/T,2*Delta^2/(T^2)];
    Mt_Est_Px(m,2)=x_est(1);%把初始估计值赋值给M仿真
    Mt_Est_Py(m,2)=x_est(3);
    %滤波开始
    for r=1:N
        z=[Zx(r,m);Zy(r,m)];
        if r<10%前10次不做IMM，只是单纯的线性卡尔曼滤波作为对比
            x_pre=phi*x_est;%预测
            p_pre=phi*p_est*phi';%预测误差协方差
            k=p_pre*h'*inv(h*p_pre*h'+q);%卡尔曼增益
            x_est=x_pre+k*(z-h*x_pre);%滤波
            p_est=(eye(4)-k*h)*p_pre;%滤波协方差
            xn(r)=x_est(1);%记录采样点滤波数据
            yn(r)=x_est(3);
            Mt_Est_Px(m,r)=x_est(1);%记录第m次仿真滤波估计数据
            Mt_Est_Py(m,r)=x_est(3);
        else
            if r==10
                X_est=[x_est;0;0];%扩充维数，引入加速度
                P_est=p_est;
                P_est(6,6)=0;
                for i=1:3
                    Xn_est{i,1}=X_est;
                    Pn_est{i,1}=P_est;
                end
                u=[0.4,0.3,0.3];%模型概率初始化
            end
            %调用IMM算法
            [X_est,P_est,Xn_est,Pn_est,u]=IMM(Xn_est,Pn_est,T,z,Delta,u);
            xn(r)=X_est(1);
            yn(r)=X_est(3);
            Mt_Est_Px(m,r)=X_est(1);
            Mt_Est_Py(m,r)=X_est(3);
        end
    end
end
%结束滤波

%滤波结果的数据分析
err_x=zeros(N,1);
err_y=zeros(N,1);
delta_x=zeros(N,1);
delta_y=zeros(N,1);
%计算滤波的误差均值及标准差
for r=1:N
    %估计误差均值
    ex=sum(Rx(r)-Mt_Est_Px(:,r));
    ey=sum(Ry(r)-Mt_Est_Py(:,r));
    err_x(r)=ex/M;
    err_y(r)=ey/M;
    eqx=sum((Rx(r)-Mt_Est_Px(:,r)).^2);
    eqy=sum((Ry(r)-Mt_Est_Py(:,r)).^2);
    %估计误差标准差
    delta_x(r)=sqrt(abs(eqx/M-(err_x(r)^2)));
    delta_y(r)=sqrt(abs(eqy/M-(err_y(r)^2)));
end
%画图
figure(1);
plot(Rx,Ry,'k',Zx,Zy,'y:',xn,yn,'r');
legend('Actual Trace','Observe Trace','Estimate Trace');

figure(2);
subplot(2,1,1);
plot(err_x);
%axis([1,N,-300,300]);
title('X mean estimate error');
subplot(2,1,2);
plot(err_y);
%axis([1,N,-300,300]);
title('Y mean estimate error');

figure(3);
subplot(2,1,1);
plot(delta_x);
%axis([1,N,0,1]);
title('X estimate error standard deviation');
subplot(2,1,2);
plot(delta_y);
%axis([1,N,0,1]);
title('Y estimate error standard deviation');

% 子函数
% X_est,P_est返回第m次仿真第r个采样点的滤波结果
% Xn_est,Pn_est记录每个模型对应的第m次仿真第r次仿真第r个采样点的滤波结果
% u为模型概率
function [X_est,P_est,Xn_est,Pn_est,u]=IMM(Xn_est,Pn_est,T,Z,Delta,u)
% 控制模型转换的马尔科夫链的转移概率矩阵
P=[0.95,0.025,0.025;0.025,0.95,0.025;0.025,0.025,0.95];
%所采用的第三个模型参数，模型一为非机动，模型二、三均为机动模型
%模型一
PHI{1,1}=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];
PHI{1,1}(6,6)=0;
PHI{2,1}=[1,T,0,0,T^2/2,0;0,1,0,0,T,0;0,0,1,T,0,T^2/2;
    0,0,0,1,0,T;0,0,0,0,1,0;0,0,0,0,0,1];%模型模型二
PHI{3,1}=PHI{2,1};%模型三
G{1,1}=[T/2,0;1,0;0,T/2;0,1];%模型一
G{1,1}(6,2)=0;
G{2,1}=[T^2/4,0;T/2,0;0,T^2/4;0,T/2;1,0;0,1];%模型二
G{3,1}=G{2,1};%模型三
Q{1,1}=zeros(2);%模型一
Q{2,1}=0.005*eye(2);%模型二
Q{3,1}=0.015*eye(2);%模型三
H=[1,0,0,0,0,0;0,0,1,0,0,0];
R=eye(2)*Delta^2;%观测噪声协方差阵
mu=zeros(3,3);%混合概率矩阵μ
c_mean=zeros(1,3);%归一化常数
for i=1:3
    c_mean=c_mean+P(i,:)*u(i);
end
for i=1:3
    mu(i,:)=P(i,:)*u(i)./c_mean;
end
%输入交互
for j=1:3
    X0{j,1}=zeros(6,1);
    P0{j,1}=zeros(6);
    for i=1:3
        X0{j,1}=X0{j,1}+Xn_est{i,1}*mu(i,j);
    end
    for i=1:3
        P0{j,1}=P0{j,1}+mu(i,j)*( Pn_est{i,1}...
            +(Xn_est{i,1}-X0{j,1})*(Xn_est{i,1}-X0{j,1})');
    end
end
%模型条件滤波
a=zeros(1,3);
for j=1:3
    %观测预测
    X_pre{j,1}=PHI{j,1}*X0{j,1};
    %协方差预测
    P_pre{j,1}=PHI{j,1}*P0{j,1}*PHI{j,1}'+G{j,1}*Q{j,1}*G{j,1}';
    %计算卡尔曼增益
    K{j,1}=P_pre{j,1}*H'*inv(H*P_pre{j,1}*H'+R);
    %状态更新
    Xn_est{j,1}=X_pre{j,1}+K{j,1}*(Z-H*X_pre{j,1});
    %协方差更新
    Pn_est{j,1}=(eye(6)-K{j,1}*H)*P_pre{j,1};
end
%模型概率更新
for j=1:3
    v{j,1}=Z-H*X_pre{j,1};%新息
    s{j,1}=H*P_pre{j,1}*H'+R;%观测协方差矩阵
    n=length(s{j,1})/2;
    a(1,j)=1/((2*pi)^n*sqrt(det(s{j,1})))*exp(-0.5*v{j,1}'...
        *inv(s{j,1})*v{j,1});%观测对于模型j的似然函数
end
c=sum(a.*c_mean);%归一化常数
u=a.*c_mean./c;%模型概率更新
%输出交互
Xn=zeros(6,1);
Pn=zeros(6);
for j=1:3
    Xn=Xn+Xn_est{j,1}.*u(j);
end
for j=1:3
    Pn=Pn+u(j).*(Pn_est{j,1}+(Xn_est{j,1}-Xn)*(Xn_est{j,1}-Xn)');
end
%返回滤波结果
X_est=Xn;
P_est=Pn;