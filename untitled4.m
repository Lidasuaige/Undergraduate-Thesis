function untitled4
%匀速运动
clc;clear;
T=1;                                %雷达扫描周期
N=100/T;                            %总的采样次数
X=zeros(4,N);                       %目标真实位置、速度
X(:,1)=[-100,2,200,20];             %目标初始位置、速度,X轴位置，X轴速度，Y轴位置，Y轴速度
S(:,1)=[-100,2,200,20];             %目标理论轨迹初始值
Z=zeros(2,N);                       %传感器对位置的观测
Z(:,1)=[X(1,1),X(3,1)];             %观测初始化
delta_w=1e-2;                       %如果增大这个参数，目标真实轨迹就是曲线了，为0.01
Q=delta_w*diag([0.5,1,0.5,1]);      %过程噪声均值矩阵
R=eye(2);                           %观测噪声均值矩阵
A=[1,T,0,0;
   0,1,0,0;
   0,0,1,T;
   0,0,0,1];                        %状态转移矩阵
H=[1,0,0,0;
   0,0,1,0];                        %观测矩阵
X_est=zeros(4,0);                   %定义估计值
for i=2:N
    S(:,i)=A*S(:,i-1);%目标理论轨迹
    X(:,i)=A*X(:,i-1)+sqrtm(Q)*randn(4,1);%目标真实轨迹，sqrtm为矩阵平方根
    Z(:,i)=H*X(:,i)+sqrtm(R)*randn(2,1);%对目标的观测
end

% Kalman 滤波
X_pre=zeros(4,N);
X_pre(:,1)=X(:,1);%卡尔曼滤波状态初始化
%M(1,:)=X_pre(:,1);
P_pre=100e-2*eye(4);% 协方差阵初始化
 

for i=2:N
    X_est=A*X_pre(:,i-1);                %估计协方差
    %M(i,:)=X_est;
    P_est=A*P_pre*A'+Q;                   %预测误差协方差
    K=P_est*H'*inv(H*P_est*H'+R);         %卡尔曼增益，inv为求逆
    X_pre(:,i)=X_est+K*(Z(:,i)-H*X_est);  %状态更新
    P_pre=(eye(4)-K*H)*P_est;             %滤波预测误差协方差更新
end
% 误差分析
for i=1:N

    Err_Observation(i)=RMS(X(:,i),Z(:,i));%观测值和真实值的误差
    Err_KalmanFilter(i)=RMS(X(:,i),X_pre(:,i));%估计值和真实值的误差
end


 figure
 hold on;box on;
 plot(S(1,:),S(3,:),'g','LineWidth',1);%理论轨迹
 plot(X(1,:),X(3,:),'b','LineWidth',1);%真实轨迹
 plot(Z(1,:),Z(2,:),'r','LineWidth',1);%观测轨迹
 plot(X_pre(1,:),X_pre(3,:),'c','LineWidth',1);%卡尔曼滤波轨迹
 legend('理论轨迹','真实轨迹','观测轨迹','滤波后轨迹');
 xlabel('横坐标 X/m');
 ylabel('纵坐标 Y/m');
 
figure
hold on;box on;
plot(Err_Observation);
plot(Err_KalmanFilter);
legend('滤波前误差','滤波后误差');
xlabel('观测时间/s');
ylabel('误差值');

% 计算欧氏距离子函数
function dist=RMS(X1,X2)
if length(X2)<=2
    dist=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
else
    dist=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2);
end
