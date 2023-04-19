%二维蛇形运动追踪

clear;clc;
%时间间隔初始化
T1 = 400;
T2 = 600;
T3 = 610;
T4 = 660;
T5 = 900;
T = 2;%观测间隔
X_init = [1000 0 0 8000 -12]';%设置初始位置，位置X，X轴速度，X轴加速度，位置Y，Y轴速度
b = 1000;
P = diag([b b b b b]);  %初始协方差

meas_sigma = 1200; %测量噪声方差
meas_noise = diag([meas_sigma   meas_sigma]);

pro_sigma = 40;%过程噪声方差
pro_noise = diag([pro_sigma pro_sigma/T pro_sigma/T.^2 pro_sigma pro_sigma/T]);

%状态转移矩阵
A = [1 T 0.5*T*T 0 0;
     0 1    T    0 0;
     0 0    1    0 0;
     0 0    0    1 T;
     0 0    0    0 1];
 %观测方程
 H =[1 0 0 0 0;
     0 0 0 1 0];

 X_true_trac = X_init;%真实轨迹
 Z_meas_trac = [1000;8000];%测量的轨迹
 X_est_trac = [1000;8000];%估计的轨迹
 X_init_copy = X_init;%copy为中间变量，便于赋值
 
 %卡尔曼滤波
for i = 1:T5/2        %判断加速度并赋值
    if i>=T1/2 &&i<=T2/2
        a      = 0.075;
        X_init(3,i) = a;
        X_init_copy(3,i) =a;
    end
    if i>=T2/2 &&i<=T3/2
        a      = 0;
        X_init(3,i) = a;
        X_init_copy(3,i) =a;

    end  
     if i>=T3/2 &&i<T4/2
        a      = -0.3;
        X_init(3,i) = a;
        X_init_copy(3,i) =a;

     end  
    if i>=T4/2
        a           = 0;
        X_init(3,i) = a;  
        X_init_copy(3,i) =a;

    end

    %产生真实的轨迹
    X_true     = A * X_init(:,i);
    X_init(:,i+1)     = X_true;  %迭代
    X_true_trac(:,i+1) = X_true; %画图用的

    Z_meas_trac(:,i+1) = [X_true(1)+wgn(1,1,30);X_true(4)+wgn(1,1,30)];%在真实值上增加1X1的高斯噪声

    %卡尔曼滤波
    %一步预测
    X_est =  A * X_init_copy(:,i);
    P_est =  A * P * A' + pro_noise;

    %滤波更新
    K = P_est * H' * inv(H *P_est* H' + meas_noise);
    X_pre = X_est + K * ( Z_meas_trac(:,i+1) - H * X_est);
    P_pre = (eye(5) - K * H) * P_est;

    X_init_copy(:,i+1) = X_pre;
    P = P_pre;

    X_est_trac(:,i+1) = [X_pre(1); X_pre(4)];%卡尔曼滤波后估计值

end
plot(X_true_trac(1,:),X_true_trac(4,:),'r',Z_meas_trac(1,:),Z_meas_trac(2,:),'y',X_est_trac(1,:),X_est_trac(2,:),'g')
title('trace');
xlabel('X(m)');
ylabel('Y(m)');
legend('Actual Trace','Observe Trace','Estimate Trace');