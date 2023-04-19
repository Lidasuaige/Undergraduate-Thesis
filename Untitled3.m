Z=(1:100);    %第一秒观测1m，第二秒观测2m，匀速运动
noise=randn(1,100);  %方差为1的高斯噪声
Z=Z+noise; 

X=[0;0];  %初始状态
P=[1 0;0 1];
F=[1 1;0 1];
Q=[0.0001, 0;0 0.0001];
H=[1 0];
R=1;       %观测噪声方差为1

figure;
hold on;

for i=1:100
    X_=F*X;
    P_=F*P*F'+Q;
    K=P_*H'/(H*P_*H'+R);
    X=X_+K*(Z(i)-H*X_);
    P=(eye(2)-K*H)*P_;
    
    plot(X(1),X(2),'b.');
end