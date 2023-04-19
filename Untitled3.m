Z=(1:100);    %��һ��۲�1m���ڶ���۲�2m�������˶�
noise=randn(1,100);  %����Ϊ1�ĸ�˹����
Z=Z+noise; 

X=[0;0];  %��ʼ״̬
P=[1 0;0 1];
F=[1 1;0 1];
Q=[0.0001, 0;0 0.0001];
H=[1 0];
R=1;       %�۲���������Ϊ1

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