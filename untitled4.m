function untitled4
%�����˶�
clc;clear;
T=1;                                %�״�ɨ������
N=100/T;                            %�ܵĲ�������
X=zeros(4,N);                       %Ŀ����ʵλ�á��ٶ�
X(:,1)=[-100,2,200,20];             %Ŀ���ʼλ�á��ٶ�,X��λ�ã�X���ٶȣ�Y��λ�ã�Y���ٶ�
S(:,1)=[-100,2,200,20];             %Ŀ�����۹켣��ʼֵ
Z=zeros(2,N);                       %��������λ�õĹ۲�
Z(:,1)=[X(1,1),X(3,1)];             %�۲��ʼ��
delta_w=1e-2;                       %����������������Ŀ����ʵ�켣���������ˣ�Ϊ0.01
Q=delta_w*diag([0.5,1,0.5,1]);      %����������ֵ����
R=eye(2);                           %�۲�������ֵ����
A=[1,T,0,0;
   0,1,0,0;
   0,0,1,T;
   0,0,0,1];                        %״̬ת�ƾ���
H=[1,0,0,0;
   0,0,1,0];                        %�۲����
X_est=zeros(4,0);                   %�������ֵ
for i=2:N
    S(:,i)=A*S(:,i-1);%Ŀ�����۹켣
    X(:,i)=A*X(:,i-1)+sqrtm(Q)*randn(4,1);%Ŀ����ʵ�켣��sqrtmΪ����ƽ����
    Z(:,i)=H*X(:,i)+sqrtm(R)*randn(2,1);%��Ŀ��Ĺ۲�
end

% Kalman �˲�
X_pre=zeros(4,N);
X_pre(:,1)=X(:,1);%�������˲�״̬��ʼ��
%M(1,:)=X_pre(:,1);
P_pre=100e-2*eye(4);% Э�������ʼ��
 

for i=2:N
    X_est=A*X_pre(:,i-1);                %����Э����
    %M(i,:)=X_est;
    P_est=A*P_pre*A'+Q;                   %Ԥ�����Э����
    K=P_est*H'*inv(H*P_est*H'+R);         %���������棬invΪ����
    X_pre(:,i)=X_est+K*(Z(:,i)-H*X_est);  %״̬����
    P_pre=(eye(4)-K*H)*P_est;             %�˲�Ԥ�����Э�������
end
% ������
for i=1:N

    Err_Observation(i)=RMS(X(:,i),Z(:,i));%�۲�ֵ����ʵֵ�����
    Err_KalmanFilter(i)=RMS(X(:,i),X_pre(:,i));%����ֵ����ʵֵ�����
end


 figure
 hold on;box on;
 plot(S(1,:),S(3,:),'g','LineWidth',1);%���۹켣
 plot(X(1,:),X(3,:),'b','LineWidth',1);%��ʵ�켣
 plot(Z(1,:),Z(2,:),'r','LineWidth',1);%�۲�켣
 plot(X_pre(1,:),X_pre(3,:),'c','LineWidth',1);%�������˲��켣
 legend('���۹켣','��ʵ�켣','�۲�켣','�˲���켣');
 xlabel('������ X/m');
 ylabel('������ Y/m');
 
figure
hold on;box on;
plot(Err_Observation);
plot(Err_KalmanFilter);
legend('�˲�ǰ���','�˲������');
xlabel('�۲�ʱ��/s');
ylabel('���ֵ');

% ����ŷ�Ͼ����Ӻ���
function dist=RMS(X1,X2)
if length(X2)<=2
    dist=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
else
    dist=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2);
end
