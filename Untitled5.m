%��ά�����˶�׷��

clear;clc;
%ʱ������ʼ��
T1 = 400;
T2 = 600;
T3 = 610;
T4 = 660;
T5 = 900;
T = 2;%�۲���
X_init = [1000 0 0 8000 -12]';%���ó�ʼλ�ã�λ��X��X���ٶȣ�X����ٶȣ�λ��Y��Y���ٶ�
b = 1000;
P = diag([b b b b b]);  %��ʼЭ����

meas_sigma = 1200; %������������
meas_noise = diag([meas_sigma   meas_sigma]);

pro_sigma = 40;%������������
pro_noise = diag([pro_sigma pro_sigma/T pro_sigma/T.^2 pro_sigma pro_sigma/T]);

%״̬ת�ƾ���
A = [1 T 0.5*T*T 0 0;
     0 1    T    0 0;
     0 0    1    0 0;
     0 0    0    1 T;
     0 0    0    0 1];
 %�۲ⷽ��
 H =[1 0 0 0 0;
     0 0 0 1 0];

 X_true_trac = X_init;%��ʵ�켣
 Z_meas_trac = [1000;8000];%�����Ĺ켣
 X_est_trac = [1000;8000];%���ƵĹ켣
 X_init_copy = X_init;%copyΪ�м���������ڸ�ֵ
 
 %�������˲�
for i = 1:T5/2        %�жϼ��ٶȲ���ֵ
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

    %������ʵ�Ĺ켣
    X_true     = A * X_init(:,i);
    X_init(:,i+1)     = X_true;  %����
    X_true_trac(:,i+1) = X_true; %��ͼ�õ�

    Z_meas_trac(:,i+1) = [X_true(1)+wgn(1,1,30);X_true(4)+wgn(1,1,30)];%����ʵֵ������1X1�ĸ�˹����

    %�������˲�
    %һ��Ԥ��
    X_est =  A * X_init_copy(:,i);
    P_est =  A * P * A' + pro_noise;

    %�˲�����
    K = P_est * H' * inv(H *P_est* H' + meas_noise);
    X_pre = X_est + K * ( Z_meas_trac(:,i+1) - H * X_est);
    P_pre = (eye(5) - K * H) * P_est;

    X_init_copy(:,i+1) = X_pre;
    P = P_pre;

    X_est_trac(:,i+1) = [X_pre(1); X_pre(4)];%�������˲������ֵ

end
plot(X_true_trac(1,:),X_true_trac(4,:),'r',Z_meas_trac(1,:),Z_meas_trac(2,:),'y',X_est_trac(1,:),X_est_trac(2,:),'g')
title('trace');
xlabel('X(m)');
ylabel('Y(m)');
legend('Actual Trace','Observe Trace','Estimate Trace');