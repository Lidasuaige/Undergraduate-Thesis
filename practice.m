%作业练习

% z=[2.347;
%    3.220;
%    10.026;
%    12.965];
% 
% H=[1,sin(1),sin(3);
%    4,sin(2),sin(6);
%    9,sin(3),sin(9);
%    16,sin(4),sin(12)];
% 
% 
% x_est1=inv(H'*H)*H'*z;

clear;
T=1;
N=100;
H=[1;1];
sigma_1=0.1;
sigma_2=0.2;
v1=sigma_1*randn(1,100);
v2=sigma_2*randn(1,100);
x_true=30*ones(1,100);
x_est=zeros(1,100);
x_est(1)=30;

z1=x_true+v1;
z2=x_true+v2;
z=[z1;z2];
error_1=zeros(1,100);
error_2=zeros(1,100);

v=[v1;v2];

R=[0.01 0;
    0 0.04];
k_g=zeros(1,2*N);

p=zeros(1,100);
p(1)=inv(H'*inv(R)*H);

for k=1:N-1

    k_g=p(k)*H'*inv(H*p(k)*H'+R);
    x_est(k+1)=x_est(k)+k_g*(z(:,k+1)-H*x_est(k));
    p(k+1)=inv(inv(p(k))+H'*inv(R)*H);
    error_1=x_true-z1;
    error_2=x_true-x_est;
    
end




figure(1);
plot(1:N,x_true,'b',1:N,x_est,'r',1:N,z1,'y',1:N,z2,'g');
legend('true','estimation','observ1','observ2');
title('Kalman course project');

figure(2);
plot(1:N,error_1,'b',1:N,error_2,'r');
legend('error-observ','error-est');
title('error different');
    
t=0:0.1:2*pi;
y1=sin(t);
y2=cos(t);
figure;
plot(t, y1, t, y2);
xlabel("zibianliang");
ylabel("yinbianliang");
legend("y1","y2");

A=[0 3 3;-1 8 6;2 -14 -10];
[V,D]=eig(A);

A=input('input A here: ');
if A>10
    disp(A)
end


s=0;
num=input('input the number you want to sum: ');
for i=1:1:num
    s=s+i;
end
disp(s);



