function self_test
clc;clear;
T=1;
N=200/T;
X=zeros(4,N);
X(:,1)=[100,5,800,-5];
Z=zeros(2,N);
Z(:,1)=[X(1,1),X(3,1)];
S(:,1)=[100,5,800,-5];
delta_w=0.01;
Q=delta_w*diag([1,1,1,1]);
R=eye(2);

A=[1,T,0,0;
   0,1,0,0;
   0,0,1,T;
   0,0,0,1];
H=[1,0,0,0;
   0,0,1,0];

for i=1:N
    S(:,i+1)=A*S(:,i);
    X(:,i+1)=A*X(:,i)+sqrtm(Q)*randn(4,1);
    Z(:,i+1)=H*X(:,i+1)+sqrtm(R)*randn(2,1);
end

X_pre=zeros(4,N);
X_pre(:,1)=X(:,1);
P_pre=eye(4);
X_est=zeros(4,0);

for i=1:N
    X_est=A*X_pre(:,i);
    P_est=A*P_pre*A'+Q;
    K=P_est*H'*inv(H*P_est*H'+R);
    X_pre(:,i+1)=X_est+K*(Z(:,i+1)-H*X_est); 
    P_pre=(eye(4)-K*H)*P_est;
end

for i=1:N
    Err_observation(i)=DIS(X(:,i),Z(:,i));
    Err_kalman(i)=DIS(X(:,i),X_pre(:,i));
end

function dis=DIS(X1,X2)
    if length(X2)<=2
        dis=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
    else
        dis=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2);
    end 
end

figure
hold on;box on;
plot(Err_observation);
plot(Err_kalman);
legend("observ-err","kal-err");
xlabel('Time');
ylabel('error');

figure
hold on;box on;
plot(S(1,:),S(3,:),'g','LineWidth',1);
plot(X(1,:),X(3,:),'k','LineWidth',1);
plot(Z(1,:),Z(2,:),'r','LineWidth',1);
plot(X_pre(1,:),X_pre(3,:),'c','LineWidth',1);
legend('theory trace','real trace','observ trace','Kalman trace');
xlabel('X distance');
ylabel('Y distance');

end
