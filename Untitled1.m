clear    
N=200;                      %采样次数为200次
w(1)=0;                     %w为过程噪声  
w=randn(1,N);   
x(1)=25;                    %房间初始温度为25°
a=1;                        %a为方程中A(k)，呈默认线性关系
for k=2:N    
x(k)=a*x(k-1)+w(k-1);    
end    
V=randn(1,N);               %V为观察噪声  
q1=std(V);                  %q1为V的误差标准差
Rvv=q1.^2;                  %Rvv为V的误差协方差                    
q3=std(w);    
Rww=q3.^2;                  %Rww为W的误差协方差   
c=0.2;                      %c为方程中H(k)  
Y=c*x+V;                    %Y为观察值  
p(1)=0;                     %初始估计误差协方差
s(1)=0;                     %s
for t=2:N    
p1(t)=a.^2*p(t-1)+Rww;     %p1为方程中p'，预测协方差  
b(t)=c*p1(t)/(c.^2*p1(t)+Rvv);  %b为卡尔曼增益  
s(t)=a*s(t-1)+b(t)*(Y(t)-a*c*s(t-1));    %s为卡尔曼估计值
p(t)=p1(t)-c*b(t)*p1(t);    %p为估计协方差
end
t=1:N;    
plot(t,s,'r',t,Y,'g',t,x,'b'); %卡尔曼估计值为红色，观察值为绿色，真实值为蓝色
legend('卡尔曼估计值','观察值','真实值');