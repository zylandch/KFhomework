%迭代将地心地固转到经纬高
clc
clear 
lamda=116;%单位度
phi=40;
h=100;

Ts=1;%采样间隔
T=100;%采样时间
N=T/Ts;
V=5;%速度

dtor=pi/180;
sigma1=[3;3;4.5];%经纬高
a=6378137;%Re
f=1/298.257223;%扁率
e=sqrt(f*(2-f));%偏心率
Rn=a/sqrt(1-e^2*sin(phi*dtor)^2);%卯酉半径
Rm=a*(1+f*sin(phi*dtor)^2);%子午半径

K1=1/(2*pi*(Rn+h)*cos(phi*dtor)/360);%距离转度
K2=1/(2*pi*Rn/360);
K3=1;

sigma2=sigma1.*[K1;K2;K3];%位置定位标准差
sigma3=[0.1;0.1;0.15];%速度定位标准差

Ve=0;
Vn=V*cos(30*dtor);
Vu=V*sin(30*dtor);

Ytrue=[116;40;100;Ve;Vn;Vu];
Y1=[116;40;100]+sigma2.*randn(3,1);
Y2=[Ve;Vn;Vu]+sigma3.*randn(3,1);%速度测量
Ysum=[Y1;Y2]; %有速度测量

x0=[lamda;phi;h;Ve;Vn;Vu];%初始状态
x1=x0;
X=x0;
X2=x0;
P1=diag([3 3 4.5 2 2 2].*[K1^2 K2^2 K3^2 1 1 1]);%初始误差方差

 Rnh=1/((Rn+h)*cos(phi*dtor));
 Rmh=1/(Rm+h);
F=[0 0 0 Rnh 0 0;
0 0 0 0 Rmh 0;
0 0 0 0 0 1;
0 0 0 0 0 0;
0 0 0 0 0 0;
0 0 0 0 0 0];%转态转移矩阵
Phi1=eye(6,6)+F.*Ts;%离散化
G=zeros(6,6);
H=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0];%观测矩阵
Q=diag([0 0 0 0.1 0.1 0.1])./Ts;%协方差阵Q
G=eye(6,6);
Qg=G*Q*G';
t=[0];

%实际位置速度与接收机测量量生成
for i=1:N
  Rn=a/sqrt(1-e^2*sin(phi*dtor)^2);
  Rm=a*(1+f*sin(phi*dtor)^2);
  Rnh=1/((Rn+h)*cos(phi*dtor));
  Rmh=1/(Rm+h);
  Ve=0+sqrt(0.1)*randn(1,1);%实际速度考虑白噪声方差[0.1 0.1 0.1]
Vn=V*cos(30*dtor)+sqrt(0.1)*randn(1,1);
Vu=V*sin(30*dtor)+sqrt(0.1)*randn(1,1);

  Y1=Y1+[Ts*Ve*Rnh/dtor;Ts*Vn*Rmh/dtor;Vu*Ts];
  Y2=[Ve;Vn;Vu];
  Y=[Y1;Y2];
  Ytrue=[Ytrue Y];%实际状态
  
    Y3=Y1+sigma2.*randn(3,1);
  Y4=Y2+sigma3.*randn(3,1);
  Yerro=[Y3;Y4];
   Ysum=[Ysum Yerro];%接收机测量状态
   
     t=[t Ts*i];%实时时间
     
   phi=Y(2);
    h=Y(3);
end

lamda=116;
phi=40;
h=100;
%KF接收机可测位置和速度
for i=1:N
    Rn=a/sqrt(1-e^2*sin(phi*dtor)^2);
    Rm=a*(1+f*sin(phi*dtor)^2);
    Rnh=1/((Rn+h)*cos(phi*dtor));
    Rmh=1/(Rm+h);
    F=[0 0 0 Rnh/dtor 0 0;
0 0 0 0 Rmh/dtor 0;
0 0 0 0 0 1;
0 0 0 0 0 0;
0 0 0 0 0 0;
0 0 0 0 0 0];
Phi1=eye(6,6)+F.*Ts;%离散化
QD=Qg*Ts+(F*Qg+Qg*F')*(Ts^2/2)+2*F*Qg*F'*(Ts^3/6);%离散化

K1=1/(2*pi*(Rn+h)*cos(phi*dtor)/360);
K2=1/(2*pi*Rn/360);
K3=1;

H=eye(6,6);
R=diag([3^2 3^2 4.5^2 0.1^2 0.1^2 0.15^2].*[K1^2 K2^2 K3^2 1 1 1])./Ts;%量测方差

%滤波
    x12=Phi1*x1;
    y=H*x12;
    P12=Phi1*P1*Phi1'+QD;
    P2=inv(inv(P12)+H'*inv(R)*H);
    km=P2*H'*inv(R);
    x1=x12+km*(Ysum(:,i+1)-y);
    X=[X x1];
    
     phi=x1(2);
    h=x1(3);
end

x1=x0;
lamda=116;%单位度
phi=40;
h=100;

%KF接收机只测位置
for i=1:N
    Rn=a/sqrt(1-e^2*sin(phi*dtor)^2);
    Rm=a*(1+f*sin(phi*dtor)^2);
    Rnh=1/((Rn+h)*cos(phi*dtor));
    Rmh=1/(Rm+h);
    F=[0 0 0 Rnh/dtor 0 0;
0 0 0 0 Rmh/dtor 0;
0 0 0 0 0 1;
0 0 0 0 0 0;
0 0 0 0 0 0;
0 0 0 0 0 0];
Phi1=eye(6,6)+F.*Ts;
QD=Qg*Ts+(F*Qg+Qg*F')*(Ts^2/2)+2*F*Qg*F'*(Ts^3/6);

K1=1/(2*pi*(Rn+h)*cos(phi*dtor)/360);
K2=1/(2*pi*Rn/360);
K3=1;

H=eye(3,6);
R=diag([3^2 3^2 4.5^2].*[K1^2 K2^2 K3^2]);

    x12=Phi1*x1;
    y=H*x12;
    P12=Phi1*P1*Phi1'+QD;
    P2=inv(inv(P12)+H'*inv(R)*H);
    km=P2*H'*inv(R);
    x1=x12+km*(Ysum(1:3,i+1)-y);
    X2=[X2 x1];
    
     phi=x1(2);
    h=x1(3);
end

figure(1)
 plot3(Ytrue(1,:),Ytrue(2,:),Ytrue(3,:));
     hold on
 plot3(X(1,:),X(2,:),X(3,:),'r');
 hold on
  plot3(X2(1,:),X2(2,:),X2(3,:),'y');
  xlabel('经度');
ylabel('纬度');
zlabel('高度');
 title('位置--时间');
 legend('真实位置','测位置和速度','只测位置');
 
  figure(2)
 plot(t(1,:),Ytrue(4,:));
     hold on
 plot(t(1,:),X(4,:),'r');
 hold on
  plot(t(1,:),X2(4,:),'y');
   xlabel('时间/s');
ylabel('东向速度');
 title('东向速度--时间');
 legend('真实速度','测位置和速度','只测位置'); 
   
  figure(3)
 plot(t(1,:),Ytrue(5,:));
     hold on
 plot(t(1,:),X(5,:),'r');
 hold on
  plot(t(1,:),X2(5,:),'y');
      xlabel('时间/s');
ylabel('北向速度');
 title('北向速度--时间');
 legend('真实速度','测位置和速度','只测位置'); 
 
  figure(4)
 plot(t(1,:),Ytrue(6,:));
     hold on
 plot(t(1,:),X(6,:),'r');
 hold on
  plot(t(1,:),X2(6,:),'y');
     xlabel('时间/s');
ylabel('天向速度');
 title('天向速度--时间');
 legend('真实速度','测位置和速度','只测位置'); 
%  figure(2)
%  plot(t(1,:),Ytrue(3,:));
%  hold on
%  plot(t(1,:),X(3,:),'r');
%  hold on
%  plot(t(1,:),X2(3,:),'y');
% 
%   figure(3)
%  plot(t(1,:),Ytrue(1,:));
%  hold on
%  plot(t(1,:),X(1,:),'r');
%  hold on
%  plot(t(1,:),X2(1,:),'y');
%  
%   figure(4)
%  plot(t(1,:),Ytrue(2,:));
%  hold on
%  plot(t(1,:),X(2,:),'r');
%  hold on
%  plot(t(1,:),X2(2,:),'y');
%  
%  figure(5)
%  Xsum=sqrt(X(4,:).*X(4,:)+X(5,:).*X(5,:)+X(6,:).*X(6,:));
%   plot(Xsum(1,:),'*r');
% 
