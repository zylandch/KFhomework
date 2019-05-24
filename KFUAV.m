%���������ĵع�ת����γ��
clc
clear 
lamda=116;%��λ��
phi=40;
h=100;

Ts=1;%�������
T=100;%����ʱ��
N=T/Ts;
V=5;%�ٶ�

dtor=pi/180;
sigma1=[3;3;4.5];%��γ��
a=6378137;%Re
f=1/298.257223;%����
e=sqrt(f*(2-f));%ƫ����
Rn=a/sqrt(1-e^2*sin(phi*dtor)^2);%î�ϰ뾶
Rm=a*(1+f*sin(phi*dtor)^2);%����뾶

K1=1/(2*pi*(Rn+h)*cos(phi*dtor)/360);%����ת��
K2=1/(2*pi*Rn/360);
K3=1;

sigma2=sigma1.*[K1;K2;K3];%λ�ö�λ��׼��
sigma3=[0.1;0.1;0.15];%�ٶȶ�λ��׼��

Ve=0;
Vn=V*cos(30*dtor);
Vu=V*sin(30*dtor);

Ytrue=[116;40;100;Ve;Vn;Vu];
Y1=[116;40;100]+sigma2.*randn(3,1);
Y2=[Ve;Vn;Vu]+sigma3.*randn(3,1);%�ٶȲ���
Ysum=[Y1;Y2]; %���ٶȲ���

x0=[lamda;phi;h;Ve;Vn;Vu];%��ʼ״̬
x1=x0;
X=x0;
X2=x0;
P1=diag([3 3 4.5 2 2 2].*[K1^2 K2^2 K3^2 1 1 1]);%��ʼ����

 Rnh=1/((Rn+h)*cos(phi*dtor));
 Rmh=1/(Rm+h);
F=[0 0 0 Rnh 0 0;
0 0 0 0 Rmh 0;
0 0 0 0 0 1;
0 0 0 0 0 0;
0 0 0 0 0 0;
0 0 0 0 0 0];%ת̬ת�ƾ���
Phi1=eye(6,6)+F.*Ts;%��ɢ��
G=zeros(6,6);
H=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0];%�۲����
Q=diag([0 0 0 0.1 0.1 0.1])./Ts;%Э������Q
G=eye(6,6);
Qg=G*Q*G';
t=[0];

%ʵ��λ���ٶ�����ջ�����������
for i=1:N
  Rn=a/sqrt(1-e^2*sin(phi*dtor)^2);
  Rm=a*(1+f*sin(phi*dtor)^2);
  Rnh=1/((Rn+h)*cos(phi*dtor));
  Rmh=1/(Rm+h);
  Ve=0+sqrt(0.1)*randn(1,1);%ʵ���ٶȿ��ǰ���������[0.1 0.1 0.1]
Vn=V*cos(30*dtor)+sqrt(0.1)*randn(1,1);
Vu=V*sin(30*dtor)+sqrt(0.1)*randn(1,1);

  Y1=Y1+[Ts*Ve*Rnh/dtor;Ts*Vn*Rmh/dtor;Vu*Ts];
  Y2=[Ve;Vn;Vu];
  Y=[Y1;Y2];
  Ytrue=[Ytrue Y];%ʵ��״̬
  
    Y3=Y1+sigma2.*randn(3,1);
  Y4=Y2+sigma3.*randn(3,1);
  Yerro=[Y3;Y4];
   Ysum=[Ysum Yerro];%���ջ�����״̬
   
     t=[t Ts*i];%ʵʱʱ��
     
   phi=Y(2);
    h=Y(3);
end

lamda=116;
phi=40;
h=100;
%KF���ջ��ɲ�λ�ú��ٶ�
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
Phi1=eye(6,6)+F.*Ts;%��ɢ��
QD=Qg*Ts+(F*Qg+Qg*F')*(Ts^2/2)+2*F*Qg*F'*(Ts^3/6);%��ɢ��

K1=1/(2*pi*(Rn+h)*cos(phi*dtor)/360);
K2=1/(2*pi*Rn/360);
K3=1;

H=eye(6,6);
R=diag([3^2 3^2 4.5^2 0.1^2 0.1^2 0.15^2].*[K1^2 K2^2 K3^2 1 1 1])./Ts;%���ⷽ��

%�˲�
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
lamda=116;%��λ��
phi=40;
h=100;

%KF���ջ�ֻ��λ��
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
  xlabel('����');
ylabel('γ��');
zlabel('�߶�');
 title('λ��--ʱ��');
 legend('��ʵλ��','��λ�ú��ٶ�','ֻ��λ��');
 
  figure(2)
 plot(t(1,:),Ytrue(4,:));
     hold on
 plot(t(1,:),X(4,:),'r');
 hold on
  plot(t(1,:),X2(4,:),'y');
   xlabel('ʱ��/s');
ylabel('�����ٶ�');
 title('�����ٶ�--ʱ��');
 legend('��ʵ�ٶ�','��λ�ú��ٶ�','ֻ��λ��'); 
   
  figure(3)
 plot(t(1,:),Ytrue(5,:));
     hold on
 plot(t(1,:),X(5,:),'r');
 hold on
  plot(t(1,:),X2(5,:),'y');
      xlabel('ʱ��/s');
ylabel('�����ٶ�');
 title('�����ٶ�--ʱ��');
 legend('��ʵ�ٶ�','��λ�ú��ٶ�','ֻ��λ��'); 
 
  figure(4)
 plot(t(1,:),Ytrue(6,:));
     hold on
 plot(t(1,:),X(6,:),'r');
 hold on
  plot(t(1,:),X2(6,:),'y');
     xlabel('ʱ��/s');
ylabel('�����ٶ�');
 title('�����ٶ�--ʱ��');
 legend('��ʵ�ٶ�','��λ�ú��ٶ�','ֻ��λ��'); 
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
