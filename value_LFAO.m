%looping for 3D plot
% Added Fragmentation

function value_LFAO
clear all
clc
n=48; % gateway 10

% P=60000000;
% x=1.65e5; % gateway 156 
% x1 =1.64e5;
% y=18.4e-1;
% y1=2.0e-1; 

% x=1.65e4; % gateway 144 30000000
% x1 =1.62e4;
% y=16.2e-1;
% y1=1.9e-1; 

% P=26000000;
% x=1.65e3; % gateway 132 
% x1 =1.52e3;
% y=14.0e-1;
% y1=1.7e-1; 

% P=26000000;                                                                                  
% x=180000e-3; % gateway 120 
% x1 =153000e-3;
% y=11.8e-1;
% y1=1.6e-1; 

% P=17000000;
% x=25000e-3; % gateway 108 
% x1 =20000e-3;
% y=8.8e-1;
% y1=1.1e-1; 

% P=17000000;
% x=2600e-3; % gateway 96  
% x1 =2000e-3;
% y=7.0e-1;
% y1=1.0e-1; 

% P=13000000;
% x=500e-3; % gateway 84  
% x1 =430e-3;
% y=5.2e-1;
% y1=0.8e-1; 

% P=6000000;
% x=130e-3; % gateway 72 
% x1 =120e-3;
% y=3.52e-1;
% y1=0.65e-1; 

% P=6000000;
% x=100e-3; % gateway 60  6000000
% x1 =40e-3;
% y=2.10e-1;
% y1=0.70e-1; 

P=120000;
x=100e-3; % gateway 60  6000000
x1 =40e-3;
y=0.40e-1;
y1=0.10e-1; 

z=3e7;
z1=5e-3;
p=5e5;
p1=5e1;
r=5e4;

u=4.5e-3;
u1 =1e-4;
v=4e7;
v1=1e-3;


A_1=0.5;
B_12=[0.0068 0.002 0.0005];
C_24=[0.0466 0.004 .00025 ];


for j=1:3
    A_12=B_12(j);
    D_24=C_24(j);

theta=[x,x1,y,y1,z,z1,p,p1,r,u,u1,v,v1]; 
Y0=zeros(1,n); 


Y0(n-9)=A_1;
Y0(1)=A_12;
Y0(13)=D_24;

t_range=linspace(0,337,337); 
[t_val,Y_val]=ode23s(@ode_LFAO_1,t_range,Y0,[],n,theta);
Y_val([1:20:300 ],[1 11 13 25 n-23 n-22  n-10 n-9 n])

signalON=Y_val(:,n)*0;


for i=14:n-23
signalON=signalON + Y_val(:,i)*(i-13);
end

signalON(1:20:100)

for i=n-22:n-10
signalON=signalON + Y_val(:,i)*P;
end

signalON=signalON +Y_val(:,n)*P;

signalON(1:50:end)
 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));

if j==1
plot(t_range, signalON,'g')
hold on;

% B= [t_range',signalON];
% fileID = fopen('model2_10_Simulated.txt','w');
% fprintf(fileID,'%6.2f %12.8f\n',B');
% fclose(fileID);


% A=[t_range',Y_val(:,1),Y_val(:,13),Y_val(:,n-23),Y_val(:,n-22),Y_val(:,n-10),Y_val(:,n-9),Y_val(:,end)];
% fileID = fopen('oligomer_10.txt','w');
% fprintf(fileID,'%12.2f %12.12f %12.12f %12.12f %12.12f %12.12f %12.12f %12.12f\n',A');
% fclose(fileID);

end

if j==2
plot(t_range, signalON,'b')
hold on;
% B= [t_range',signalON];
% fileID = fopen('model2_01_Simulated.txt','w');
% fprintf(fileID,'%6.2f %12.8f\n',B');
% fclose(fileID);

% A=[t_range',Y_val(:,1),Y_val(:,13),Y_val(:,n-23),Y_val(:,n-22),Y_val(:,n-10),Y_val(:,n-9),Y_val(:,end)];
% fileID = fopen('oligomer_01.txt','w');
% fprintf(fileID,'%12.2f %12.12f %12.12f %12.12f %12.12f %12.12f %12.12f %12.12f\n',A');
% fclose(fileID);

end

if j==3
plot(t_range, signalON,'r')
hold on;

% B= [t_range',signalON];
% fileID = fopen('model2_00001_Simulated.txt','w');
% fprintf(fileID,'%6.2f %12.8f\n',B');
% fclose(fileID);

% A=[t_range',Y_val(:,1),Y_val(:,13),Y_val(:,n-23),Y_val(:,n-22),Y_val(:,n-10),Y_val(:,n-9),Y_val(:,end)];
% fileID = fopen('oligomer_0001.txt','w');
% fprintf(fileID,'%12.2f %12.12f %12.12f %12.12f %12.12f %12.12f %12.12f %12.12f\n',A');
% fclose(fileID);

end

if (j==1)
     load 'LFAO_DATA.txt';
    Data=LFAO_DATA;
    plot(Data(:,1),Data(:,2),'-*')
    hold on
elseif (j==2)
    load 'LFAO_DATA_01.txt';
    Data=LFAO_DATA_01;
    plot(Data(:,1),Data(:,2),'-*')
    hold on
 else
    load 'LFAO_DATA_00001.txt';
    Data=LFAO_DATA_00001;
    plot(Data(:,1),Data(:,2),'-*')
    hold on
 end

X=Data(:,2);
Y=signalON(Data(:,1)+1);
mdl = fitlm(Y,X)

signalON (192)/signalON (144)
signalON (300)/signalON (225)
 end
% ratio50=sum(Y_val(50,1:13))./Y_val(50,17)
% 
% ratio100=sum(Y_val(100,1:13))./Y_val(100,17)
% 
% ratio200=sum(Y_val(200,1:13))./Y_val(200,17)
% 
% ratio300=sum(Y_val(300,1:13))./Y_val(300,17)
% 
% B= [t_range',signalON];
% fileID = fopen('model1_010_Simulated.txt','w');
% fprintf(fileID,'%6.2f %12.8f\n',B');
% fclose(fileID);

% A=[t_range',Y_val(:,17)];
% fileID = fopen('Fibril_50.txt','w');
% fprintf(fileID,'%12.12f %12.12f \n',A');
% fclose(fileID);

end
