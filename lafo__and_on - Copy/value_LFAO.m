%looping for 3D plot
% Added Fragmentation
function value_LFAO

n=48; % gateway 10

x=30e-3;
x1 =0.5e-3;
y=2.23e-1;
y1=1.6e-1; 
z=3e5;
z1=5e-3;


u=6e-3;
u1 =1e-4;
v=5e6;
v1=1e-3;


A_1=0.5;
B_12=[0.0068 0.002 0.0005];
C_24=[0.0466 0.004 .00025 ];


for j=1:3
    A_12=B_12(j);
    D_24=C_24(j);

theta=[x,x1,y,y1,z,z1,u,u1,v,v1]; 
Y0=zeros(1,n); 


Y0(39)=A_1;
Y0(1)=A_12;
Y0(13)=D_24;

t_range=linspace(0,337,337); 
[t_val,Y_val]=ode23s(@ode_LFAO_1,t_range,Y0,[],n,theta);
Y_val([1:20:300 ],[1 13 25 38 48])

signalON=Y_val(:,n)*0;

for i=13:37
signalON=signalON + Y_val(:,i)*(i-1);
end

signalON=signalON + Y_val(:,38)*10000000;

signalON=signalON +Y_val(:,48)*10000000;

signalON(end)
 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));

plot(t_range, signalON)
hold on;

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
