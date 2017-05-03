%looping for 3D plot
% Added Fragmentation
function value_LFAO

n=30; %gateway 8

x=46e-3; %nu
x1 =10e-3;
y=5e6; %ilafo
y1=5e-2; 
z=5e5; %plafo
z1=5e-3;
p=5e5; %fib
p1=5e-3; 
q=5e3; % fag
q1=0;

A_1=0.5;
A_12=0.01;

theta=[x,x1,y,y1,z,z1,p,p1,q,q1]; 
Y0=zeros(1,n); 

Y0(n)=A_1;
Y0(1)=A_12;
t_range=linspace(0,337,337); 
[t_val,Y_val]=ode23s(@ode_LFAO_1,t_range,Y0,[],n,theta);
Y_val([1:20:300 ],[1 13 16 17  n-1 n])

signalON=Y_val(:,n)*0;


for i=2:13
signalON=signalON + Y_val(:,i)*(i-1);
end

for i=14:16
signalON=signalON + 12*(i-12);
end

for i=17:29
signalON=signalON +Y_val(:,i)*1700000;
end


signalON(end)
 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));

plot(t_range, signalON)
hold on;
load 'LFAO_DATA_01.txt';
Data=LFAO_DATA_01;
plot(Data(:,1),Data(:,2),'-*')

% X=Data(:,2);
% Y=signalON(Data(:,1)+1);
% mdl = fitlm(Y,X)

% X1=X(5:7);
% Y1=signalON(Data(5:7,1)+1);
% md2 = fitlm(Y1,X1)
% 
% signalON (192)/signalON (144)
% 
% signalON (300)/signalON (225)

% ratio50=sum(Y_val(50,1:13))./Y_val(50,17)
% 
% ratio100=sum(Y_val(100,1:13))./Y_val(100,17)
% 
% ratio200=sum(Y_val(200,1:13))./Y_val(200,17)
% 
% ratio300=sum(Y_val(300,1:13))./Y_val(300,17)
% 
% B= [t_range',signalON];
% fileID = fopen('LAFO_001_Simulated.txt','w');
% fprintf(fileID,'%6.2f %12.8f\n',B');
% fclose(fileID);

% A=[t_range',Y_val(:,17)];
% fileID = fopen('Fibril_50.txt','w');
% fprintf(fileID,'%12.12f %12.12f \n',A');
% fclose(fileID);

end