%looping for 3D plot
% New model with monomer elogation for pLAFO formation.
function M_value_LFAO

n=87; %gateway 10

x=30e-3;
x1 =10e-3;
y=5e2;
y1=5e-3; 
z=1e5;
z1=5e-3;


A_1=0.5;
A_12=0.1;

theta=[x,x1,y,y1,z,z1]; 
Y0=zeros(1,n); 

Y0(n)=A_1;
Y0(1)=A_12;
t_range=linspace(0,337,337); 
[t_val,Y_val]=ode23s(@M_ode_LFAO_2,t_range,Y0,[],n,theta);
Y_val([1:20:300 ],[1 4  13 85 n-1 n ])

signalON=Y_val(:,n)*0;


for i=2:85
signalON=signalON + Y_val(:,i)*(i-1);
end

signalON=signalON + Y_val(:,n-1)*2500000;

signalON(end)
 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));

plot(t_range, signalON)
hold on;
load 'LFAO_DATA.txt';
Data=LFAO_DATA;
plot(Data(:,1),Data(:,2),'-*')

X=Data(:,2);
Y=signalON(Data(:,1)+1);
mdl = fitlm(Y,X)

signalON (175)/signalON (125)
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
