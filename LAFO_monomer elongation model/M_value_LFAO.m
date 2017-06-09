%looping for 3D plot
% New model with monomer elogation for pLAFO formation.
function M_value_LFAO

n=39; %gateway 10

x=30e-3;
x1 =0.5e-3;
y=2.23e-1;
y1=1.6e-1; 
z=3e5;
z1=5e-3;
p=0;
p1=0;


A_1=0.5;
B_12=[0.0068 0.002 0.0005];
C_24=[0.0466 0.004 .00025 ];

for j=1:3
    A_12=B_12(j);
    D_24=C_24(j);

theta=[x,x1,y,y1,z,z1,p,p1]; 
Y0=zeros(1,n); 

Y0(n)=A_1;
Y0(1)=A_12;
Y0(13)=D_24;

t_range=linspace(0,337,337); 
[t_val,Y_val]=ode23s(@M_ode_LFAO_2,t_range,Y0,[],n,theta);
Y_val([1:20:300 ],[1 4 12 13 14 n-2  n-1 n ])

signalON=Y_val(:,n)*0;


for i=13:n-2
signalON=signalON + Y_val(:,i)*(i-1);
end

signalON=signalON + Y_val(:,n-1)*10000000;

signalON(end)
 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));

plot(t_range, signalON)
hold on;

if (j==1)
     load 'LFAO_DATA.txt';
    Data=LFAO_DATA;
    plot(Data(:,1),Data(:,2),'-*')
%     s=(Data(7,2)-Data(5,2))./(signalON(193)-signalON(145))
%   B= [t_range',signalON];
%  fileID = fopen('Model2_10_Simulated.txt','w');
%  fprintf(fileID,'%6.2f %12.8f\n',B');
%  fclose(fileID);
elseif (j==2)
    load 'LFAO_DATA_01.txt';
    Data=LFAO_DATA_01;
    plot(Data(:,1),Data(:,2),'-*')
%     s=(Data(8,2)-Data(5,2))./(signalON(217)-signalON(145))
%     B= [t_range',signalON];
% fileID = fopen('Model2_1_Simulated.txt','w');
% fprintf(fileID,'%6.2f %12.8f\n',B');
% fclose(fileID);
 else
    load 'LFAO_DATA_00001.txt';
   Data=LFAO_DATA_00001;
     plot(Data(:,1),Data(:,2),'-*')
%     s=(Data(11,2)-Data(8,2))./(signalON(241)-signalON(169))
%     B= [t_range',signalON];
% fileID = fopen('Model2_01_Simulated.txt','w');
% fprintf(fileID,'%6.2f %12.8f\n',B');
% fclose(fileID);
 end

X=Data(:,2);
Y=signalON(Data(:,1)+1);
mdl = fitlm(Y,X)

signalON (end)/signalON (250)


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
% fileID = fopen('Model2_001_Simulated.txt','w');
% fprintf(fileID,'%6.2f %12.8f\n',B');
% fclose(fileID);

% A=[t_range',Y_val(:,17)];
% fileID = fopen('Fibril_50.txt','w');
% fprintf(fileID,'%12.12f %12.12f \n',A');
% fclose(fileID);

end
