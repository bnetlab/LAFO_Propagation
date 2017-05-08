%looping for 3D plot
% Added Fragmentation
function value_LFAO

n=16; % gateway 10

x=26e-3;
x1 =2e-3;
y=4e6;
y1=5e-2; 
z=5e5;
z1=5e-3;
p=3e6;
p1=5e-3; 



A_1=0.5;
B_12=[0.1 0.01 0.001];

for j=1:3
     A_12=B_12(j);

theta=[x,x1,y,y1,z,z1,p,p1]; 
Y0=zeros(1,n); 

Y0(n)=A_1;
Y0(1)=A_12;
t_range=linspace(0,337,337); 
[t_val,Y_val]=ode23s(@ode_LFAO_1,t_range,Y0,[],n,theta);
Y_val([1:20:300 ],[1 4  13 n-1 n ])

signalON=Y_val(:,n)*0;

for i=2:13
signalON=signalON + Y_val(:,i)*(i-1);
end

for i=14:n-2
 signalON=signalON + 12*(i-12);
end

signalON=signalON +Y_val(:,n-1)*20000000;

signalON(end)
 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));

plot(t_range, signalON)
hold on;

if (j==1)
     load 'LFAO_DATA.txt';
    Data=LFAO_DATA;
    plot(Data(:,1),Data(:,2),'-*')
    s=(Data(7,2)-Data(5,2))./(signalON(193)-signalON(145))
elseif (j==2)
    load 'LFAO_DATA_01.txt';
    Data=LFAO_DATA_01;
    plot(Data(:,1),Data(:,2),'-*')
    s=(Data(8,2)-Data(5,2))./(signalON(217)-signalON(145))
 else
    load 'LFAO_DATA_00001.txt';
   Data=LFAO_DATA_00001;
     plot(Data(:,1),Data(:,2),'-*')
    s=(Data(11,2)-Data(8,2))./(signalON(241)-signalON(169))
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
% fileID = fopen('model1_010_Simulated.txt','w');
% fprintf(fileID,'%6.2f %12.8f\n',B');
% fclose(fileID);

% A=[t_range',Y_val(:,17)];
% fileID = fopen('Fibril_50.txt','w');
% fprintf(fileID,'%12.12f %12.12f \n',A');
% fclose(fileID);

end
