function value100
% estimate error of one run
   
n=12;
x=3.2e-3;
y =1e-4;
z=1e12;
zz=1e-3;

A_1=0.5;
theta=[x,y,z,zz]; 
Y0=zeros(1,n); 
Y0(1)=A_1;


t_range=linspace(0,337,337); 
[t_val,Y_val]=ode23s(@lee_ode100,t_range,Y0,[],n,theta);
Y_val([1:20:337 ],[1 4 11 n])

signalON=Y_val(:,n)*0;
for i=n
signalON=signalON + Y_val(:,i);
end

signalON([1 25 end],1) 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));
plot(t_range,signalON)
hold on

load control_LAFO.txt;
Data=control_LAFO;
plot(Data(:,1),Data(:,2),'-*')

X=Data(:,2);
Y=signalON(Data(:,1)+1);
mdl = fitlm(Y,X)
sum(Y_val(end,2:11))