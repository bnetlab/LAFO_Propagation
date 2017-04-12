%looping for 3D plot
function value_LFAO

n=17;

x=67e-3;
% x=70e-3;  %10um
x1 =10e-3;
y=5e3;
y1=0.5e-1; 
z=1e5;
z1=5e-3;
p=5e3;
p1=1e-3;

A_1=0.5;
A_12=0.01;

theta=[x,x1,y,y1,z,z1,p,p1]; 
Y0=zeros(1,n); 
Y0(1)=A_12;
Y0(n)=A_1;
t_range=linspace(0,300,300); 
[t_val,Y_val]=ode23s(@ode_LFAO_1,t_range,Y0,[],n,theta);
Y_val([ 20:20:300 ],[1 4  12 13 14 16 17])

signalON=Y_val(:,n)*0;


for i=2:12
signalON=signalON + Y_val(:,i)*i;
end

signalON=signalON + Y_val(:,13)*24+Y_val(:,14)*36+Y_val(:,15)*48+Y_val(:,16)*65000;

signalON(end)
 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));

plot(t_range, signalON)
hold on;
load 'LFAO_DATA_01.txt';
Data=LFAO_DATA_01;
plot(Data(:,1),Data(:,2),'-*')

end
