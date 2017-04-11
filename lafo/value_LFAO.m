%looping for 3D plot
function value_LFAO

n=17;

x=63e-3;
x1 =4e-3;
y=10000e-1;
y1=1e-6; 
z=1e6;
z1=5e-1;
p=1e3;
p1=5e-3;

A_1=0.5;
A_12=0.1;

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

signalON=signalON + Y_val(:,13)*24+Y_val(:,14)*36+Y_val(:,15)*48+Y_val(:,16)*14000;

signalON(end)
 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));

plot(t_range, signalON)
hold on;
load 'LFAO_DATA.txt';
Data=LFAO_DATA;
plot(Data(:,1),Data(:,2),'-*')

end
