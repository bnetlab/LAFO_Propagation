function dA_dt=ode_LFAO_1(t,A ,n,theta)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Added Fragmentation
% Defines the set of ODEs to be solved to simulate insulin fibrillation %
% A(1~n-1) is the vector of i-mer concentrations (i=1~n-1) %
% A(n) is fibril conc. and A(n+1) is natural hexamer conc. %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
dA_dt=zeros(size(A)); % First order derivatives of i-mer concentrations
% Flux of i-th fibrillation reaction
Jla=zeros(1,n);
Jfb3=zeros(1,n);

knu=theta(1);
knu_=theta(2);
kla=theta(3); % First forward fibrillation rate constant
kla_=theta(4);
kfb1=theta(5);
kfb1_=theta(6);
kfb4=theta(7);
kfb4_=theta(8);
kfb3=theta(9);
kfb3_=theta(10);
kfag=theta(11);
kfag_=theta(12);

% Definitions of reaction fluxes Jfb
for i=1:12
 Jnu(i)=knu*A(n)*A(i)-knu_*A(i+1); % The flux of i-mer nucleation rxn
end
for i=13:n-15 %
 Jla(i)=kla*A(13)*A(i)-kla_*A(i+1); % The flux of i-mer nucleation rxn
end

% Jfag=kfag1 * A(13) - kfag1_ * A(1).^2;
Jfag=0;
Jfb1=kfb1 * A(n-14)*A(n-14)- kfb1_ * A(n-13); %
Jfb2= kfb1 *A(n-14) *A(n-13)-kfb1_ * A(n-13);
Jfb4a= kfb4 *A(13) *A(n-14)-kfb4_ * A(n-13);%
Jfb4b= kfb4/2 *A(13) *A(n-13)-kfb4_ * A(n-13);

for i=n-13:n-2
 Jfb3(i)=kfb3*A(i)*A(n)-kfb3_*A(i+1); % The flux of i-mer nucleation rxn
end

Jfb4=kfag*A(n-1);

% There are n equations representing the conc. change of n species

dA_dt(1)=-Jnu(1)+2*Jfag+Jfb4; % Derivative of monomer conc.
for i=2:12 % from dimer to (n-1)-mer
 dA_dt(i)=-Jnu(i)+Jnu(i-1); % Derivatives of oligomer concentrations
end
dA_dt(13)=Jnu(12)-sum(Jla)-Jla(13)-Jfag-Jfb4a-Jfb4b; %
for i=14:n-15 % from dimer to (n-1)-mer %
 dA_dt(i)=-Jla(i)+Jla(i-1); % Derivatives of oligomer concentrations
end
dA_dt(n-14)= Jla(n-15)-2*Jfb1-Jfb2-Jfb4a;
dA_dt(n-13)=Jfb1-Jfb3(n-13)+Jfb4+Jfb4a;

for i=n-12:n-2% from dimer to (n-1)-mer
 dA_dt(i)=-Jfb3(i)+Jfb3(i-1); % Derivatives of oligomer concentrations
end
dA_dt(n-1)=Jfb3(n-2)-Jfb4;
dA_dt(n)=-sum(Jnu)-sum(Jfb3);
end