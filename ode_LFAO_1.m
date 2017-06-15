function dA_dt=ode_LFAO_1(t,A ,n,theta)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Added Fragmentation
% Defines the set of ODEs to be solved to simulate insulin fibrillation %
% A(1~n-1) is the vector of i-mer concentrations (i=1~n-1) %
% A(n) is fibril conc. and A(n+1) is natural hexamer conc. %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
dA_dt=zeros(size(A)); % First order derivatives of i-mer concentrations
% Flux of i-th fibrillation reaction
n1=n-9;
Jla=zeros(1,n);
Jnuon=zeros(1,n);
Jfbon=zeros(1,n);

knu=theta(1);
knu_=theta(2);
kla=theta(3); % First forward fibrillation rate constant
kla_=theta(4);
kfb1=theta(5);
kfb1_=theta(6);
kfb3=theta(7);
kfb3_=theta(8);
kfag=theta(9);

knuon=theta(10); 
knuon_=theta(11); % Reverse nucleation constants
kfbon=theta(12); % First forward fibrillation rate constant
kfbon_=theta(13); % Reverse fibrillation rate constant


% Definitions of reaction fluxes Jfb
for i=1:12
 Jnu(i)=knu*A(n1)*A(i)-knu_*A(i+1); % The flux of i-mer nucleation rxn
end
for i=13:n1-15 %
 Jla(i)=kla*A(i)*A(n1)-kla_*A(i+1); % The flux of i-mer nucleation rxn
end

Jfb1=kfb1 * A(n1-14)*A(n1-14)- kfb1_ * A(n1-13); %
Jfb2= kfb1/10 *A(n1-14) *A(n1-13)-kfb1_ * A(n1-13); %

for i=n1-13:n1-2
 Jfb3(i)=kfb3*A(i)*A(n1)-kfb3_*A(i+1); % The flux of i-mer nucleation rxn
end

Jfb4=kfag*A(n1-1);

for i=n1:n-1
  Jnuon(i)=knuon*A(n1)*A(i)-knuon_*A(i+1); % The flux of i-mer nucleation rxn
  Jfbon(i)=kfbon*A(n)*A(i)-kfbon_*A(n); % The flux of i-mer elongation rxn
end


% There are n equations representing the conc. change of n species

dA_dt(1)=-Jnu(1)+Jfb4; % Derivative of monomer conc.
for i=2:12 % from dimer to (n-1)-mer
 dA_dt(i)=-Jnu(i)+Jnu(i-1); % Derivatives of oligomer concentrations
end
dA_dt(13)=Jnu(12)-Jla(13); %
for i=14:n1-15 % from dimer to (n-1)-mer %
 dA_dt(i)=-Jla(i)+Jla(i-1); % Derivatives of oligomer concentrations
end
dA_dt(n1-14)= Jla(n1-15)-2*Jfb1-Jfb2; %
dA_dt(n1-13)=Jfb1-Jfb3(n1-13)+Jfb4; 

for i=n1-12:n1-2% from dimer to (n-1)-mer
 dA_dt(i)=-Jfb3(i)+Jfb3(i-1); % Derivatives of oligomer concentrations
end
dA_dt(n1-1)=Jfb3(n1-2)-Jfb4;

dA_dt(n1)=-sum(Jnu)-sum(Jla)-sum(Jfb3)-sum(Jnuon)-Jnuon(n1)-Jfbon(n1);
for i=n1+1:n-1 % from dimer to (n-1)-mer
dA_dt(i)=-Jnuon(i)+Jnuon(i-1)-Jfbon(i); % Derivatives of oligomer concentrations
end
dA_dt(n)=Jnuon(n-1);

end