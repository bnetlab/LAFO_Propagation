function dA_dt=ode_LFAO_1(t,A ,n,theta)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Added Fragmentation
% Defines the set of ODEs to be solved to simulate insulin fibrillation %
% A(1~n-1) is the vector of i-mer concentrations (i=1~n-1) %
% A(n) is fibril conc. and A(n+1) is natural hexamer conc. %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
dA_dt=zeros(size(A)); % First order derivatives of i-mer concentrations
% Flux of i-th fibrillation reaction
n1=39;
Jla=zeros(1,n);
Jnuon=zeros(1,n);
Jfbon=zeros(1,n);

knu=theta(1);
knu_=theta(2);
kla=theta(3); % First forward fibrillation rate constant
kla_=theta(4);
kfb1=theta(5);
kfb1_=theta(6);

knuon=theta(7); 
knuon_=theta(8); % Reverse nucleation constants
kfbon=theta(9); % First forward fibrillation rate constant
kfbon_=theta(10); % Reverse fibrillation rate constant


% Definitions of reaction fluxes Jfb
for i=1:12
 Jnu(i)=knu*A(n1)*A(i)-knu_*A(i+1); % The flux of i-mer nucleation rxn
end
for i=13:n1-3 %
 Jla(i)=kla*A(i)*A(n1)-kla_*A(i+1); % The flux of i-mer nucleation rxn
end

Jfb1=kfb1 * A(n1-2)*A(n1)- kfb1_ * A(n1-1); %
Jfb2= kfb1 *A(n1-1) *A(n1).^0.8-kfb1_ * A(n1-1); %

for i=n1:n-1
  Jnuon(i)=knuon*A(n1)*A(i)-knuon_*A(i+1); % The flux of i-mer nucleation rxn
  Jfbon(i)=kfbon*A(n)*A(i)-kfbon_*A(n); % The flux of i-mer elongation rxn
end


% There are n equations representing the conc. change of n species

dA_dt(1)=-Jnu(1); % Derivative of monomer conc.
for i=2:12 % from dimer to (n-1)-mer
 dA_dt(i)=-Jnu(i)+Jnu(i-1); % Derivatives of oligomer concentrations
end
dA_dt(13)=Jnu(12)-Jla(13); %
for i=14:n1-3 % from dimer to (n-1)-mer %
 dA_dt(i)=-Jla(i)+Jla(i-1); % Derivatives of oligomer concentrations
end
dA_dt(n1-2)= Jla(n1-3)-Jfb1; %
dA_dt(n1-1)=Jfb1; %
dA_dt(n1)=-sum(Jnu)-sum(Jla)-Jfb1-Jfb2-sum(Jnuon(18:28))-Jnuon(18)-Jfbon(18);
for i=n1+1:n-1 % from dimer to (n-1)-mer
dA_dt(i)=-Jnuon(i)+Jnuon(i-1)-Jfbon(i); % Derivatives of oligomer concentrations
end
dA_dt(n)=Jnuon(n-1);

end