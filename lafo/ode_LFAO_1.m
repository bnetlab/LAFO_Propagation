function dA_dt=ode_LFAO_1(t,A ,n,theta)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Defines the set of ODEs to be solved to simulate insulin fibrillation %
% A(1~n-1) is the vector of i-mer concentrations (i=1~n-1) %
% A(n) is fibril conc. and A(n+1) is natural hexamer conc. %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
dA_dt=zeros(size(A)); % First order derivatives of i-mer concentrations
% Flux of i-th fibrillation reaction
Jla=zeros(1,n);

knu=theta(1);
knu_=theta(2);
kla=theta(3); % First forward fibrillation rate constant
kla_=theta(4);
kfb1=theta(5);
kfb1_=theta(6);
kfb2=theta(7);
kfb2_=theta(8);


% Definitions of reaction fluxes Jfb
for i=1:n-1
 Jnu(i)=knu*A(n)*A(i)-knu_*A(i+1); % The flux of i-mer nucleation rxn
end
for i=12:14
 Jla(i)=kla*A(12)*A(i)-kla_*A(i+1); % The flux of i-mer nucleation rxn
end

Jfb1=kfb1 * A(15)*A(n)- kfb1_ * A(16);
Jfb2= kfb2 *A(16) *A(n) -kfb2_ * A(16);
% There are n equations representing the conc. change of n species

dA_dt(1)=-Jnu(1); % Derivative of monomer conc.
for i=2:11 % from dimer to (n-1)-mer
 dA_dt(i)=-Jnu(i)+Jnu(i-1); % Derivatives of oligomer concentrations
end
dA_dt(12)=Jnu(11)-sum(Jla(12:14))-Jla(12);
for i=13:14 % from dimer to (n-1)-mer
 dA_dt(i)=-Jla(i)+Jla(i-1); % Derivatives of oligomer concentrations
end
dA_dt(15)= Jla(14)-Jfb1;
dA_dt(16)=Jfb1;
dA_dt(n)=-sum(Jnu)-Jfb1-Jfb2;
