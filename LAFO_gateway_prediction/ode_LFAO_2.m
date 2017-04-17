function dA_dt=ode_LFAO_1(t,A ,n,theta)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%Gateway 6

% Added Fragmentation
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
kfag=theta(9);
kfag_=theta(10);
kfb3=theta(11);
kfb3_=theta(12);


% Definitions of reaction fluxes Jfb
for i=1:11
 Jnu(i)=knu*A(n)*A(i)-knu_*A(i+1); % The flux of i-mer nucleation rxn
end
for i=12:16 %
 Jla(i)=kla*A(12)*A(i)-kla_*A(i+1); % The flux of i-mer nucleation rxn
end

Jfag=kfag * A(12) - kfag_ * A(1).^2;
Jfb1=kfb1 * A(17)*A(n)- kfb1_ * A(18); %
Jfb2= kfb2 *A(18) *A(n)-kfb2_ * A(18); %

Jfb3=kfb3 * A(17)*A(12)- kfb3_ * A(18); %


% There are n equations representing the conc. change of n species

dA_dt(1)=-Jnu(1)+2 * Jfag; % Derivative of monomer conc.
for i=2:11 % from dimer to (n-1)-mer
 dA_dt(i)=-Jnu(i)+Jnu(i-1); % Derivatives of oligomer concentrations
end
dA_dt(12)=Jnu(11)-sum(Jla(12:16))-Jla(12)-Jfag-Jfb3; %
for i=13:16 % from dimer to (n-1)-mer %
 dA_dt(i)=-Jla(i)+Jla(i-1); % Derivatives of oligomer concentrations
end
dA_dt(17)= Jla(16)-Jfb1; %
dA_dt(18)=Jfb1+Jfb3; %
dA_dt(n)=-sum(Jnu)-Jfb1-Jfb2;