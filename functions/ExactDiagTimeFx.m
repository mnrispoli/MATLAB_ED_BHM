function [psiAll, PhiN, En, Cn] = ExactDiagTimeFx(psiInit,Ts,Ham)

%this function just takes the Hamiltonian "Ham" and calculates the dynamics
%of the initial state PsiInit

%find column eigenvectors PhiN from  V(:,index) and eigenvalues En from diag(D)
[V,D]=eig(Ham);
%eigenvalues
En=diag(D);
%eigenvectors
PhiN=V;

%overlap amplitudes with initial state
Cn=conj(psiInit)*PhiN;  
%make it a matrix
cnsMat=Cn'*ones(1,length(Ts));

%compute time dynamics from phase wrapping of eigenvalues for each
%eigenstate and taking care of the projection with the matrix cnsMat

psiAll(:,:)=(PhiN* ...
            (cnsMat.* ...
                exp(i.*2.*pi.*(En*Ts)) ...
            ) ...
            );

end