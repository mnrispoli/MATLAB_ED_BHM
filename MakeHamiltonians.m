function [Hi,Hj,Basis] = MakeHamiltoniansAndBasis(HilbD,NSites,NPart)


% creates a basis set that is written as a vector that is the hilbert space
% dimension x number of sites (way of writing down the many-body fock
% state)

basis=zeros(HilbD,NSites);

% initial state for the function to make all other basis states from
% starts with all the bosons on the edge-most site
basis(1,1)=NPart;

% tunneling matrix
Hj=zeros(HilbD,HilbD);

cc=1;
for k=1:NSites-1
    
    %find states where site-k has an atom
    kset=find(basis(:,k)~=0)';
    
    %iterate through all states where site-k has an atom
    for kk=kset
        %number of atoms on site-k in basis state kk
        bb=basis(kk,k);
        
        TiInds=[kk cc+1:cc+maskY-1];
        TfInds=[cc+1:cc+maskY];
        
        %mask for moving all of the atoms off of this site to the
        %neighboring one
        mask=[bb-1:-1:0;1:bb]';
        
        %set the y size from the mask
        maskY=size(mask,1);
        
        %add this basis of moved particles by one to the basis to the left
        %of site-k
        basis(cc+1:cc+maskY,1:(k-1))=ones(maskY,1)*basis(kk,1:(k-1));
        
        %add on the particles moved to the right of the basis
        basis(cc+1:cc+maskY,k:k+1)=mask;
        
        Hj(sub2ind(size(Hj),TfInds,TiInds))=1
        %number of basis states moved down (unique ones were made above)
        cc=cc+bb;
    end
        
end



end