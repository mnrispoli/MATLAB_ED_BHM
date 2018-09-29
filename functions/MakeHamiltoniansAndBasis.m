function [Hi,Hj,basis] = MakeHamiltoniansAndBasis(NSites,NPart,jn,isPB)

%NSites = Number of sites in bose hubbard chain
%NPart = Number of bosons in total
%jn order of the tunneling, nearest neighbor = 1
%isPB determines boundary conditions: isPB = 1 means periodic boundary
%conditions

%find hilbert space size
HilbD=HilbDim(NSites,NPart);

% tunneling matrix
Hj=zeros(HilbD,HilbD);
Hi=Hj;

% creates a basis set that is written as a vector that is the hilbert space
% dimension x number of sites (way of writing down the many-body fock
% state)

basis=BasisMake(NSites,NPart);


%find elements coupled between sites for tunneling hamiltonian
%depends on boundary conditions
if isPB == 1
    %periodic boundary conditions
    EndSite=NSites;
else
    %open boundary conditions
    EndSite=NSites-jn;
end

for ss=1:EndSite
    
    %look at the pairs of sites coupled by the order (jn) of tunneling we
    %care about. Will make two 2D planes that are the hilbert space size
    %and find all pairs of sites i and j that differ by +/- 1 and also have
    %all other sites the same
    
    %site i for coupling
    si=ss;
    %site j for coupling
    sj=mod(si+jn-1,NSites)+1;
    
    %array initial states for sites i and j
    vbi=basis(:,[si,sj]);
    vi=cat(3,vbi(:,1)*ones(1,HilbD),vbi(:,2)*ones(1,HilbD));
    
    %array of final states for states i and j
    vbj=basis(:,[si,sj]);
    vj=cat(3,ones(HilbD,1)*vbj(:,1)',ones(HilbD,1)*vbj(:,2)');
    dvij=vi-vj;
    temp=dvij(:,:,1).*dvij(:,:,2);
    temp(find(temp~=-1))=0;
    dfij=abs(temp); % make all values 1

    % check now if all other sites are the same
    NSList=1:NSites;
    inds=intersect(find(NSList~=si),find(NSList~=sj)); % all remaining sites

    % go through remaining sites and make masks for when they don't change
    for ii=inds
        vii=basis(:,ii);
        dfii=abs(vii*ones(1,HilbD)-ones(HilbD,1)*vii');
        dfi=zeros(size(dfii));
        dfi(find(dfii==0))=1;
        dfij=dfij.*dfi;
    end
    
    %zero out bases elements that have more than one change
    dfij=abs(dfij);
    Hj=Hj+dfij;
end

%find linear list of elements
ijIndList=find(Hj>0);
for ijInd=1:length(ijIndList);
    ijElem=ijIndList(ijInd);
    [ii,jj]=ind2sub(size(Hj),ijElem);
    
    %only do upper right of matrix
    if ii>jj
        statei=basis(ii,:);
        statej=basis(jj,:);
        dstate=statei-statej;
        Hj(ii,jj)=-sqrt(statej(find(dstate==-1))).*sqrt(statei(find(dstate==1)));
    else
        
        statei=basis(ii,:);
        statej=basis(jj,:);
        dstate=statei-statej;
        Hj(ii,jj)=-sqrt(statej(find(dstate==-1))).*sqrt(statei(find(dstate==1)));
    end

end

%fix lower half with the transpose;
%Hj=Hj+ctranspose(Hj);

% make interaction hamiltonian (it's diagonal so it's super easy)
for hh=1:HilbD
    Hi(hh,hh)=sum(basis(hh,:).*(basis(hh,:)-1));
end

end