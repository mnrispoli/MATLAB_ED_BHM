function [sn2k,S2kInds,indsNKStore,corrn,gn] = dynBldGnFxFsPPMapPDF(pdf,basis,n,inds,PP)

NPart=mean(sum(basis,2),1);
%PP=1

S2kInds=cell(n,n); % create cell structure to find all groupings as following the Stirling Numbers of the second kind
gn=cell(1,n); % create a cell to save high n-order correlation functions
sn2k=zeros(n,n); % just the number of unique partitions for nparts into k groups

for ni=1:n
    %tic %timing check
    sn2k(ni,1)=1; %initialize for 1 unique grouping of n parts into 1 group
    sn2k(ni,ni)=1; % initialize for 1 unique grouping of n parts into n groups
    
    S2kInds{ni,1}{1,1}=[1:ni];
    
    for ki=2:ni
        %number of permutations at given Stirling's number of the second
        %kind      
        
        sn2k(ni,ki)=sn2k(ni-1,ki-1)+ki*sn2k(ni-1,ki);
        
        temp=cell(sn2k(ni,ki),ki);
        
        ref1=S2kInds{ni-1,ki-1};
        [nx1,ny1]=size(ref1);
        ref2=S2kInds{ni-1,ki};
        [nx2,ny2]=size(ref2);
        
        for ii=1:nx1
            for jj=1:ki-1
                temp{ii,jj}=ref1{ii,jj};
            end
            temp{ii,ki}=[ni];
        end
        
        for ii=1:nx2
            for rr=1:ki
            for jj=1:ki
                if jj==rr
                    temp{nx1+(ii-1)*ki+rr,jj}=[ref2{ii,jj} ni];
                else
                    temp{nx1+(ii-1)*ki+rr,jj}=ref2{ii,jj};
                end
            end
            end
        end
        

        
        S2kInds{ni,ki}=temp;
        
        
        
    end
    
    indsNK=sort(nchoosek(inds,ni),2);
    indsNKStore{ni}=indsNK;
    
    if ni==1
        %gntemp=zeros([1 n]);
        %ontemp=zeros([1 n]);
        gnTemp=containers.Map();
        onTemp=containers.Map();
        
        for aa=1:size(indsNK,1)

            cellIndsStr=num2str(indsNK(aa,:));

            %gnTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1);
            %onTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1);
            
            if PP==1
                gnTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1)./NPart;
                onTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1)./NPart;
            else
                gnTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1);
                onTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1);
            end
        end
        %gn{ni}=gntemp;
        %corrn{ni}=ontemp;
        
        gn{ni}=gnTemp;
        corrn{ni}=onTemp;
        
    else
        %gntemp=zeros(ones(1,ni).*n);
        %ontemp=zeros(ones(1,ni).*n);
        
        gnTemp=containers.Map();
        onTemp=containers.Map();
        
        
        for aa=1:size(indsNK,1)

            cellIndsStr=num2str(indsNK(aa,:));
            
            %gnTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1);
            %onTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1);
            
            ppF=factorial(NPart)/factorial(NPart-length(indsNK(aa,:)));
            if PP==1
                gnTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1)./ppF;
                onTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1)./ppF;
            else
                gnTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1);
                onTemp(cellIndsStr)=sum(pdf.*prod(basis(:,indsNK(aa,:)),2),1);
            end
            
            for s2ni=2:ni
                tempCell=S2kInds{ni,s2ni};
                [nterms, nmult]=size(tempCell);
                for nti=1:nterms
                    tempVal=1;
                    %smult=13;
                    for nmi=1:nmult
                        indTerms=tempCell{nti,nmi};
                        refGn=gn{length(indTerms)};
                        
                        cellIndsStr2=num2str(indsNK(aa,indTerms(:)));

                        tempVal=tempVal.*refGn(cellIndsStr2);
                    end

                    %gntemp(cellInds{:})=gntemp(cellInds{:})-tempVal;
                    gnTemp(cellIndsStr)=gnTemp(cellIndsStr)-tempVal;
                end

            end
        end

        gn{ni}=gnTemp;
        corrn{ni}=onTemp;
    end
    
    %toc
end

end