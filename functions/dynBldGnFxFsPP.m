function [sn2k,S2kInds,indsNKStore,corrn,gn] = dynBldGnFxFsPP(data,n,inds,PP)

%inds=1:nsites; % set the inds over which to find the correlations from the ground up
%basis=BasisMake(nsites,nsites); % example data set that is the infinite temperature ensemble for BH
%{
nsep=6;
if nsites>nsep
basis=BasisMake(nsep,nsep); % example data set that is the infinite temperature ensemble for BH
basis2=BasisMake(nsites-nsep,nsites-nsep);
data=[];
    for bb=1:size(basis,1)
        data=[data; ones(size(basis2,1),1)*basis(bb,:) basis2];
    end
else
    data=BasisMake(nsites,nsites); 
end
%}
%{
data=zeros(3,n);
data(1,1:2:end)=2;
data(2,2:2:end)=2;
data(3,:)=1;
%}
%{
nconsv=[2E1 1E2 1E3 1E4 1E5];
for nc=1:length(nconsv)
    clearvars -except nconsv nc gnmean
test=randi([1 12],12,nconsv(nc));
size(test)
datatest=zeros(size(test));
for jj=1:size(test,2)
for ii=1:size(test,1)
datatest(test(ii,jj),jj)=datatest(test(ii,jj),jj)+1;
end
end


%data=datanew;
data=datatest';
n=8
inds=1:12
%}
NPart=mean(sum(data,2),1);
%PP=1

S2kInds=cell(n,n); % create cell structure to find all groupings as following the Stirling Numbers of the second kind
gn=cell(1,n); % create a cell to save high n-order correlation functions
sn2k=zeros(n,n); % just the number of unique partitions for nparts into k groups

for ni=1:n
    tic %timing check
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
        gntemp=zeros([1 n]);
        ontemp=zeros([1 n]);
        
        for aa=1:size(indsNK,1)
            cellInds=cell(1,ni);
            for ii=1:ni
                cellInds{ii}=indsNK(aa,ii);
            end
            
                    
            
            if PP==1
                gntemp(cellInds{:})=mean(prod(data(:,indsNK(aa,:)),2),1)./NPart;
                ontemp(cellInds{:})=mean(prod(data(:,indsNK(aa,:)),2),1)./NPart;
            else
                gntemp(cellInds{:})=mean(prod(data(:,indsNK(aa,:)),2),1);
                ontemp(cellInds{:})=mean(prod(data(:,indsNK(aa,:)),2),1);
            end

        end
        gn{ni}=gntemp;
        corrn{ni}=ontemp;
        
    else
        gntemp=zeros(ones(1,ni).*n);
        ontemp=zeros(ones(1,ni).*n);
        
        for aa=1:size(indsNK,1)
            cellInds=cell(1,ni);
            for ii=1:ni
                cellInds{ii}=indsNK(aa,ii);
            end

            ppF=factorial(NPart)/factorial(NPart-length(indsNK(aa,:)));
            if PP==1
                gntemp(cellInds{:})=mean(prod(data(:,indsNK(aa,:)),2),1)./ppF;
                ontemp(cellInds{:})=mean(prod(data(:,indsNK(aa,:)),2),1)./ppF;
            else
                gntemp(cellInds{:})=mean(prod(data(:,indsNK(aa,:)),2),1);
                ontemp(cellInds{:})=mean(prod(data(:,indsNK(aa,:)),2),1);
            end
            %gn{ni}=gntemp;

            for s2ni=2:ni
                tempCell=S2kInds{ni,s2ni};
                [nterms, nmult]=size(tempCell);
                for nti=1:nterms
                    tempVal=1;
                    %smult=13;
                    for nmi=1:nmult
                        indTerms=tempCell{nti,nmi};
                        refGn=gn{length(indTerms)};
                        %smult=smult-length(indTerms); % from philipp's paper if I'm doing it right

                                cellInds2=cell(1,length(indTerms));
                                for ii=1:length(indTerms)
                                    cellInds2{ii}=indsNK(aa,indTerms(ii));
                                end


                        tempVal=tempVal.*refGn(cellInds2{:});

                            %tempVal=tempVal.*refGn(cellInds2{:});

                    end

                    gntemp(cellInds{:})=gntemp(cellInds{:})-tempVal;
                end

            end
        end

        gn{ni}=gntemp;
        corrn{ni}=ontemp;
    end
    
    toc
end


%end
%%
%{
for nn=1:n
    ppF=factorial(NPart)/factorial(NPart-nn);
    if PP==1
    temp=gn{nn}.*ppF;
    else
        temp=gn{nn}.*1;
    end
%    tempPP=gnPP{nn}.*ppF;

    gnmean(nc,nn)=mean(abs(temp(find(abs(temp)>0))));
%   gnmeanPP(ff,nn)=mean(abs(tempPP(find(abs(tempPP)>0))));
end
%}
end