% create subsystems n shit
load('8site_loadme.mat')

halfS=size(basis,2)/2;
fullS=size(basis,2);
hD=size(basis,1);


masksA=zeros(halfS,fullS);
for ii=1:halfS
    masksA(ii,1:ii)=ones(size(masksA(ii,1:ii)))
end

masksAt=masksA;

masksA=masksAt(end,:);

EofPInd=size(masksA,1)+1;

 
masksB=ones(size(masksA))-masksA;

expon=(fullS:-1:1)-1;
bases=ones(size(masksA,1),1)*(10.^(expon));

masksD=size(masksA,1);

subSys=cell(1,masksD);
uniqueRedSys=cell(1,masksD);
psiInds=cell(1,masksD);

subA=zeros(masksD,hD);
subB=subA;

subN=zeros(masksD,hD);
for ii=1:hD
    subA(:,ii)=sum(masksA.*bases.*(ones(masksD,1)*basis(ii,:)),2);
    subB(:,ii)=sum(masksB.*bases.*(ones(masksD,1)*basis(ii,:)),2);
    
    subN(:,ii)=sum(masksA.*(ones(masksD,1)*basis(ii,:)),2);
end

allUSubA=cell(1,masksD);
allUSubB=cell(1,masksD);

allUindsSubA=cell(1,masksD);
allUindsSubB=cell(1,masksD);

% create mapping matrices
AMap=cell(masksD,1);
BMap=cell(masksD,1);
BAMap=cell(masksD,1);
rhoMap=cell(masksD,1);

for ii=1:masksD
        %unrestricted n subspace

        allUSubA{ii}=unique(subA(ii,:));
        tempA=allUSubA{ii};

        allUSubB{ii}=unique(subB(ii,:));
        tempB=allUSubB{ii};

        BAMap{ii}=zeros(length(allUSubB{ii}),length(allUSubA{ii}));

        AMap{ii}=zeros(hD,1);
        BMap{ii}=zeros(hD,1);

        tempBA=BAMap{ii};

        tempAM=AMap{ii};
        tempBM=BMap{ii};

        for jj=1:length(allUSubA{ii})
            indsA=find(tempA(jj)==subA(ii,:));
            tempAM(indsA)=jj;
        end
        
        
        for jj=1:length(allUSubB{ii})
            indsB=find(tempB(jj)==subB(ii,:));
            tempBM(indsB)=jj;
            
            for kk=1:length(indsB)
                indsA=find(subA(ii,indsB(kk))==allUSubA{ii});
                tempBA(jj,indsA)=indsB(kk);
                %tempAM(jj,indsA)=indsB(kk);
            end
        end
        BMap{ii}=tempBM;
        AMap{ii}=tempAM;
        BAMap{ii}=tempBA;
        
end  

%okay go through all the psis now?
NW=size(psiAllW,1);
ND=size(psiAllW,2);
NT=size(psiAllW{1,1},1);

SvnSubs=cell(NW,ND,masksD);
SvnSubsAvg=cell(NW,masksD);
rhoSubs=cell(NW,ND,masksD);
%%
%SPMTemp=
h=waitbar(0,'tracing')
for ii=1:masksD
    ii
    rdRhoMap=BAMap{ii};
    hdA=size(rdRhoMap,2);
    tempBA=BAMap{ii};
    
    for ww=1:length(Ws)
        
        
        SvnSubsAvg{ww,ii}=zeros(1,NT);

        saveVs=cell(ND,masksD);
        saveDs=cell(ND,masksD);
        
        for dd=1:ND
            ww
            dd
            
            rhoA=zeros(hdA,hdA,NT);
            psiT=psiAllW{ww,dd};
            
            for bb=1:size(tempBA,1)

                indsA=find(tempBA(bb,:)>0);
                
                clear rho
                rho=zeros(hdA,hdA);
                
                for tt=1:NT
                    rho(indsA,indsA)=ctranspose(psiT(tt,tempBA(bb,indsA)))* ...
                        psiT(tt,tempBA(bb,indsA));
                    
                    rhoA(:,:,tt)=rhoA(:,:,tt)+ ...
                        rho;

                    
                end
                
            end
            

            
            %clear rhoiiSave
            rhoiiSave=zeros(hdA,NT);
            for tt=1:NT
                               
                temprho=reshape(rhoA(:,:,tt),[size(rhoA,1) size(rhoA,2)]);
                [V,D]=eig(temprho);
                rhoii=diag(abs(ctranspose(V(:,:))*temprho*V(:,:)));
                
                indrho=find(rhoii>0);
                Svn(tt)=-sum(rhoii(indrho).*log(rhoii(indrho)));
                
                rhoiiSave(:,tt)=rhoii;
                
            end

            
            rhoiiSaveW{ww,dd,ii}=rhoiiSave;
            SvnSubs{ww,dd,ii}=Svn;
            SvnSubsAvg{ww,ii}=SvnSubsAvg{ww,ii}+SvnSubs{ww,dd,ii}./ND;
            
            waitbar(((dd-1)+(ww-1)*(ND)+(ii-1)*(NW*ND))/(masksD*NW*ND),h)
        end

        
    end
    
end
close(h)


%%
svnM=zeros(NT,length(Ws));
for ww=1:length(Ws)
    svnM(:,ww)=SvnSubsAvg{ww,1};
end
%%
figure(22)
subplot(1,3,1)
hold on
semilogx(tj,SPM)
subplot(1,3,2)
hold on
plot(tj,svnM)
subplot(1,3,3)
hold on
plot(tj,svnM-SPM)

