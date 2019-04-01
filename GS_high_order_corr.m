% general ED bose-hubbard code
% last cleaned up 9/24/2018
% author: Matthew Rispoli

addpath('functions')

close all
clear all
clc


% tunneling
J=38.1911/1000; %8Er tunneling rate; defined in Hz/(2 pi)

% disorder strengths
NW=1;
Ws=linspace(0.3,22,NW).*0.*J % units of J
Ws=[1 4 8 12 16 20]
% tupe of disorder
isQP=1; % 0 = uniform random, 1 = FFT of disorder from exp, 2 = actual QP

% interaction strengths
Uo=2.7*J % units of J
Us=ones(size(Ws)).*Uo; %make them all the same
Us=linspace(0,10,NW).*J;

% system size parameters
NPart=8 % number of bosons
NSites=8 % number of sites

% num of disorders and time scans
ND=51; %disorder number
DisordN=ND
NT=151; % time steps
tau=(1./(2.*pi.*J));
Ts=linspace(0,1000,NT) %actual times for ED evaluation 
Ts=logspace(-1,4,NT).*tau;
% periodic boundary conditions
isPB = 0; % 0 = open boundary; 1 = periodic boundary
jn = 1; % order of neighbor tunnelings

%%

% folder for storing all matrix files
folderName = 'Hamiltonians/'
% hamiltonian and basis parameters
edFileName=[folderName sprintf('%i_sites_%i_bosons_%i_pb_%i_jn_loadme.mat',NSites,NPart,isPB,jn)] %format for Hamiltonian filename

% find hilbert space dimension
HilbD=HilbDim(NSites,NPart)

% check if we've already made the basis and hamiltonians
if exist(edFileName,'file')==2
    % if it exists then load
    load(edFileName);
else
    % if it doesn't exist then make them all
    [Hi,Hj,basis] = MakeHamiltoniansAndBasis(NSites,NPart,jn,isPB);
    save(edFileName,'Hi','Hj','basis','NPart','NSites','jn','isPB')
end

%%
% make giant dataStruct to save wavefunction for all disorders
psiAllW=cell(NW,ND);

% find initial state in basis
psiInit=zeros(1,HilbD); 
psiInit(find(ismember(basis,ones(1,NSites),'rows')==1))=1;
%%
gnPDFStoreStd=zeros(NT,ND,NSites,NW);
gnPDFStoreMAbs=zeros(NT,ND,NSites,NW);
for ww=1:length(Ws)
    ww
for dd=1:ND
    dd

%bose hubbard parameters

%Disorder Hamiltonian
isQP=2;
beta=1.618;
Hd=MakeDisorderHam(basis,isQP,2*pi*dd/ND,beta);

W=Ws(ww)*J;
%W=2*J;
U=Us(1).*J;
Ham = Hj.*J+Hi.*U+Hd.*W;

[psiAll, PhiN, En, Cn] = ExactDiagTimeFx(psiInit,Ts,Ham);

%{
PhiNSave{ww}=PhiN;
EnSave{ww}=En;
CnSave{ww}=Cn;
%}


pdf = abs(psiAll(:,:)).^2; 
%NSamp=1E6
%[dataSmpl] = MonteCarloSmp(pdf,tind,NSamp);

nord=NSites
inds=1:nord;
%[sn2k,S2kInds,indsNKStore,corrn,gn] = dynBldGnFxFsPPMap(basis(dataSmpl,:),nord,inds,1)
for tind=1:NT
    pn=pdf(:,tind);
    %gnDat=gn;
    [sn2k,S2kInds,indsNKStore,corrn,gn] = dynBldGnFxFsPPMapPDF(pn,basis,nord,inds,0);
    
    for nn=1:length(gn)
        gnPDFStoreStd(tind,dd,nn,ww)=std(((cell2mat(values(gn{nn})))));
        gnPDFStoreMAbs(tind,dd,nn,ww)=mean(abs((cell2mat(values(gn{nn})))));
    end
end



%[tt,psiAll,Deig]=mbl_8site_function_red(J,U,W,DisordN,NT,isQP);
%{
TT{ww}=tt;
DeigS{ww}=Deig;

for dd=1:DisordN
    psiAllW{ww,dd}=psiAll{dd};
end
%}

end
end

%%

figure(11)
for ww=1:length(Ws)
for nn=1:NSites
    subplot(3,3,nn)
    semilogx(Ts./tau,mean(gnPDFStoreStd(:,:,nn,ww),2))
    hold on
end
end

%%

figure(22)
for ww=NW
for nn=1:NSites
    subplot(3,3,nn)
    semilogx(Ts./tau,mean(gnPDFStoreStd(:,:,nn,ww),2))
    hold on
end
end

%%
%{

for ww=1:NW
    tempE(ww,:)=EnSave{ww}./J;
    gapsW(ww,:)=tempE(ww,2:end)-tempE(ww,1:end-1);
    tempC(ww,:)=CnSave{ww};
    pn=tempC(ww,:).^2;
    pns(ww,:)=(pn(2:end).*pn(1:end-1)).*4;
    HL=round(HilbD/3);
    HH=round(HilbD.*2./3);
    tempGap=[gapsW(ww,HL:HH+1); gapsW(ww,HL-1:HH)];
    rspec(ww)=mean(min(tempGap,[],1)./max(tempGap,[],1),2);
    
    [ii,jj]=sort(pn,'descend');
    
end
%}
%%
%load('8site_loadme.mat');
%%
%{
nhalf=sum(basis(:,1:size(basis,2)/2),2);
ND=DisordN;


Sp=zeros(NT,length(Ws),ND);
SPM=zeros(NT,length(Ws));

NS=size(basis,2);
rhon=zeros(NT,NS+1,length(Ws),ND);
nbar=zeros(NS,NT,ND,length(Ws));

for ww=1:length(Ws)
    for dd=1:ND
        rho=abs(psiAllW{ww,dd}).^2;
        nbar(:,:,dd,ww)=(rho*basis)';
        
        for nn=1:NS+1
            indsn=find(nhalf==(nn-1));
            rhon(:,nn,ww,dd)=sum(rho(:,indsn),2);
        end
        temp=reshape(rhon(:,:,ww,dd),[NT NS+1]);
        
        Sp(:,ww,dd)=-sum(temp.*log(temp),2);
    end
    SPM=mean(Sp,3);
    nbarM=reshape(mean(nbar,3),[NS,NT,length(Ws)]);
end
%}

%%
%tj=TT{1}.*J.*2.*pi;

%cab_corr_red;

%subSystemsEP_small_half_only;

%% Correlation length section
%{

% Okay so this is looking at the old data

%NN=size(prSitesW{1},1);
NN=8;
%NS=size(prSitesW{1},2);
NS=9;
NT=size(psiAllW{1},1);
ND=size(psiAllW,2);
% find long time indices
indsLongTime=find(tt.*2.*pi.*J>100);


% create g2matrices. Need to redo simulation to keep psi
g2mats=zeros(size(basis,2),size(basis,2),size(basis,1));
nImats=zeros(size(basis,2),size(basis,2),size(basis,1));
nJmats=zeros(size(basis,2),size(basis,2),size(basis,1));
disN=zeros(NSites,size(basis,1));
g2dMats=zeros(NSites,size(basis,1));

for zz=1:size(basis,1)
    zz
    for ii=1:size(basis,2)
        for jj=ii:size(basis,2)
            if ii==jj
                g2mat(ii,jj,zz)=basis(zz,ii).^2-basis(zz,ii);

            else
                g2mat(ii,jj,zz)=basis(zz,ii).*basis(zz,jj);
            end
            nImat(ii,jj,zz)=basis(zz,ii);
            nJmat(ii,jj,zz)=basis(zz,jj);
            
            dd=abs(ii-jj);
            disN(dd+1,zz)=disN(dd+1,zz)+1;
            g2dMats(dd+1,zz)=g2mat(ii,jj,zz)+g2dMats(dd+1,zz);
            
        end
    end
end




%%
tempOp=g2dMats./disN;

for ww=1:length(Ws)
    temp=zeros(size(psiAllW{ww,1}));
    ww
    %psiAllW{ww,:
    
    
    for dd=1:DisordN
        g2dWD{ww,dd}=(abs(psiAllW{ww,dd}).^2)*tempOp';
        temp=temp+abs(psiAllW{ww,dd}).^2;
        nbarWD{ww,dd}=(abs(psiAllW{ww,dd}).^2)*basis;
        nbar2WD{ww,dd}=(abs(psiAllW{ww,dd}).^2)*(basis.^2);
    end
    
    tempG2=zeros(size(g2dWD{ww,1}));
    for dd=1:DisordN
        tempG2=g2dWD{ww,dd}+tempG2;
    end
    tempG2=tempG2./DisordN;
    g2dWDMean{ww}=tempG2;
    
    temp=temp./DisordN;
    g2dW{ww}=temp*tempOp'
    nbarW{ww}=temp*basis;
    nbar2W{ww}=temp*(basis.^2)
end

%%

avgRho=cell(length(Ws),1);
for ww=1:length(Ws)
    temp=zeros(size(psiAllW{ww,1}));
    
    for dd=1:DisordN    
        temp=temp+abs(psiAllW{ww,dd}).^2;
    end
    avgRho{ww}=temp./DisordN;
end


%%

G2DNt=cell(size(basis,2),length(Ws));
G2IJN=cell(size(basis,2),length(Ws));
G2DNW=cell(length(Ws),1);
for ww=1:length(Ws)
    ww
    temp=zeros(size(basis,2),NT);
    for ii=1:size(basis,2)
        %for ii=1:size(basis,2)-dd+1
        for jj=ii:size(basis,2)
            dd=abs(jj-ii)+1
            %jj=ii+dd-1;
            G2DNt{dd,ww}=[G2DNt{dd,ww}; reshape(g2mat(ii,jj,:),[1 size(basis,1)])*avgRho{ww}';];
            G2IJN{dd,ww}=[G2IJN{dd,ww}; (reshape(nImat(ii,jj,:),[1 size(basis,1)])*avgRho{ww}').*(reshape(nJmat(ii,jj,:),[1 size(basis,1)])*avgRho{ww}');];
        end
        %G2DNt{dd,ww}=mean(G2DNt{dd,ww},1);
        %G2IJN{dd,ww}=mean(G2IJN{dd,ww},1);
        %temp(dd,:)=G2DNt{dd,ww}./G2IJN{dd,ww};
        
        G2DNtt{dd,ww}=mean(G2DNt{dd,ww}./G2IJN{dd,ww},1);
        G2IJN{dd,ww}=mean(G2IJN{dd,ww},1);
        %temp(dd,:)=G2DNt{dd,ww}./G2IJN{dd,ww};
        temp(dd,:)=G2DNtt{dd,ww};
    end
    G2DNW{ww}=temp;
end


%%
tt=TT{1};

czz_un=zeros(NSites,NSites,length(Ws),length(tt));
czz=zeros(NSites,NSites,length(Ws),length(tt));

for ii=1:size(basis,2)
    for jj=ii:size(basis,2)
        if ii==jj
            optr=basis(:,ii).*basis(:,jj)-basis(:,ii);
            normrii=basis(:,ii);
            normrjj=basis(:,jj);            
            
        else
            optr=basis(:,ii).*basis(:,jj);
            normrii=basis(:,ii);
            normrjj=basis(:,jj);            
            
        end
        
        for ww=1:length(Ws)
            czz(ii,jj,ww,:)=avgRho{ww}*optr-(avgRho{ww}*normrii).*(avgRho{ww}*normrjj);
            
            czz_un(ii,jj,ww,:)=avgRho{ww}*optr;
            
        end
        
    end
end

%%
czz_r=zeros([size(czz,1) size(czz,4) size(czz,3)]);
for ww=1:length(Ws)
    for tti=1:length(tj)
        for dd=1:size(czz_r,1)
            czz_r(dd,tti,ww)=mean(diag( ...
                reshape(czz(:,:,ww,tti),[size(czz_r,1) size(czz_r,1)]), ...
                dd-1));
        end
        xd(tti,ww)=-sum(czz_r(:,tti,ww).*([0:7]'));
    end
end
%%
save('corrzz.mat','czz','czz_r','czz_un','Ws','tj','xd')
xdLow=xd;
tjLow=tj;
WsLow=Ws;
xd100jtLow=xd(67,:);
save('LowWXDat.mat','WsLow','tjLow','xdLow','xd100jtLow')
%}