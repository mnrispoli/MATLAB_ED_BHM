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
Ws=[0.1 2 6 12 18 24 30 40]
Ws=[30 30 30 30 30]
Ws=[2 4 6 8 10]
% tupe of disorder
isQP=1; % 0 = uniform random, 1 = FFT of disorder from exp, 2 = actual QP

% interaction strengths
Uo=2.7*J % units of J
Us=ones(size(Ws)).*Uo; %make them all the same
%Us=[6];

%Us=linspace(0,10,NW).*J;

% system size parameters
NPart=6 % number of bosons
NSites=6 % number of sites

% num of disorders and time scans
ND=51; %disorder number
DisordN=ND
NT=45; % time steps
tau=(1./(2.*pi.*J));
Ts=linspace(0,1000,NT) %actual times for ED evaluation 
Ts=logspace(-1.5,3,NT).*tau;
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

NW=length(Ws);

gnPDFStoreStd=zeros(NT,ND,NSites,NW);
gnPDFStoreMAbs=zeros(NT,ND,NSites,NW);
gnPDFStoreMn=zeros(NT,ND,NSites,NW);

gnPDFAvgStoreStd=zeros(NT,ND,NSites,NW);
gnPDFAvgStoreMAbs=zeros(NT,ND,NSites,NW);
gnPDFAvgStoreMn=zeros(NT,ND,NSites,NW);

PDFStore=zeros(HilbD,NT,NW);

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
U=Us(ww).*J;
Ham = Hj.*J+Hi.*U+Hd.*W;

[psiAll, PhiN, En, Cn] = ExactDiagTimeFx(psiInit,Ts,Ham);

%{
PhiNSave{ww}=PhiN;
EnSave{ww}=En;
CnSave{ww}=Cn;
%}


pdf = abs(psiAll(:,:)).^2; 

PDFStore(:,:,ww)=PDFStore(:,:,ww)+pdf./ND;
%NSamp=1E6
%[dataSmpl] = MonteCarloSmp(pdf,tind,NSamp);

nord=NSites;
inds=1:nord;
%[sn2k,S2kInds,indsNKStore,corrn,gn] = dynBldGnFxFsPPMap(basis(dataSmpl,:),nord,inds,1)
for tind=1:NT
    pn=pdf(:,tind);
    %gnDat=gn;
    [sn2k,S2kInds,indsNKStore,corrn,gn] = dynBldGnFxFsPPMapPDF(pn,basis,nord,inds,0);
    
    for nn=1:length(gn)
        gnPDFStoreStd(tind,dd,nn,ww)=std(((cell2mat(values(gn{nn})))));
        gnPDFStoreMAbs(tind,dd,nn,ww)=mean(abs((cell2mat(values(gn{nn})))));
        gnPDFStoreMn(tind,dd,nn,ww)=mean(((cell2mat(values(gn{nn})))));
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

for tind=1:NT
    pn=reshape(PDFStore(:,tind,ww),size(pn));
    %gnDat=gn;
    [sn2k,S2kInds,indsNKStore,corrn,gn] = dynBldGnFxFsPPMapPDF(pn,basis,nord,inds,0);
    
    for nn=1:length(gn)
        gnPDFAvgStoreStd(tind,dd,nn,ww)=std(((cell2mat(values(gn{nn})))));
        gnPDFAvgStoreMAbs(tind,dd,nn,ww)=mean(abs((cell2mat(values(gn{nn})))));
        gnPDFAvgStoreMn(tind,dd,nn,ww)=mean(((cell2mat(values(gn{nn})))));
    end
end

end

%%

figure(11)
hold off
for ww=1:length(Ws)
for nn=1:NSites
    subplot(2,4,nn)
    loglog(Ts./tau,mean(gnPDFStoreStd(:,:,nn,ww),2))
    hold on
    title(sprintf('W/J=30; N=%i',nn))
    grid on
    xlim([Ts(1)./tau,Ts(end)./tau])
    xlabel('time(\tau)')
    ylabel('Std[G^{(N)}]')
    Ws
    
end
end
set(gcf,'color','white')
legend('U/J=0.001','U/J=0.01','U/J=0.1','U/J=1','U/J=2')



figure(12)
hold off
for ww=1:length(Ws)
for nn=1:NSites
    subplot(2,4,nn)
    loglog(Ts./tau,mean(gnPDFStoreMAbs(:,:,nn,ww),2))
    hold on
    title(sprintf('W/J=30; N=%i',nn))
    grid on
    xlim([Ts(1)./tau,Ts(end)./tau])
    xlabel('time(\tau)')
    ylabel('Mean|G^{(N)}|')
    Ws
    
end
end
set(gcf,'color','white')
legend('U/J=0.001','U/J=0.01','U/J=0.1','U/J=1','U/J=2')



figure(13)
hold off
for ww=1:length(Ws)
for nn=1:NSites
    subplot(2,4,nn)
    loglog(Ts./tau,abs(mean(gnPDFStoreMn(:,:,nn,ww),2)))
    hold on
    title(sprintf('W/J=30; N=%i',nn))
    grid on
    xlim([Ts(1)./tau,Ts(end)./tau])
    xlabel('time(\tau)')
    ylabel('Mean(G^{(N)})')
    Ws
    
end
end
set(gcf,'color','white')
legend('U/J=0.001','U/J=0.01','U/J=0.1','U/J=1','U/J=2')

%%


figure(22)
hold off
for ww=1:length(Ws)
for nn=1:NSites
    subplot(2,4,nn)
    loglog(Ts./tau,mean(gnPDFAvgStoreStd(:,:,nn,ww),2))
    hold on
    title(sprintf('W/J=30; N=%i',nn))
    grid on
    xlim([Ts(1)./tau,Ts(end)./tau])
    xlabel('time(\tau)')
    ylabel('Std[G^{(N)}]')
    Ws
    
end
end
set(gcf,'color','white')
legend('U/J=0.001','U/J=0.01','U/J=0.1','U/J=1','U/J=2')



figure(23)
hold off
for ww=1:length(Ws)
for nn=1:NSites
    subplot(2,4,nn)
    loglog(Ts./tau,mean(gnPDFAvgStoreMAbs(:,:,nn,ww),2))
    hold on
    title(sprintf('W/J=30; N=%i',nn))
    grid on
    xlim([Ts(1)./tau,Ts(end)./tau])
    xlabel('time(\tau)')
    ylabel('Mean|G^{(N)}|')
    Ws
    
end
end
set(gcf,'color','white')
legend('U/J=0.001','U/J=0.01','U/J=0.1','U/J=1','U/J=2')



figure(24)
hold off
for ww=1:length(Ws)
for nn=1:NSites
    subplot(2,4,nn)
    loglog(Ts./tau,abs(mean(gnPDFAvgStoreMn(:,:,nn,ww),2)))
    hold on
    title(sprintf('W/J=30; N=%i',nn))
    grid on
    xlim([Ts(1)./tau,Ts(end)./tau])
    xlabel('time(\tau)')
    ylabel('Mean(G^{(N)})')
    Ws
    
end
end
set(gcf,'color','white')
legend('U/J=0.001','U/J=0.01','U/J=0.1','U/J=1','U/J=2')

