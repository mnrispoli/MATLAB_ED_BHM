% general ED bose-hubbard code
% last cleaned up 9/24/2018
% author: Matthew Rispoli

%close all
clear all
clc

addpath('functions');
%load('plot_colormaps\alex_cmap.mat');
%cmapAlex=ccs; clear ccs;

% tunneling
J=38.1911/1000; %8Er tunneling rate; defined in Hz/(2 pi)

% disorder strengths
NW=1;
Ws=linspace(0,1,NW).*J.*0 % units of J
% tupe of disorder
isQP=1; % 0 = uniform random, 1 = FFT of disorder from exp, 2 = actual QP
NB=1;
GR=(1+sqrt(5))./2;
BetaS=linspace(GR,GR,NB);


% interaction strengths
Uo=1.*2.7*J % units of J
Us=ones(size(Ws)).*Uo; %make them all the same


% system size parameters
NPart=8 % number of bosons
NSites=8 % number of sites

% num of disorders and time scans
ND=1; %disorder number
NT=151; % time steps
Ts=logspace(-1,3,NT) %actual times for ED evaluation 
tj=Ts.*2.*pi.*J;

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
psiAllSave=cell(NW,ND,NB);
nbarAllSave=cell(NW,ND,NB);
IPRSave=zeros(NW,ND,NB);

% find initial state in basis
psiInit=zeros(1,HilbD); 

%indInit = floor(NSites/2); %single-particle initial state
indInit = find(ismember(basis,ones(1,NSites),'rows')==1) % uniform density initial state
psiInit(indInit) = 1; 


%% Dynamics from ED
for ww=1:NW
    ww
for dd=1:ND
for bb=1:NB
    

%bose hubbard parameters
W=Ws(ww);
U=Us(ww);

%Disorder Hamiltonian
isQP=1;
beta=BetaS(bb);
dd.*2.*pi./ND
Hd=MakeDisorderHam(basis,isQP,dd,beta);

Ham=J.*Hj+(U/2).*Hi+W.*Hd;
try
    [psiAll, PhiN, En] = ExactDiagTimeFx(psiInit,Ts,Ham);
catch
    %incase diagonalization fails it tries again with new phase
    Hd=MakeDisorderHam(basis,isQP,2*pi*dd/ND+0.1,beta);
    Ham=J.*Hj+(U/2).*Hi+W.*Hd;
    [psiAll, PhiN, En] = ExactDiagTimeFx(psiInit,Ts,Ham);
end
psiAll(:,end)=PhiN(:,1);
%save psi's from all parameters
psiAllSave{ww,dd,bb}=psiAll;

%evaluate for on-site density, variance, and fluctuations
[nbar,nvar,nfluc]=DensityEval(psiAll,basis);

%save density stuff
nbarAllSave{ww,dd,bb}=nbar;
%inverse participation ratio
IPRSave(ww,dd,bb)=sum(nbar(:,end).^2,1);

end
end
end
%%
figure(3)

hold on
CAB=cabFx([ND NW NB NT],[NPart NSites],basis,psiAllSave);
semilogx(tj,CAB)
set(gca,'xscale','log')
grid on
ylim([0 0.25])

%%
temp=@(T) (abs(PhiN).^2)*(exp(-En./T)./(sum(exp(-En./T))))
Ttemp=logspace(-3,2,100);
vinds=1:1:size(PhiN,2);
for vv=1:length(vinds)-1
vv
%psiAllSave{1,1,1}=sqrt(temp(Ttemp(tt)));
%CAB(tt)=cabFx([ND NW NB 1],[NPart NSites],basis,psiAllSave);
psiAllSave{1,1,1}=PhiN(:,vv);
CABeig(vv)=cabFx([ND NW NB 1],[NPart NSites],basis,psiAllSave);
pn=mean(abs(PhiN(:,[vinds(vv):vinds(vv+1)])).^2,2);
%pn=temp(Ttemp(tt));
for dd=1:7
%g2(tt,dd)=sum((basis(:,dd).*basis(:,dd+1)).*pn)-sum(basis(:,dd).*pn)*sum(basis(:,dd+1).*pn);
g2Eig(vv,dd)=sum((basis(:,dd).*basis(:,dd+1)).*pn)-sum(basis(:,dd).*pn)*sum(basis(:,dd+1).*pn);
end
end
%%
for tt=1:length(Ttemp)
tt
%psiAllSave{1,1,1}=sqrt(temp(Ttemp(tt)));
%CAB(tt)=cabFx([ND NW NB 1],[NPart NSites],basis,psiAllSave);
%pn=mean(abs(PhiN(:,[vinds(vv):vinds(vv+1)])).^2,2);
pn=temp(Ttemp(tt));
for dd=1:7
g2Temp(tt,dd)=sum((basis(:,dd).*basis(:,dd+1)).*pn)-sum(basis(:,dd).*pn)*sum(basis(:,dd+1).*pn);
%g2Eig(vv,dd)=sum((basis(:,dd).*basis(:,dd+1)).*pn)-sum(basis(:,dd).*pn)*sum(basis(:,dd+1).*pn);
end
end


%%
figure()
semilogx(Ttemp./J,mean(g2Temp,2))
%%

%{
figure(1)

imagesc(BetaS,Ws./J,reshape(mean(IPRSave,2),[NW NB]))

set(gcf,'color','white')
xlabel('\beta (sites)')
ylabel('Disorder W(J)')
colormap('jet')
title('IPR: Single-Particle, Anderson Localization Example')
colorbar()
%}


%% monte carlo testing and correlator estimation

%produces probability distribution function since there certainly isn't a
%way to get the coherences back out after this process
ww=1; bb=1;
pdfAllDisAvg = DisorderAvgPsi(psiAllSave,ww,bb);

tind=NT;
NSamp=10^5;
%get monte carlo sampled basis state indices from pdf determined above
dataSmpInds = MonteCarloSmp(pdfAllDisAvg,tind,NSamp);

%create sampled dataset from basis and dataSmp indices
MCData=basis(dataSmpInds,:);
%{
for ll=1:NSites
    %{
    [gnConn,gnDisConn]=GnConn(MCData,[1:ll]);
    gnConnSave(ll)=gnConn;
    gnDisConnSave(ll)=gnDisConn;
    %}
    %{
    [gnConn]=GnConn(MCData,[1:ll]);
    gnConnSave(ll)=gnConn;
    %}
    gncon(ll)=GnConn(MCData,1:ll);
    gndis(ll)=GnCorrDatFx(MCData,1:ll)-gncon(ll);
    gn(ll)=GnCorrDatFx(MCData,1:ll);
end
%}

%gnConnSave
%gnDisConnSave

%{
%find two point correlations
g2=zeros(NSites,NSites);
for ii=1:NSites
    for jj=ii+1:NSites
        g2(ii,jj)=GnCorrDatFx(MCData,[ii,jj]);   
    end
end

%find 3 point correlations
g3=zeros(NSites,NSites,NSites);
for ii=1:NSites
    for jj=ii+1:NSites
        for kk=jj+1:NSites
            g3(ii,jj,kk)=GnCorrDatFx(MCData,[ii,jj,kk]); 
        end
    end
end

%find 4 point correlations
g4=zeros(NSites,NSites,NSites,NSites);
for ii=1:NSites
    for jj=ii+1:NSites
        for kk=jj+1:NSites
            for ll=kk+1:NSites
                g4(ii,jj,kk,ll)=GnCorrDatFx(MCData,[ii,jj,kk,ll]);
            end
        end
    end
end
%}