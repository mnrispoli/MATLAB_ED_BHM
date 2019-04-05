% general ED bose-hubbard code
% last cleaned up 9/24/2018
% author: Matthew Rispoli

addpath('functions')

close all
clear all
clc
NSites=2:8

gnSave=zeros(NSites(end),NSites(end));
rhoSave=zeros(NSites(end),NSites(end));

for ns=1:length(NSites)
    ns

    basis=BasisMake(NSites(ns),NSites(ns));
    nf=mean(mean(basis.^2,1)-mean(basis,1).^2,2);
    pn=1./size(basis,1);
    nord=NSites(ns);
    inds=1:nord;

    [sn2k,S2kInds,indsNKStore,corrn,gn] = dynBldGnFxFsPPMapPDF(pn,basis,nord,inds,1);

    for nn=1:nord
        gnSave(ns,nn)=mean(cell2mat(values(gn{nn})));
        rhoSave(ns,nn)=mean(cell2mat(values(gn{nn}))./(sqrt(nf).^nn));
    end
end
%%
figure(1)
subplot(1,2,1)
imagesc(2:12,1:12,abs(gnSave(1:end-1,:)'))
axis('square')
colorbar()

subplot(1,2,2)
imagesc(2:12,1:12,abs(rhoSave(1:end-1,:)'))
axis('square')
colorbar()


