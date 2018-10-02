function [gn] = GnCorrDatFx(varargin)

% This is to calculate the correlation of non-overlapping indices in the
% data structure. This is, as of now, just for the density-density...
% correlations and so overlapping indices could be fixed with commutation
% rules. Fix this stuff later if it comes up

%define varagin
%sampled data (in fock basis, but generally whatever the basis states are)
dataStruct=varargin{1};
%indices -- this determines the "order" of the correlator 
inds=varargin{2};

%find the dimensions of sampled monte carlo data
[NSmp NSites]=size(dataStruct);

%on-site mean values to find differences from
meanVals=mean(dataStruct(:,inds),1);
deltaValsInds=dataStruct(:,inds)-ones(NSmp,1)*meanVals;
gn=mean(prod(deltaValsInds,2),1);

end