function [pdfAllDisAvg] = DisorderAvgPsi(varargin)

%This exists to simply go through all the disorder realizations and average
%them. It should (physically anyway) inherently produce a probability
%distribution function since we will lose all the phase information upon
%averaging the different disorders.

% psiAllSave is the saved wave function outputs
psiAllSave = varargin{1};
% ww is the disorder depth index value
ww=varargin{2};
% bb is the incommensurate ratio index value 
bb=varargin{3};

% get all dimensions
[NW,ND,NB]=size(psiAllSave);

%initialize pdf
pdfAllDisAvg=zeros(size(psiAllSave{ww,1,bb}));

% average over disorder realizations
for dd=1:ND
    pdfAllDisAvg=pdfAllDisAvg+(abs(psiAllSave{ww,dd,bb}).^2)./ND;
end


end