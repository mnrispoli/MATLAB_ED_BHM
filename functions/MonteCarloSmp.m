function [dataSmpl] = MonteCarloSmp(pdf,tind,NSamp)

% both in the spirit of how data is taken in the experiment and sometimes
% for the sake of not saving the entire computed wave functions for all
% parameter regimes, we can do some montecarlo sampling for fast and
% experimentally similar methods of computation

%input is assumed to be the probability distribution function in the basis
%of interest. So in this case pdf is pdf=abs(psi).^2;

%determine dimension sizes
[HilbD NT]=size(pdf);

%make sure sampled time point is valid
assert(tind<=NT && tind>0,'Time point invalid for sampling')

%initialize cumulative distribution function
cdfMsk=tril(ones(HilbD));
cdf=cdfMsk*pdf;

%to use the cdf backwards to randomly sample data points
cdfIntp=@(x) ceil(interp1(cdf(:,tind),1:HilbD,x));

%random vector for sampling
randVec=rand(1,NSamp).*(1-cdf(1,tind))+cdf(1,tind);

dataSmpl=cdfIntp(randVec);


end