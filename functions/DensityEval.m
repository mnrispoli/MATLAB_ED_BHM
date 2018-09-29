function [nbar,nvar,nfluc] = DensityEval(psiAll,basis);


nbar=basis'*abs(psiAll).^2;
nvar=(basis.^2)'*abs(psiAll).^2;
nfluc=nvar-nbar;


end