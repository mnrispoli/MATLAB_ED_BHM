function [Hd] = MakeDisorderHam(varargin)

% arg 1 is the basis for the disorder hamiltonian
% arg 2 is the QP disorder boolean that determines the type of disorder
%       isQP = 0 -  random, 1 - exp disorder, 2 - quasiperiodic, o.w.
%       nothing else
% arg 3 either offset in experimental samples or phi offset for qp
% arg 4 quasiperiodicity

basis=varargin{1};
isQP=varargin{2};

NSites=size(basis,2);
HilbD=size(basis,1);


if isQP == 0
    
    % make random disorder vector
    Rvec=2.*rand(NSites,1)';
    Hd=diag( ...
        sum((ones(HilbD,1)*Rvec).*basis,2) ...
        );
    
elseif isQP == 1
    
    % load experimental disorder
    load('Hamiltonians/ExpDisorder.mat');
    
    % normalize disorder
    siteSave=savedData{3};
    siteSave=siteSave./max(max(siteSave));
    
    % create vector to multiply disorder by
    PhiOff=varargin{3};
    
    Rvec=siteSave(8-(round(NSites/2)-1):9+(round(NSites/2)-1),round(PhiOff))';
    Hd=diag( ...
        sum((ones(HilbD,1)*Rvec).*basis,2) ...
        );
            
elseif isQP == 2
    
    beta=varargin{4};
    
    % make quasiperiodic disorder vector
    PhiOff=varargin{3};

    Rvec=cos((beta*2*pi).*[0:NSites-1]+PhiOff);
    Hd=diag( ...
        sum((ones(HilbD,1)*Rvec).*basis,2) ...
        );
    
else
    
    Hd=zeros(size(basis,1));
    
end