function [jterms] = FindTunnelingTerms(NSites,isPB)

jterm=zeros([1 NSites]);
jterm(1:2)=[1,-1];
jterm=jterm(toeplitz(1:numel(jterm),[1 numel(jterm):-1:2]));

if isPB == 0
    %makes open boundary condition
    jterms=jterm(2:end,:);
else
    %makes periodic boundary conditions
    jterms=jterm;
end

%make hermitian terms
jterms=vertcat(jterms,fliplr(jterms));


end