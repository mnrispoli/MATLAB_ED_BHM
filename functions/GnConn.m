function [gnConn, gnDisConn] = GnConn(varargin)

% find connected vs disconnected correlation functions 
% monte carlo data
dataRaw=varargin{1};
inds=varargin{2};
NInd=length(inds);

% the basis produced for unity filling is actually a convenient way to find
% all unique integer partitions. Will fail for very large system sizes due
% to tryting to make the enormous basis
intParts=flipud( ...
                unique( ...
                    sort(  BasisMake(NInd,NInd),2,'descend'  ) ...
                ,'rows') ...
                );

% initialize disconnected and connected correlators
gnDisConn=GnCorrDatFx(dataRaw,inds);
gnConn=gnDisConn;

% find number of unique integer partitions. we skip the first one below
% because it is already the disconnected part+connected part
intParts=intParts(2:end,:);
NIntPrt=size(intParts,1);

for ni=1:NIntPrt
    
    gnConn=gnConn-permIndsProd(dataRaw,1,intParts(ni,:),inds);
    
end

gnDisConn=gnDisConn-gnConn;
gnDisConn=mean(prod(dataRaw(:,inds),2),1)-gnConn;

end