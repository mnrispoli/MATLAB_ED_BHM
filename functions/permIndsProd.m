function [prodOut] = permIndsProd(dataIn,prodIn,knum,inds)
    
    if size(inds,2)~=0
        nkInds=nchoosek(inds,knum(1));
        
        for nk=1:size(nkInds,1)            
            indsOut=setdiff(inds,nkInds(nk,:));
            prodOut=GnCorrDatFx(dataIn,nkInds(nk,:))*permIndsProd(dataIn,prodIn,knum(2:end),indsOut);
        end
        
    else
        prodOut=prodIn;
    end

end