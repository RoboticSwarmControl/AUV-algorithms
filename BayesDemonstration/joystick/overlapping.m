%% range search 
function inmatrix=overlapping(X,c,r,mode)
    if nargin<4
    mode=2;
    end
    if mode==2
        idx = (X(:,1)-c(1)).^2+(X(:,2)-c(2)).^2 <= r^2;
    %[idx,dist]=range2(c,r,X);
    else
    [idx,dist]=range1(c,r,X);
    end
    inmatrix = X(idx,:);
end




    function [idx,dist]=range1(c,r,X)
    [nPoints,nVariables]=size(X);
    s=zeros(nPoints,1);
    for d=1:nVariables
    x=abs(X(:,d)-c(d));
    id=x>s;
    s(id)=x(id);
    end
    fidx=s<r;
    idx=find(fidx);
    dist=s(fidx);
    end

    function [idx,dist]=range2(c,r,X)
    nVariables=numel(c);
    r2=r*r;
    s=0;
    for d=1:nVariables
    s=s+(X(:,d)-c(d)).^2;
    end
    fidx=s<r2;
    idx=find(fidx);
    dist=s(fidx);
    end
        
