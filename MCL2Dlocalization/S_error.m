%% error distance
    function distance_error = S_error(a,b,Sd)
    [r,~]=size(b);
    distance_error=zeros(1,r);
    for i= 1:r
        m=dist(a,b(i,:));
        distance_error(1,i)=m+randn(1,1)*Sd;  
        
%         distance_error(1,i)=m*(1+randn(1,1)*0.01);   
    end
    end
