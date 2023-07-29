function [dS] = Channel_rank(N,array,r,lambda)
for n1=1:N
    for n2=1:N
        vc1z=array(n1);
        vc1y=0;
        vc2z=array(n2);
        vc2y=r;
        dd=sqrt((vc2z-vc1z)^2+(vc2y-vc1y)^2);
        HN(n1,n2)=sqrt(1/N/N)*exp(-1i*2*pi/lambda*dd);
    end
end
[~,S,~]=svd(HN);
dS=diag(S);
end

