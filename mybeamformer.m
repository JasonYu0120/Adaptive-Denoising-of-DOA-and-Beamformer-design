function y_hat = mybeamformer(matX,sfilter,ifilter,delta0,lamdb0)

N=length(matX(:,1)); %number of weight(N=10)     
L=length(matX(1,:)); %length of time(L=2000)

%initial condition
y_hat=zeros(L,1);
P=(1/delta0)*eye(N);
for m=1:L
    x=matX(:,m);
    C=[steeringvec(sfilter(m),N),...
       steeringvec(ifilter(m)-1,N),...
       steeringvec(ifilter(m),N),steeringvec(ifilter(m)+1,N),];
    g=[1;10^-4;10^-6;10^-4];
    gK=((1/lamdb0)*P*x)/(1+(1/lamdb0)*x'*P*x);
    P=(1/lamdb0)*P-(1/lamdb0)*gK*x'*P;
    w_hat=P*C/(C'*P*C)*g;
    y_hat(m)=w_hat'*x;
end
end




