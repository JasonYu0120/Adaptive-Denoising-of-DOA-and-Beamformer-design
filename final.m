clear all
clc
load('ASP_Final_Data.mat') 

%% DOA estimation by EMA

L=length(theta_s_noisy(1,:));
W=1:20;
sigma1=1;
sigma2=2;
ker1=9*sigma1;% ker1*2+1=window size
ker2=9*sigma2;
sfilter=zeros(1,L);
ifilter=zeros(1,L);

aa=0.06;
fcut=1/(2*pi)*acos(1-(aa^2)/2*(1-aa));
% Recursive moving average filter
ssum1=theta_s_noisy(1,1);
for i=1:L
    ssum1=(1-aa)*ssum1+aa*(theta_s_noisy(1,i));
    sfilter(1,i)=ssum1;
end

bb=0.01;
ssum2=theta_i_noisy(1,1);
for i=1:L
    ssum2=(1-bb)*ssum2+bb*(theta_i_noisy(1,i));
    ifilter(1,i)=ssum2;
end

figure
plot(sfilter)
title('filtered DOA')
hold on
plot(ifilter)
legend('filtered theta_s noisy','filtered theta_i noisy')
xlabel('time')
ylabel('DOA in degree')
savefig('ASP_Final_DOA.fig')
theta_s_hat=sfilter;
theta_i_hat=ifilter;
%% My beamformer (LCMV)

delta0=0.1;
lamdb0=0.99;
s_t_hat=mybeamformer(matX,sfilter,ifilter,delta0,lamdb0);
s_t_hat=s_t_hat.';
figure
subplot(2,1,1)
plot(real(s_t_hat))
title('real part of estimated s(t)')
xlabel('time (n)')
ylabel('amplitude')
subplot(2,1,2)
plot(imag(s_t_hat))
title('imaginary part of estimated s(t)')
xlabel('time (n)')
ylabel('amplitude')
savefig('ASP_Final_Source.fig')

save('Estimated.mat',"theta_i_hat","theta_s_hat",'s_t_hat')