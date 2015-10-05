clear all;
clc;
n=1000;
N=1000;
T=linspace(0,10^8,N);
%f=linspace(10^-7,10^-9,n);
f=10^-8
for i=1
    s=sin(2*pi*f(i)*T);
    c=cos(2*pi*f(i)*T);
    ns(i)=sqrt(s*s'/N);
    nc(i)=sqrt(c*c'/N);
    nm(i)=sqrt(c*s'/N);
end
hold on;
grid on;
set(gca,'XScale','log');
plot(f,ns,'b')
plot(f,nc,'g')
plot(f,nm,'r')
hold off;