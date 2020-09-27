%{
For f
o=1/50 and N=200, write a program to quantize the signal x(n), using
truncation, to 64, 128 and 256 quantization levels. Plot x(n), xq
(n) and e(n)
and compute SQNR
%}

N= 200;
fo=1/50;
n=[0:N]';
x=sin(2*pi*fo*n)
plot(n,x)

%Q=64,128,256 
Q=64;
q=2/(Q+1);
xq=q*(floor(x/q));
plot(xq)
xlabel('n')
ylabel('xq(n) - x(n)')


plot(xq-x);
 axis([0 200 -9e-3 0])
xlabel('n')
ylabel('xq(n) - x(n)')
Px=sum((x).^2)/N;
Pq=sum((xq-x).^2)/N;
SQNR=10*log10(Px/Pq)