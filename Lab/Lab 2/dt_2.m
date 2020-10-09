function y=dt_2(x)
y= zeros(1,length(x));

for n=1:length(x);
    if n<2
        y(n) = x(n);
    elseif n<3
        y(n) = x(n) - 2*y(n-1);
    else
        y(n) = x(n) - 2*y(n-1) - 2*y(n-2)
    end
end