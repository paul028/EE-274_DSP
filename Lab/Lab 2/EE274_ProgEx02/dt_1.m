function y=dt_1(x)
y= zeros(1,length(x)); % zero initialization
for n=1:length(x);
    if n<2
        y(n) = x(n);
    else
        y(n) = 0.5*x(n) + 0.5*x(n-1)
    end
end
end