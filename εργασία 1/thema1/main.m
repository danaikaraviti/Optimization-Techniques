syms x;
%f(x) = (x-2)^2 + x*log(x+3);       %f1    
%f(x) = 5^x + (2-cos(x))^2 ;         %f2      
f(x) = exp(x)*(x^3-1) + (x-1)*sin(x);  %f3

l = 0.01;
e = 0.001;
a = -1;
b = 3;
eRangeIt(f,a,b,l);
lRangeIt(f,a,b,e);
figure;
dichMethIt(f, e,  0.01, a, b, 1, [0, 0.4470, 0.7410])
dichMethIt(f, e,  0.1, a, b, 1, [0.8500, 0.3250, 0.0980])
dichMethIt(f, e,  1, a, b, 1, [0.9290, 0.6940, 0.1250])
title('variation of edges a, b');
legend('l = 0.01', '', 'l = 0.1', '', 'l = 1', '');
xlabel('k values', 'fontweight', 'bold');
hold off;