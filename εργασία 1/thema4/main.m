syms x;
%  f(x) = (x-2)^2 + x*log(x+3);       %f1    
% f(x) = 5^x + (2-cos(x))^2 ;         %f2      
f(x) = exp(x)*(x^3-1) + (x-1)*sin(x);  %f3
a = -1;
b = 3;
lRange(f,a,b);
figure;
derDechMeth(f, 0.01, a, b, 1, [0, 0.4470, 0.7410])
derDechMeth(f, 0.1, a, b, 1, [0.8500, 0.3250, 0.0980])
derDechMeth(f, 1, a, b, 1, [0.9290, 0.6940, 0.1250])
title('variation of edges a, b');
legend('a', 'b');
legend('l = 0.01', '', 'l = 0.1', '', 'l = 1', '');
xlabel('k values', 'fontweight', 'bold');
hold off;
