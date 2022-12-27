syms x y
f = x^5 * exp(-x^2-y^2);
fsurf(f, 'red');
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Plot of function: f = x^5 * exp(-x^2-y^2)');