syms x y
f = (x^2)/3 + 3*(y^2);
fsurf(f, 'red');
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Plot of function: f = (x^2)/3 + 3*(y^2);');