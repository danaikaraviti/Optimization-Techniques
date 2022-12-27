clear all;
clc


 elaxistoMeGammaStathero(1,-1, 0.2, 0);
 elaxistoMeGammaStathero(1,-1, 0.5, 1);

function elaxistoMeGammaStathero(x0, y0, gamma, flag)
    
    epsilon = gamma / 100;
    syms x y
    f = x^5 * exp(-x^2-y^2);
    klisi = gradient(f, [x,y]);
    
    k = 1;
    xList = [];             yList = [];
    xList(1) = x0;          yList(1) = y0;
    fList = [];             normKlisisList = [];
    fList(1) = subs(f, {x,y}, {xList(length(xList)), yList(length(yList))});
    normKlisisList(1) = norm(subs(klisi, {x,y}, {xList(length(xList)), yList(length(yList))}));

    
    while normKlisisList(length(normKlisisList)) > epsilon
        k = k + 1;
        dk = -subs(klisi, {x,y}, {xList(length(xList)), yList(length(yList))});
       
        xList(k) = xList(k-1) + gamma * dk(1);
        yList(k) = yList(k-1) + gamma * dk(2);
        fList(k) = subs(f, {x,y}, {xList(length(xList)), yList(length(yList))});
        normKlisisList(k) = norm(subs(klisi, {x,y}, {xList(length(xList)), yList(length(yList))}));
    end
    
    
    xList
    yList
    fList
    normKlisisList
    k
    
    display('**********************************************************')

    xx = xList(length(xList))
    yy = yList(length(yList))
    F_xx_yy = fList(length(fList))
    NORM_KLISIS = normKlisisList(length(normKlisisList))
    if flag == 0
        plot(fList, 'ro');
        title('Red for gamma = 0.2, blue for gamma = 0.5');
        xlabel('Steps k');
        ylabel('Values of f');
        hold on
    else
        plot(fList, 'bx');
        title('Red for gamma = 0.2, blue for gamma = 0.5');
        xlabel('Steps k');
        ylabel('Values of f');
        hold on
    end
    
    display('**********************************************************')
end