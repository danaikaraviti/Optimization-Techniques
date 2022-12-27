clear all
clc
 
% elaxistoMeGammaMetavlhto(1,-1, 0.8, 0);
 elaxistoMeGammaMetavlhto(1,-1, 1, 1);
 
 
function elaxistoMeGammaMetavlhto(x0, y0, s, flag)
    
    epsilon = 0.002;
    syms x y
    f = x^5 * exp(-x^2-y^2);
    klisi = gradient(f, [x,y]);
    
    
    k = 1;
    xList = [];             yList = [];
    xList(1) = x0;          yList(1) = y0;
    fList = [];             normKlisisList = [];        gammaList = [];
    fList(1) = subs(f, {x,y}, {xList(length(xList)), yList(length(yList))});
    normKlisisList(1) = norm(subs(klisi, {x,y}, {xList(length(xList)), yList(length(yList))}));
    
    % Armijo Constants
    a = 1/1000;             b = 0.3;         
    gammaList(1) = s;
    
    while normKlisisList(length(normKlisisList)) > epsilon
        k = k + 1;
        dk = -subs(klisi, {x,y}, {xList(length(xList)), yList(length(yList))});     % ç êëßóç ðñéí
        
        xPrin = xList(k-1);
        yPrin = yList(k-1);
        

        for mk = 0:100
            k
            mk
            gk = s * b^mk;
            xNew = xPrin + gk * dk(1);
            yNew = yPrin + gk * dk(2);
            F = matlabFunction(f);
            DF = matlabFunction(klisi);
            aristeroMelos = F(xPrin, yPrin) - F(xNew, yNew);
            grad = DF(xPrin, yPrin);
            dexioMelos = a * b^mk * s * grad' * grad;
            if aristeroMelos > dexioMelos
                gammaList(k) = gk;
                xList(k) = xNew;
                yList(k) = yNew;
                fList(k) = F(xNew, yNew);
                normKlisisList(k) = norm(subs(klisi, {x,y}, {xList(length(xList)), yList(length(yList))}));
                display('Success for this step k with mk in value of:')
                mk
                break;
            end
        end
            
            
    end
    
    xList
    yList
    fList
    normKlisisList
    k
    gammaList

    display('**********************************************************')
    xx = xList(length(xList))
    yy = yList(length(yList))
    F_xx_yy = fList(length(fList))
    NORM_KLISIS = normKlisisList(length(normKlisisList))
    if flag == 0
        plot(fList, 'ro');
        title('Red for s = 0.8, blue for s = 1');
        xlabel('Steps k');
        ylabel('Values of f');
        hold on
    else
        plot(fList, 'bx');
        title('Red for s = 0.8, blue for s = 1');
        xlabel('Steps k');
        ylabel('Values of f');
        hold on
    end
    display('**********************************************************')
end