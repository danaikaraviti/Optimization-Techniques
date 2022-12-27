clear all
clc

xxx = 1;
yyy = -1;
elaxistoMeGammaMetavlhto(xxx, yyy, 0.8, 0);
elaxistoMeGammaMetavlhto(xxx, yyy, 1, 1);
 
 
function elaxistoMeGammaMetavlhto(x0, y0, s, flag)
    
    epsilon = 0.002;
    syms x y
    f = x^5 * exp(-x^2-y^2);
    klisi = gradient(f, [x,y]);
    KLISI = matlabFunction(klisi);
    essianos = jacobian(klisi)
    ESSIANOS = matlabFunction(essianos);        
    
    
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
        
        % ************************ Levenberg - Marquardt **********************************
        
        
        ESS = ESSIANOS(xList(k-1), yList(k-1))
        EIG = eig(ESS);
        MIN = min(EIG);
        dk = [0;0];                 
        
        if MIN > 0
            
            disp('Thetika orismenos');
            inv_essianos = inv(ESS)
            dk = - inv_essianos * KLISI(xList(k-1), yList(k-1));        
        else
            
            disp('Oxi thetika orismenos');
            mk = - MIN * 1.1
            A = ESS + mk * eye(2);
            inv_A = inv(A)
            dk = - inv_A * KLISI(xList(k-1), yList(k-1));
        end
        
        
        for mk = 0:10
            k
            mk
            gk = s * b^mk;
            xNew = xList(k-1) + gk * dk(1);
            yNew = yList(k-1) + gk * dk(2);
            F = matlabFunction(f);
            DF = matlabFunction(klisi);
            aristeroMelos = F(xList(k-1), yList(k-1)) - F(xNew, yNew);
            grad = DF(xList(k-1), yList(k-1));
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