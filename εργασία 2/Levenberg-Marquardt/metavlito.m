clear all
clc

xArxiko = -1;       yArxiko = 1;
elaxistoMeGammaMetavlhto(xArxiko, yArxiko, 0.1, 0.3, 0);
elaxistoMeGammaMetavlhto(xArxiko, yArxiko, 0.3, 0.5, 1);
 
 
function elaxistoMeGammaMetavlhto(x0, y0, akro1, akro2, flag)
    
    
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
    gammaList(1) = 0;
    
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
        
        
        gamma = internalOptimization(f, xList(k-1), yList(k-1), dk, akro1, akro2);
        
       
        xList(k) = xList(k-1) + gamma * dk(1);
        yList(k) = yList(k-1) + gamma * dk(2);
        fList(k) = subs(f, {x,y}, {xList(length(xList)), yList(length(yList))});
        normKlisisList(k) = norm(subs(klisi, {x,y}, {xList(length(xList)), yList(length(yList))}));
        gammaList(k) = gamma;
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
        title('Red for 0.1 <= gamma <= 0.3, blue for 0.3 <= gamma <= 0.5');
        xlabel('Steps k');
        ylabel('Values of f');
        hold on
    else
        plot(fList, 'bx');
        title('Red for 0.1 <= gamma <= 0.3, blue for 0.3 <= gamma <= 0.5');
        xlabel('Steps k');
        ylabel('Values of f');
        hold on
    end
    display('**********************************************************')
end
 
 
 
 
 
 
 
 
function gamma = internalOptimization(f, xPrin, yPrin, dk, akro1, akro2)
    syms x y G
    X = xPrin + G * dk(1);
    Y = yPrin + G * dk(2);
   
    F = subs(f, {x,y}, {X, Y});
    DF = diff(F, 'G');
    l = 0.005;
    counter = 0;
    while akro2 - akro1 > l
        counter = counter + 1;
        kentro = (akro1 + akro2) / 2;
        paragwgos = subs(DF, kentro);
        if paragwgos > 0
            akro2 = kentro;
        else
            akro1 = kentro;
        end
     end
     
     gamma = (akro2 + akro1) / 2;
end