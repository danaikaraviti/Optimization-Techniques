clear all
clc
 
 
syms x y
f = x^5 * exp(-x^2-y^2)
klisi = gradient(f, [x,y])
essianos = jacobian(klisi)
dk = - inv(essianos) * klisi
DK = matlabFunction(dk);

xArxiko = 1 ;       yArxiko = -1;
elaxistoMeGammaStathero(xArxiko, yArxiko, 0.2, 0);
elaxistoMeGammaStathero(xArxiko, yArxiko, 0.3, 1);
elaxistoMeGammaStathero(xArxiko, yArxiko, 0.35, 2);
 
function elaxistoMeGammaStathero(x0, y0, gamma, flag)
    
    
    epsilon = gamma/100;
    syms x y
    f = x^5 * exp(-x^2-y^2);
    klisi = gradient(f, [x,y]);
    KLISI = matlabFunction(klisi);
    essianos = jacobian(klisi)
    ESSIANOS = matlabFunction(essianos);        
    
   
    k = 1;
    xList = [];             yList = [];
    xList(1) = x0;          yList(1) = y0;
    fList = [];             normKlisisList = [];
    fList(1) = subs(f, {x,y}, {xList(length(xList)), yList(length(yList))});
    normKlisisList(1) = norm(subs(klisi, {x,y}, {xList(length(xList)), yList(length(yList))}));
 
    
    while normKlisisList(length(normKlisisList)) > epsilon
        k = k + 1
        
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
        title('Red for gamma = 0.2, blue for gamma = 0.3, green = 0.35');
        xlabel('Steps k');
        ylabel('Values of f');
        hold on
    elseif flag == 1
        plot(fList, 'bx');
        title('Red for gamma = 0.2, blue for gamma = 0.3, green = 0.35');
        xlabel('Steps k');
        ylabel('Values of f');
        hold on
    else
        plot(fList, 'g^');
        title('Red for gamma = 0.2, blue for gamma = 0.3, green = 0.35');
        xlabel('Steps k');
        ylabel('Values of f');
        hold on
    end
    display('**********************************************************')
end