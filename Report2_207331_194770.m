%% Digital Data Transmission - Patricia Martinez 207331 - Sofia Martinez 194770
%% Informacion mutua para una constelacion continua

clear all
close all

for i=1:4
nMc=1000;
sigma2 = 1;
enV = linspace(0,20,100);
%enV = 10*log10(enV); %Energia en dB
mInf = zeros(size(enV));
C = zeros(size(enV));
SL = zeros(size(enV));
Z = zeros(1,nMc);

    switch i
        case 1 %QAM_ continous
            for iE =1:length(enV)
                E = enV(iE);
                delta = sqrt(3*E);
                %Entrada tiene valores 0 o Delta (parte real e imaginaria)
                %Condicion para que los puntos se generen dentro de un cuadrado
                X_0 = -delta + (delta+delta)*rand(1,nMc); 
                X_1 = -delta + (delta+delta)*rand(1,nMc); 
                X = X_0+X_1*1j;
                [mInf_c(iE), C(iE), SL(iE)]=continous_IM(X, nMc,sigma2, Z, delta, iE); 
                mInf(iE)= mInf_c(iE);
            end
        case 2 %4-PAM
            X =  generate(4);
            [mInf, enV]= discrete_MI(X, nMc); 
        case 3 %4-PSK
            X = [-1-1j -1+1j 1-1j 1+1j]/sqrt(2);
            [mInf, enV]= discrete_MI(X, nMc); 
        case 4 %9-QAM
            X = [-3-3*1j -1-3*1j 1-3*1j 3-3*1j -3-1j -1+1*1j 1+1*1j 3+1*1j -3+3*1j -1+3*1j 1+3*1j 3+3*1j];
            X_energy = X * energy(X);
            [mInf, enV]= discrete_MI(X, nMc); 
    end



if i == 1 
    figure()
    B=C/SL; %bandwidth
    plot(enV,B*SL);
    xlabel('SNR');
    ylabel('Mutual information bits/channel use');
end

hold on;
plot(enV,mInf);
end
legend('Capacidad','Modulacion continua','4-PAM','4-PSK','9-QAM', 'Location','southeast');
hold off;
figure()
scatter(X_0,X_1,10,'ro','filled')
xlabel('Real value of constellation X');
ylabel('Imaginary value of constellation X');
axis([-delta-2 delta+2 -delta-2 delta+2]);
%% Aproximacion de la informacion mutua
HQ = zeros(1,length(enV));
I_aprox = zeros(1,length(enV));
HW = log2(pi*exp(1)*sigma2);
for i=1:length(enV)  
 E = enV(i);
 delta = (E*3);
 HQ(i) = (delta*2)^2;
 I_aprox(i) = log2(HQ(i)*HW);
end

figure()
plot(enV,I_aprox,'b');
hold on
plot(enV,mInf_c,'r');
xlabel('Average signal energy');
ylabel('Mutual information bits/channel use');
legend('Información mutua aproximada','Información mutua Monte Carlo','Location','southeast');

%% FUNCIONES PARA GENERAR CONSTELACIONES Y INFORMACION  MUTUA

function [mInf, enV]= discrete_MI(X, nMc)
%Funcion que calcula la informacion mutua de una constelaicon discreta
    sigma2 =  1;
    M = length(X);
    Q = ones(size(X))/M;
    enV = linspace(0,20,100);
    mInf =zeros(size(enV));

    for iE =1:length(enV)
        E = enV(iE); 
        I=zeros(1,nMc);
        for iZ=1:nMc
            z = normrnd(0, sqrt(sigma2/2),1,nMc)+ 1j*normrnd(0,sqrt(sigma2/2),1,nMc);
            for iX=1:M
                 x = X(iX)*sqrt(E);
                 y =x+z;
                 Py=0;
                 for iM = 1:M
                   Py = Py + Q(iM) * (1/(pi*sigma2))*exp((-abs(y- sqrt(E)*X(iM)).^2)/sigma2); 
                 end
                 I(iZ) = I(iZ) + Q(iX) * log2(((1/(pi * sigma2))*exp(-abs(z).^2/sigma2)) / Py);
            end
        end
        mInf(iE) = mean(I);
    end
end

function [mInf, C, SL]=continous_IM(X, nMc,sigma2, Z, delta, iE)
%Funcion que calcula la informacion mutua de una constelaicon continua
    P = 0;
    I=zeros(1,nMc); 
    for iZ=1:nMc
        Z(iZ) = sqrt(sigma2/2)*(randn+ 1j*randn);
        Y = delta*X(iZ)+Z(iZ);
        
        X_0_bar = -delta + (delta+delta)*rand(1,nMc); 
        X_1_bar = -delta + (delta+delta)*rand(1,nMc); 
        X_bar = X_0_bar+X_1_bar*1j;
                
        for iX= 1:nMc
            P =  P + exp((-abs(Y-delta*X_bar(iX)).^2)/sigma2);   
        end
        P = P/nMc;
        I(iZ) =(log2(exp(-abs(Z(iZ)/sigma2)^2)/P));
    end
    
    C = max(I);             %Capacidad
    SL= log2(1+ iE/sigma2); %Capacidad * Bandwith
    mInf = mean(I);
end

 

function DIS= generate(M)
%Funcion que genera la constelacion PAM y QAM
    mitad= M/2; 
    DIS_POS=1:1:mitad; 
    DIS_NEG=-mitad:1:-1; 
    %PAM
        DIS=[DIS_NEG DIS_POS]; 
        DIS=DIS/sqrt(energy(DIS)); 
end 

function Es=energy(DIS)
%Funcion que calcula la energia de una constelacion
    Es=0; 
    for i= 1:length(DIS)
        Es= Es+(abs(DIS(i))^2*(1/length(DIS))); 
    end
end
