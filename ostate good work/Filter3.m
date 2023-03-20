%ECKF-RD for 1d (Nonstationary growth model)
function [xp3, gamt] = Filter3(n, ny, Q, R, cp, W, col, phi, gam, y, xp3, pp3, kmax, xi1, p, z,del)

%----- Event triggering -------
gamt = zeros(1, kmax); %gamma triggering
gamt(1) = 1;
% zbar - last transmitted measurement from sensor
zbar(1) = z(1);
% zbarc - The current transmitted measurement reaches to the estimator
zbarc = zeros(ny, kmax);
zbarc(1) = z(1);
%The current transmitted measurement reaches to the estimator after delay
ybarc = zeros(ny, kmax);
ybarc(1) = z(1);
%Initial value how I will take? I think zbar is more suitable, becuase
%on sensor, we know zbar

alp1 = 0.3;  % 0.02 In Li Li paper
alp2 = 0.35; % 0.02;
alp3 = 0.3;
alp4 = 0.35;
alp5 = 0.3;
alp6 = 0.35;
alp7 = 0.3;


%----- Delayed measurements -----
xprim = zeros(n, kmax);
xprim(1) = xp3(1);
ppric = zeros(n, kmax);
ppric(1) = pp3;

for k = 1:kmax-1
    s = sqrt(pp3);
    kipred = zeros(1, col); kiupd = zeros(1, col);
    for j =1:col
        kipred(j) = s*cp(j) + xp3(k);
        kiupd(j) = phi(kipred(j), k);
    end

    xpri = 0;
    for j = 1:col
        xpri = xpri + W(j)*kiupd(j);
    end
    xprim(k+1) = xpri;

    ppri = 0;
    for j = 1:col
        ppri = ppri + W(j)*(kiupd(j) - xpri)^2;
    end
    ppri = ppri + Q;
    ppric(k+1) = ppri;


    %% - Measurement update
    sp = sqrt(ppric(k+1));
    kiipred = zeros(n, col); zpred = zeros(ny, col);
    for j = 1:col
        kiipred(j) = sp*cp(j) + xprim(k+1);
        zpred(j) = gam(kiipred(j));
    end
    %-------------------
    B1 = W.*zpred;

    B = B1*cp';

    %----------------------
    zhat = 0;
    for j = 1:col
        zhat = zhat + W(j)*zpred(j);
    end
    Pzz = 0; Pxz = 0;
    for j = 1:col
        Pzz = Pzz + W(j)*(zpred(j) - zhat)^2;
        Pxz = Pxz + W(j)*(kiipred(j) - xprim(k+1))*(zpred(j) - zhat);
    end
    Pzz = Pzz+R;


    %--- one step previous measurement
    ssp = sqrt(ppric(k));
    kiiipred = zeros(n, col); zzpred = zeros(ny, col);
    for j = 1:col
        kiiipred(j) = ssp*cp(j) + xprim(k);
        zzpred(j) = gam(kiiipred(j));
    end

    %-------------------
    B2 = W.*zzpred;

    BB = B2*cp';
    %----------------------

    zzhat = 0;
    for j = 1:col
        zzhat = zzhat + W(j)*zzpred(j);
    end

    Pzzz = 0; Pxzz = 0;
    for j = 1:col
        Pzzz = Pzzz + W(j)*(zzpred(j) - zzhat)^2;
        Pxzz = Pxzz + W(j)*(kiiipred(j) - xprim(k))*(zzpred(j) - zzhat);
    end
    Pzzz = Pzzz + R;

    %--- Combining both the measurements ----


    yhat = (1-p)*zhat + p*zzhat;

    Pyy = (1-p)*Pzz + p*Pzzz + p*(1-p)*(zhat - zzhat)^2; %need to be checked

    Pxy = (1-p)*Pxz + p*Pxzz;


    %-----------------------------------------

    if (z(k+1) - zbar)'*(z(k+1) - zbar) > del
        gamt(k+1) = 1;
        zbar = z(k+1);
    else
        gamt(k+1) = 0;
    end

    zbarc(k+1) =  gamt(k+1)*z(k+1) + (1-gamt(k+1))*zbar;


    ybarc(k+1) = (1-xi1(k+1))*zbarc(k+1) + xi1(k+1)*zbarc(k);

    L = Pxy/Pyy; % Normal Kalman gain (measurement is being sent) gam1 = 1, triggering case

    H1 = B/sp; %H_{k/k-1}
    H2 = BB/ssp; %H_{k-1/k-2}




    %Basically multiplier in the Kalman gain, K
    K1 = (1+alp1+alp2)*(1-p)^2*H1*ppric(k+1)*H1' + (1 + alp3)*(1-p)*R + (1+alp4)*p*(1-p)*zhat^2 - p*(1-p)*zhat*zzhat...
          + (1 + 1/alp1 + alp5)*p*H2*ppric(k)*H2' + (1 + alp6)*p*R - p*(1-p)*zzhat*zhat' + p*(1-p)*(1+alp7)*zzhat^2 ...
          + (1 + 1/alp2 +1/alp3 + 1/alp4 + 1/alp5 + 1/alp6 + 1/alp7)*del*eye(ny); 

    K = (1+alp1+alp2)*(1-p)*ppric(k+1)*H1'/K1;  %Kalman gain when gam1 = 0

    xp3(k+1) = xpri + gamt(k+1)*L*(y(k+1) - yhat) + (1-gamt(k+1))*K*(ybarc(k+1) - yhat);

    pp1 = ppri - L*Pyy*L';
    pp2 = (1 + alp1 + alp2)*(eye(n) - (1-p)*K*H1)*ppric(k+1)*(eye(n) - (1-p)*K*H1)' + (1 + alp3)*(1-p)*K*R*K' + p*(1-p)*(1+alp4)*K*zhat^2*K' ...
           -p*(1-p)*K*zhat*zzhat*K' + (1 + 1/alp1 + alp5)*p*K*H2*ppric(k)*H2'*K' + (1 + alp6)*p*K*R*K' - p*(1-p)*K*zzhat*zhat*K' + p*(1-p)*(1+alp7)*K*zzhat^2*K' ...
           + (1 + 1/alp2 +1/alp3 + 1/alp4 + 1/alp5 + 1/alp6 + 1/alp7)*del*K*eye(ny)*K'; 

    pp3 = gamt(k+1)*pp1 + (1-gamt(k+1))*pp2;

end