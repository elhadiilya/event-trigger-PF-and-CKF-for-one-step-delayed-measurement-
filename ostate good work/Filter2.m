%% CKF-RD for 1d (Nonstationary growth model)

function xp2 = Filter2(n, ny, Q, R, cp, W, col, phi, gam, y, xp2, pp2, kmax, p)

xprim = zeros(n, kmax);
xprim(1) = xp2(1);
ppric = zeros(n, kmax);
ppric(1) = pp2;

for k = 1:kmax-1
    %% Time update
    s = sqrt(pp2);
    kipred = zeros(n, col); kiupd = zeros(n, col);
    for j =1:col
        kipred(j) = s*cp(j) + xp2(k);
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

    %% -- Measurement update
    sp = sqrt(ppric(k+1));
    kiipred = zeros(n, col); zpred = zeros(ny, col);
    for j = 1:col
        kiipred(j) = sp*cp(j) + xprim(k+1);
        zpred(j) = gam(kiipred(j));
    end

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

    Pyy = (1-p)*Pzz + p*Pzzz + p*(1-p)*(zhat - zzhat)^2; %checked, it's correct

    Pxy = (1-p)*Pxz + p*Pxzz;

    %-----------------------------------------

    K = Pxy/Pyy;
    xp2(k+1) = xpri + K*(y(k+1) - yhat);
    pp2 = ppri - K*Pyy*K';

end