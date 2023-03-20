%Standard CKF for nonstationary growth model
function xp1 = Filter1(n, ny, Q, R, cp, W, col, phi, gam, y, xp1, pp1, kmax)

for k = 1:kmax-1
    s = sqrt(pp1);
    kipred = zeros(n, col); kiupd = zeros(n, col);
    for j =1:col
        kipred(j) = s*cp(j) + xp1(k);
        kiupd(j) = phi(kipred(j), k);
    end

    xpri = 0;
    for j = 1:col
        xpri = xpri + W(j)*kiupd(j);
    end
    ppri = 0;
    for j = 1:col
        ppri = ppri + W(j)*(kiupd(j) - xpri)^2;
    end
    ppri = ppri + Q;



    sp = sqrt(ppri);
    kiipred = zeros(n, col); zpred = zeros(ny, col);
    for j = 1:col
        kiipred(j) = sp*cp(j) + xpri;
        zpred(j) = gam(kiipred(j));
    end
    zhat = 0;
    for j = 1:col
        zhat = zhat + W(j)*zpred(j);
    end
    Pzz = 0; Pxz = 0;
    for j = 1:col
        Pzz = Pzz + W(j)*(zpred(j) - zhat)^2;
        Pxz = Pxz + W(j)*(kiipred(j) - xpri)*(zpred(j) - zhat);
    end
    Pzz = Pzz+R;

    K1 = Pxz/Pzz;
    xp1(k+1) = xpri + K1*(y(k+1) - zhat); %Because measurement is randomly delayed
    pp1 = ppri - K1*Pzz*K1';

end