%Standard CKF for nonstationary growth model
function [er gamt] = Filter(n, ny, Q, R, cp, W, col,M, phi, gam, y, xp, pp, kmax,xi1,z, x,del,p)

 %% --- CKF-RD -------------
    xp1 = Filter2(n, ny, Q, R, cp, W, col, phi, gam, y,xp, pp, kmax, p); 
    e1 = (xp1 - x).^2;
      
    %% --- ECKF-RD -------------
    [xp2, gamt1] = Filter3(n, ny, Q, R, cp, W, col, phi, gam, y, xp, pp, kmax, xi1, p, z,del); 
     e2 = (xp2 - x).^2;

      %% ----- PF-RD -----
    xp3  = Filter5(n, ny, Q, R, M, phi, gam, xp, pp, kmax, xi1, p, z,0); 
    e3 = (xp3 - x).^2;

    %% ----- EPF-RD -----
    [xp4 gamt2] = Filter4(n, ny, Q, R, M, phi, gam, xp, pp, kmax, xi1, p, z,del); 
    e4 = (xp4 - x).^2;
    er=[e1; e2 ;e3 ;e4];
    gamt=[gamt1 gamt2]';
end