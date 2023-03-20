function gm= sysmodel(m,phi, gam,p)
    %% Truth generation
    x = zeros(m.n, m.kmax); z = zeros(m.ny, m.kmax);
    y = zeros(m.ny, m.kmax);
    x(1) = normrnd(m.xp3(1), sqrt(m.pp3));
    z(1) = gam(x(1)) + normrnd(0, sqrt(m.R));
    y(1) = z(1);
    itr = 0;

    xi1 = rand(1, m.kmax) <= p;
    for k = 1:m.kmax-1
        x(k+1) = phi(x(k), k) + normrnd(0, sqrt(m.Q));
        z(k+1) = gam(x(k+1)) + normrnd(0, sqrt(m.R));
        y(k+1) = (1-xi1(k+1))*z(k+1) + xi1(k+1)*z(k);
    end
 gm.x=x;
 gm.xi1=xi1;
 gm.z=z;
 gm.y=y;
end