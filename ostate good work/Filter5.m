%ECKF-RD for 1d (Nonstationary growth model)
function xhat= Filter4(n, ny, Q, R, M, phi, gam, xp3, pp3, kmax, xi1, p, z,del);
% xi1=Z.xi1; z=Z.z;
m.n=n; m.ny=ny;m.M=M; m.Q=Q; m.R=R; m.kmax=kmax; 
%----- Event triggering -------
gamt = zeros(1, m.kmax); %gamma triggering
gamt(1) = 1;
% zbar - last transmitted measurement from sensor
zbar(1) = z(1);
% zbarc - The current transmitted measurement reaches to the estimator
zbarc = zeros(m.ny, m.kmax);
zbarc(1) = z(1);
%The current transmitted measurement reaches to the estimator after delay
ybarc = zeros(m.ny, m.kmax);
ybarc(1) = z(1);
%Initial value how I will take? I think zbar is more suitable, becuase
%on sensor, we know zbar

%Initial particle drawn as x1 ~~ p(x1|x0) k=0
 x = repmat(xp3(1),1,m.M) + chol(pp3)*randn(1,m.M);
%----- Delayed measurements -----
% xprim = zeros(m.n, m.kmax);
% xprim(1) = m.xp3(1);
xhat          = zeros(m.n,m.kmax);
for k = 1:m.kmax-1
     %%  time update------------------------------
     prouviosMesparticle = gam(x);
     x = phi(x, k) + chol(m.Q)*randn(1,m.M); %  Predict the particles 
     Mesparticle = gam(x);
     %%---delay----
      ybarc(k+1) = (1-xi1(k+1))*z(k+1) + xi1(k+1)*z(k);
     %% - Measurement update -----------------------------------
    weights = ((1-p)*likelihood(ybarc( k+1),Mesparticle,m.R)+ (p)*likelihood(ybarc( k+1),prouviosMesparticle,m.R));
    if sum(weights) ==0
       disp("Bad conditioned weights for EBPF! Resetting to uniform")
    %  weights = (1/M).* ones(1, M);
      weights = likelihood(z(:, k+1),Mesparticle,40*m.R);
     end
     weights=  weights/sum(weights) ;
     xhat(k+1) = sum(repmat(weights,m.n,1).*x,2);                   % Compute state estimate
     index = sysresample(weights);                              % 4. Resample
     x = x(:,index); 
    end
  g.xp4=xhat;
  g.gamt=gamt;
end