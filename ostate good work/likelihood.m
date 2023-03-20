function W = likelihood(z,zhat,R)
    M=size(zhat,2);
         d=(exp(-0.5*sum((repmat(z,1,M) - zhat).*(inv(R)*(repmat(z,1,M) - zhat)),1)));% 3. Evaluate the likelihood
         W=((d));
end