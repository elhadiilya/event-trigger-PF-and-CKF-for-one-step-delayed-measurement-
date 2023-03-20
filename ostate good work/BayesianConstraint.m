function W  = BayesianConstraint(z,zparticle,R,del)
     M=size(zparticle,2);W=(1/M).* ones(1, M);;
     for i=1:20
      yt = zparticle-repmat(z,1,M); % Here y_k is not available, so we compute y_k-E[y_k|I_{k-1}] using predicted particles
        y_diff = yt+(chol(R)*randn(2,M)); % y_k-E[y_k|I_{k-1}] for each of the particle x_k^i stacked columnwise
        y_diff =(sum((y_diff).^2,1).^(1/2));
        y_diff(find(y_diff==0)) = 1;
        y_diff(find(abs(y_diff)>del)) = 0;
        y_diff(find(y_diff~=0)) = 1;
   W=W+(y_diff); 
     end
   W=W/20;
%    dd=sum(y_diff/(M));
   LH =sum(log(y_diff))/M^2;