function [cp_ghf, W_ghf] = gh_pt(m)
J = zeros(m); %please verify it
for i=1:m-1
    J(i,i)=0; % not necessary automatically also it takes zero value
    J(i,i+1)=sqrt(i/2);
    J(i+1,i)=sqrt(i/2);
end
%J is symmetric tridiagonal matrix with zero diagonal elements

[V,E]=eig(J);
%%E eigen value matrix
%%V eigen vector matrix
cp = zeros(1, m);
for i=1:m
    cp(i)=sqrt(2)*E(i,i);
end
cp_ghf=cp;
norm = zeros(1, m); W =zeros(1, m);
for j=1:m
    sum_1=0;
    for i=1:m
        sum_1=sum_1+V(i,j)^2;
    end
    norm(j)=sqrt(sum_1);
    w=V(1,j)/norm(j);
    W(j)=w^2;
end
W_ghf=W;
end