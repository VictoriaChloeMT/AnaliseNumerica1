Accuracy Test:
n=100; B=cholhilb(n); A=hilb(n); max(max(B*B'-A))
peed Test:
% tic;B=cholhilb(1000);toc
% tic;for k=1:100,B=cholhilb(k);end;toc
%
B=zeros(n);
b=1;
for i=1:n,
    c=b;
    d=b;
    for j=1:i,
        c=c*(i+j-1);
        B(i,j)=sqrt(2*j-1)*b*b/(c*d);
        d=d/max(1,(i-j));
    end
    b=b*i;
end
