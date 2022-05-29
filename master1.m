function f = master1(t,x,v0,beta,g,d,rho0,rho1,sigma0,sigma1,N)
num = 2*N;
Q = zeros(num);
one = ones(num,1);
for i = 1:N-1
    Q(i+1,i) = i*d;
    Q(i+N+1,i+N) = i*d;
    Q(i,i+1) = rho0*power(v0,beta)*exp(beta*g*t);
    Q(i+N,i+N+1) = rho1*power(v0,beta)*exp(beta*g*t);
end
for i = 1:N
    Q(i,i+N) = sigma1;
    Q(i+N,i) = sigma0;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = Q'*x;