function f = master2(t,x,v0,beta,g,d,rho0,rho1,sigma0,sigma1,N)
num = 4*N;
Q = zeros(num);
one = ones(num,1);
for i = 1:N-1
    Q(i+1,i) = i*d;
    Q(i+N+1,i+N) = i*d;
    Q(i+2*N+1,i+2*N) = i*d;
    Q(i+3*N+1,i+3*N) = i*d;
    Q(i,i+1) = 2*rho0*power(v0,beta)*exp(beta*g*t);
    Q(i+N,i+N+1) = (rho0+rho1)*power(v0,beta)*exp(beta*g*t);
    Q(i+2*N,i+2*N+1) = (rho0+rho1)*power(v0,beta)*exp(beta*g*t);
    Q(i+3*N,i+3*N+1) = 2*rho1*power(v0,beta)*exp(beta*g*t);
end
for i = 1:N
    Q(i,i+N) = sigma1;
    Q(i,i+2*N) = sigma1;
    Q(i+N,i) = sigma0;
    Q(i+N,i+3*N) = sigma1;
    Q(i+2*N,i) = sigma0;
    Q(i+2*N,i+3*N) = sigma1;
    Q(i+3*N,i+N) = sigma0;
    Q(i+3*N,i+2*N) = sigma0;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = Q'*x;