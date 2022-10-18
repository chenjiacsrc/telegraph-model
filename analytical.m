warning off all
tic

% model parameters
v0 = 1; beta = 1; g = 1; p = 0.5; q = 1-p; T = log(1/p)/g; w = 0.4;
eta = 4; d = eta*g; deff = d+g; rho0 = 0*deff; rho1 = 10*deff;
sigma1 = 3; sigma0 = 0.5*sigma1; sigma1prime = 0.8*sigma1;

% FSP algorithm
N = 81; maxnum = 61; t = 0.6*T;
s0 = min(round(rho0/deff),N-1); s1 = min(round(rho1/deff),N-1); 
birth = zeros(1,2*N);
birth(1+s0) = sigma0/(sigma0+sigma1); birth(N+1+s1) = sigma1/(sigma0+sigma1);
if t <= w*T
    [~,sol] = ode45(@master1,[0,t],birth,[],v0,beta,g,d,rho0,rho1,sigma0,sigma1,N); 
    prob = sol(end,:);
    dist = prob(1:N)+prob(N+1:2*N);
else
    [~,sol] = ode45(@master1,[0,w*T],birth,[],v0,beta,g,d,rho0,rho1,sigma0,sigma1,N);
    prerepli = sol(end,:);
    repli = zeros(1,4*N);
    repli(1:N) = prerepli(1:N);
    repli(3*N+1:4*N) = prerepli(N+1:2*N);
    [~,sol] = ode45(@master2,[0,t-w*T],repli,[],v0*2^w,beta,g,d,rho0,rho1,sigma0,sigma1prime,N);
    prob = sol(end,:);
    dist = prob(1:N)+prob(N+1:2*N)+prob(2*N+1:3*N)+prob(3*N+1:4*N);
end
plot(0:maxnum-1,dist(1:maxnum),'b'); hold on

% derived parameters
itv = 3;
r = sigma0+sigma1-beta*g;
a = sigma1/(d+beta*g);
b = (sigma0+sigma1)/(d+beta*g);
u = (rho1-rho0)*power(v0,beta)/(d+beta*g);
v = rho0*power(v0,beta)/(d+beta*g);
rprime = sigma0+sigma1prime-beta*g;
aprime = sigma1prime/(d+beta*g);
bprime = (sigma0+sigma1prime)/(d+beta*g);
uprime = (rho1-rho0)*power(v0,beta)*2^(beta*w)/(d+beta*g);
vprime = rho0*power(v0,beta)*2^(beta*w)/(d+beta*g);
z = exp(-2*pi*1i*(0:N-1)/N)-1;
prop = exp((v*exp(beta*g*t)-(u+v)*exp(-d*t))*z);

% analytical solution
L0 = (hypergeom(1+a-b,1-b,u*exp(-d*t)*z).*hypergeom(a,b,u*exp(beta*g*t)*z)+...
    a*u*z/b/(b-1).*exp(-r*t).*hypergeom(1+a,1+b,u*exp(-d*t)*z).*hypergeom(1+a-b,2-b,u*exp(beta*g*t)*z)).*prop;
L1 = (hypergeom(a-b,1-b,u*exp(-d*t)*z).*hypergeom(a,b,u*exp(beta*g*t)*z)-...
    (b-a)*u*z/b/(b-1).*exp(-r*t).*hypergeom(a,1+b,u*exp(-d*t)*z).*hypergeom(1+a-b,2-b,u*exp(beta*g*t)*z)).*prop;
tprime = t-w*T;
propprime = exp((vprime*exp(beta*g*tprime)-(uprime+vprime)*exp(-d*tprime))*z);
L0prime = (hypergeom(1+aprime-bprime,1-bprime,uprime*exp(-d*tprime)*z).*hypergeom(aprime,bprime,uprime*exp(beta*g*tprime)*z)+...
    aprime*uprime*z/bprime/(bprime-1).*exp(-rprime*tprime).*hypergeom(1+aprime,1+bprime,uprime*exp(-d*tprime)*z).*...
    hypergeom(1+aprime-bprime,2-bprime,uprime*exp(beta*g*tprime)*z)).*propprime;
L1prime = (hypergeom(aprime-bprime,1-bprime,uprime*exp(-d*tprime)*z).*hypergeom(aprime,bprime,uprime*exp(beta*g*tprime)*z)-...
    (bprime-aprime)*uprime*z/bprime/(bprime-1).*exp(-rprime*tprime).*hypergeom(aprime,1+bprime,uprime*exp(-d*tprime)*z).*...
    hypergeom(1+aprime-bprime,2-bprime,uprime*exp(beta*g*tprime)*z)).*propprime;
F0 = zeros(1,N); F1 = zeros(1,N);
if birth(1) == 1
    F0 = ones(1,N);
elseif birth(N+1) == 1
    F1 = ones(1,N);
else
    for i=1:N
        F0(i) = sum(birth(1:N).*power(exp(-d*t)*z(i)+1,0:N-1));
        F1(i) = sum(birth(N+1:2*N).*power(exp(-d*t)*z(i)+1,0:N-1));
    end
end
if t <= w*T
    gen = F0.*L0+F1.*L1;
else
    z = exp(-d*tprime)*z;
    K00 = (b-a)/b*(hypergeom(1+a-b,1-b,u*exp(-d*w*T)*z).*hypergeom(a,1+b,u*exp(beta*g*w*T)*z)+...
        a/(b-a)*exp(-(r+beta*g)*w*T)*hypergeom(1+a,1+b,u*exp(-d*w*T)*z).*hypergeom(a-b,1-b,u*exp(beta*g*w*T)*z)).*prop;
    K01 = (b-a)/b*(hypergeom(a-b,1-b,u*exp(-d*w*T)*z).*hypergeom(a,1+b,u*exp(beta*g*w*T)*z)-...
        exp(-(r+beta*g)*w*T)*hypergeom(a,1+b,u*exp(-d*w*T)*z).*hypergeom(a-b,1-b,u*exp(beta*g*w*T)*z)).*prop;
    K10 = a/b*(hypergeom(1+a-b,1-b,u*exp(-d*w*T)*z).*hypergeom(1+a,1+b,u*exp(beta*g*w*T)*z)-...
        exp(-(r+beta*g)*w*T)*hypergeom(1+a,1+b,u*exp(-d*w*T)*z).*hypergeom(1+a-b,1-b,u*exp(beta*g*w*T)*z)).*prop;
    K11 = a/b*(hypergeom(a-b,1-b,u*exp(-d*w*T)*z).*hypergeom(1+a,1+b,u*exp(beta*g*w*T)*z)+...
        (b-a)/a*exp(-(r+beta*g)*w*T)*hypergeom(a,1+b,u*exp(-d*w*T)*z).*hypergeom(1+a-b,1-b,u*exp(beta*g*w*T)*z)).*prop;
    F0rep = F0.*K00+F1.*K01;
    F1rep = F0.*K10+F1.*K11;
    gen = F0rep.*power(L0prime,2)+F1rep.*power(L1prime,2);
end
dist = ifft(gen);
plot(0:itv:maxnum-1,dist(1:itv:maxnum),'ro'); hold on

toc