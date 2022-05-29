warning off all
tic

% model parameters
v0 = 1; beta = 1; g = 1; p = 0.5; q = 1-p; T = log(1/p)/g; w = 0.4;
eta = 4; d = eta*g; deff = d+g; rho0 = 0*deff; rho1 = 10*deff;
sigma1 = 3; sigma0 = 0.5*sigma1; sigma1prime = 0.8*sigma1;

% agent-based model
N = 81; maxnum = 61; itv = 4;
s0 = min(round(rho0/deff),N-1); s1 = min(round(rho1/deff),N-1); 
birth = zeros(1,2*N);
birth(1+s0) = sigma0/(sigma0+sigma1); birth(N+1+s1) = sigma1/(sigma0+sigma1);
err = 1;
while err > 1e-4
    % compute the distribution before replication
    [~,sol] = ode45(@master1,[0,w*T],birth,[],v0,beta,g,d,rho0,rho1,sigma0,sigma1,N);
    prerepli = sol(end,:);
    repli = zeros(1,4*N);
    repli(1:N) = prerepli(1:N);
    repli(3*N+1:4*N) = prerepli(N+1:2*N);
    
    % compute the distribution after replication
    [~,sol] = ode45(@master2,[0,(1-w)*T],repli,[],v0*2^w,beta,g,d,rho0,rho1,sigma0,sigma1prime,N);
    predivsn = sol(end,:);
    divsn = zeros(1,2*N);
    divsn(1:N) = predivsn(1:N)+predivsn(N+1:2*N);
    divsn(N+1:2*N) = predivsn(2*N+1:3*N)+predivsn(3*N+1:4*N);
    
    % compute the distribution at birth for the next generation
    z = exp(-2*pi*1i*(0:N-1)/N)-1;
    gen0 = zeros(1,N); gen1 = zeros(1,N);
    for i = 1:N
        gen0(i) = sum(divsn(1:N).*power(p*z(i)+1,0:N-1));
        gen1(i) = sum(divsn(N+1:2*N).*power(p*z(i)+1,0:N-1));
    end
    prob0 = real(ifft(gen0));
    prob1 = real(ifft(gen1));
    birthnext = [prob0,prob1];
    err = 1/sqrt(2)*norm(sqrt(birthnext)-sqrt(birth));
    birth = birthnext;
end

% compute the distribution before replication
[time,sol] = ode45(@master1,[0,w*T],birth,[],v0,beta,g,d,rho0,rho1,sigma0,sigma1,N);
distbirth = sol(1,1:N)+sol(1,N+1:2*N);
distrepli = sol(end,1:N)+sol(end,N+1:2*N);
len = length(time); distlin = zeros(1,N); distpop = zeros(1,N);
for i = 1:len-1
    distlin = distlin+(sol(i,1:N)+sol(i,N+1:2*N))*(time(i+1)-time(i))/T;
    distpop = distpop+2*g*exp(-g*time(i))*(sol(i,1:N)+sol(i,N+1:2*N))*(time(i+1)-time(i));
end
prerepli = sol(end,:);
repli = zeros(1,4*N);
repli(1:N) = prerepli(1:N);
repli(3*N+1:4*N) = prerepli(N+1:2*N);

% compute the distribution after replication
[time,sol] = ode45(@master2,[0,(1-w)*T],repli,[],v0*2^w,beta,g,d,rho0,rho1,sigma0,sigma1prime,N);
distdivsn = sol(end,1:N)+sol(end,N+1:2*N)+sol(end,2*N+1:3*N)+sol(end,3*N+1:4*N);
len = length(time);
for i = 1:len-1
    soltotal = sol(i,1:N)+sol(i,N+1:2*N)+sol(i,2*N+1:3*N)+sol(i,3*N+1:4*N);
    distlin = distlin+soltotal*(time(i+1)-time(i))/T;
    distpop = distpop+2^(1-w)*g*exp(-g*time(i))*soltotal*(time(i+1)-time(i));
end

% plot various distributions
plot(0:maxnum-1,distlin(1:maxnum),'r'); hold on
% plot(0:maxnum-1,distpop(1:maxnum),'r'); hold on
% plot(0:maxnum-1,distbirth(1:maxnum),'b'); hold on
% plot(0:maxnum-1,distrepli(1:maxnum),'c'); hold on
% plot(0:maxnum-1,distdivsn(1:maxnum),'c'); hold on

toc