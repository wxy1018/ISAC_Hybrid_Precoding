% close all;
p0 = 30; %功率约束（dBm）
p0_w = (10^(p0/10)-3); % 功率约束（w）
N = 64;    % BS antenna number
pd = zeros(1,181);
pd(89:93) = 1;
cvx_clear
cvx_begin sdp
    variable u1
    variable t
    variable rd(N,N) hermitian semidefinite
    expressions sum_a(1,181)  x(N,1)
    minimize t
    subject to
    for i=1:181
        x = af((i-91), N); 
        sum_a(i) = ( u1*pd(i) - real(x' * rd * x) ) ^2 ;
    end
    sum(sum_a)<=t;
    trace(rd) == p0_w; 
    u1 >= 0;
cvx_end
Rd = rd;
save('Rd.mat', 'Rd');
for i = 1 : 181
    x = af((i-91), N);
    aRa(i) = real(x'*Rd*x);
end
figure;
aRa_norm = 10*log10(aRa/trace(Rd));
plot([-90:90], aRa_norm);
grid on; xlim([-90 90]);