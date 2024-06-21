clear all; close all;
global Pi_bold Pi_ki lambda gamma_th K
opt_FBB = 1;    %1: SDR;  2: SOCP
max_iteration = 10;
epsilon1 = 1e-3;
epsilon2 = 1;
%% communication related parameters
P_dB = 30;
P_BS = db2pow(P_dB);
N_t = 64;
N_RF = 32;
M = N_t/N_RF;
noise_dB_user = -100;
sigma2_user = db2pow(noise_dB_user);
psi_user = [-75, -50, 30, 45];
K = length(psi_user);   % number of users
% gamma_th_dB_set = [0:1:8];
gamma_th_dB_set = [5:2:15];
CRB = zeros(size(gamma_th_dB_set));
n_NL = 4;   % mmWave NLOS path number
psi_nlos = zeros(n_NL, K);
alpha_nlos = zeros(n_NL, K);
h_user = zeros(N_t, K);
for k = 1: K
    psi_nlos(:,k) = psi_user(k) + rand(n_NL,1)*4-2;
    alpha_nlos(1:n_NL,k) = 1./db2pow(5*rand(n_NL,1)+5); % nlos path lossï¼?5-10 dB
    h_user(:,k) = af(psi_user(k),N_t);
    for i = 1:n_NL
        h_user(:,k) = h_user(:,k) + alpha_nlos(i,k)*af(psi_nlos(i,k), N_t);
    end
end
rho0_dB = -60; % channel gain at the reference distance
rho0 = db2pow(rho0_dB);
dist_comm = 5e1; % distance between BS and users
H = h_user*sqrt(rho0*dist_comm^(-2))/sqrt(sigma2_user); %normalized by noise_power
% H = h_user;
load channel.mat;
H = H_total{log2(N_t)-3};

%% radar related parameters
psi_0 = 0;
noise_dB_radar = -100;
sigma2_radar = db2pow(noise_dB_radar);
psi_0 = 0;  %DOA of target
psi_clutter = [-40, 30, 60];    %DOA of clutters
n_clutter = length(psi_clutter);  %clutter number
dist_radar = 5e2; % distance between BS and target/clutter (assuming target and clutters in the same range bin)
xi0 = 1*sqrt(rho0*dist_radar^(-4)); % target RCS & path loss
xi_clutter = 0.1*sqrt(rho0*dist_radar^(-2))*ones(1, n_clutter); % clutter RCS
R_N = eye(N_t,N_t);

% F_BB = randn(N_RF, K)+ 1j*randn(N_RF, K);
F_RF = zeros(N_t, N_RF);
for i = 1:N_RF
    i_user = mod(i,K);
    if i_user == 0
        i_user = K;
    end
    F_RF(:,i) = af(psi_user(i_user),N_t);
end
F_RF = (randn(N_t, N_RF) + 1j*randn(N_t, N_RF));
F_RF = F_RF./abs(F_RF);
% load channel.mat;
    F_RF = zeros(N_t, N_RF);
    for m = 1:N_RF
        F_RF((M*m-M+1):(M*m), m) = ones(M,1);
    end

a_psi = af(psi_0, N_t);
a_dot = [0:(N_t-1)].'.*1j*pi*cos(psi_0) .* a_psi;
A_dot = a_dot*a_psi.' + a_psi*a_dot.';
A0 = af(psi_0, N_t)*af(psi_0, N_t).';
% R_in = R_N;
% for i = 1 : n_clutter
%     A = af(psi_clutter(i), N_t)*af(psi_clutter(i), N_t).';
%     R_in = R_in + xi_clutter(i)^2/sigma2_radar*A*F_RF*F_BB*F_BB'*F_RF'*A';
% end
% Xi = xi0^2/sigma2_radar*A0'*inv(R_in)*A0;

Z = A_dot'*inv(R_N)*A_dot;

for k = 1:K
    H_tilde{k} = F_RF'*H(:,k) *H(:,k)'*F_RF;
end
x = ones(N_t,1);
% load F_RF_32_8.mat;

for n_gamma = 1:length(gamma_th_dB_set)
    tic
    gamma_th_dB = gamma_th_dB_set(n_gamma);
    gamma_th = db2pow(gamma_th_dB);
    F_RF = zeros(N_t, N_RF);
    for m = 1:N_RF
        F_RF((M*m-M+1):(M*m), m) = ones(M,1);
    end
    delta = 1; n_ite = 1;
    objective_old = 0;
%     if opt_FBB == 1
%         load('F_RF/F_RF_32_16.mat');
%     end
    F_BB = randn(N_RF, K)+ 1j*randn(N_RF, K);
    F_BB = F_BB/norm(F_BB, 'fro')*sqrt(P_BS/M);
    F_BB_tilde = zeros(N_RF, K);
    Phi = F_RF'*A_dot'*inv(R_N)*A_dot*F_RF;
    convergence = 1/real(1/trace(F_BB'*Phi*F_BB)/ xi0^2 *sigma2_radar);
%     for n_ite = 1:4
    while (real(delta) > epsilon1) && (n_ite <= max_iteration)
        n_ite = n_ite + 1;
        %% optimization of digital beamformer
        Phi = F_RF'*A_dot'*inv(R_N)*A_dot*F_RF;
        for k = 1:K
            H_tilde{k} = F_RF'*H(:,k) *H(:,k)'*F_RF;
        end
        if opt_FBB == 1
            cvx_clear
            cvx_solver sedumi
            cvx_begin 
            cvx_precision high
            variable T(N_RF, N_RF, K) hermitian semidefinite;
            expression v(K, K+1);
            maximize real(trace(Phi*sum(T,3))/1e8)
            
            subject to
                for i_u = 1:K
                    v(i_u,1) = 0;
                    for j_u = 1:K
                        v(i_u,j_u+1) = v(i_u,j_u) + trace(H_tilde{i_u}*T(:,:,j_u));
                    end
                    real((1+gamma_th)*trace(H_tilde{i_u}*T(:,:,i_u)) - gamma_th*v(i_u,K+1)) >= gamma_th;
            %         real((1+gamma_th)*trace(R_bu(:,:,i_u)*T(:,:,i_u)) - gamma_th*v(i_u,K+1)) >= gamma_th * (sigma2_user);
                end
                trace(sum(T,3)) <= N_RF/N_t *P_BS;
            cvx_end
            
            for k = 1:K
                [eigenvector, eigenvalue] = eig(T(:,:,k));
                F_BB(:,k) = sqrt(eigenvalue(end,end))*eigenvector(:,end);
                gamma(k) = trace(H_tilde{k}*T(:,:,k))/(trace(H_tilde{k}*sum(T,3)) - trace(H_tilde{k}*T(:,:,k)) + 1);
            %     gamma1(k) = F_BB(:,k)'*H_tilde{k}*F_BB(:,k)/(trace(H_tilde{k}*sum(T,3)) - trace(H_tilde{k}*T(:,:,k)) + sigma2_user);
            end
            eigenvalue;
        else
            h_tilde = F_RF'*H;
            obj1 = 0; obj2 = 1e5;
%             num_ite = 0;
%             for nn = 1:2
            while (obj2 - obj1)/obj1 >= 1e-2
                obj1 = obj2;
                cvx_clear
                cvx_begin
                variable F_BB(N_RF, K) complex
                variable t(K)
                expression u(K,K-1)
                maximize sum(t)
            
                subject to
                    norm(vec(F_BB),2) <= sqrt(N_RF/N_t *P_BS);
                    for k = 1:K
                        n_u = 0;
                        for i = 1:K
                            if i ~= k
                                n_u = n_u + 1;
                                u(k,n_u) = 2*h_tilde(:,k)'*F_BB(:,i);
                            end
                        end
                        norm([u(k,:), 2, real(h_tilde(:,k)'*F_BB(:,k))/sqrt(gamma_th)-1]) <= real(h_tilde(:,k)'*F_BB(:,k))/sqrt(gamma_th)+1;
                        1e8*t(k) <= 2*real(F_BB_tilde(:,k)'*Phi*F_BB(:,k)) - real(F_BB_tilde(:,k)'*Phi*F_BB_tilde(:,k));
                    end
                cvx_end
                obj2 = trace(F_BB'*Phi*F_BB);
                F_BB_tilde = F_BB;
            end
            for k = 1:K
                gamma(k) = trace(h_tilde(:,k)*h_tilde(:,k)'*F_BB(:,k)*F_BB(:,k)')/(trace(h_tilde(:,k)*h_tilde(:,k)'*F_BB*F_BB') - trace(h_tilde(:,k)*h_tilde(:,k)'*F_BB(:,k)*F_BB(:,k)') + 1);
            %     gamma1(k) = F_BB(:,k)'*H_tilde{k}*F_BB(:,k)/(trace(H_tilde{k}*sum(T,3)) - trace(H_tilde{k}*T(:,:,k)) + sigma2_user);
            end
        end
        objective_d = trace(F_BB'*Phi*F_BB);
        
        %% optimization of analog beamformer
    %     F_RF = zeros(N_t, N_RF);
    %     for m = 1:N_RF
    %         F_RF((2*m-1), m) = 1;
    %         F_RF((2*m), m) = -1;
    %     end
        B = F_BB*F_BB';
        Pi_bold = zeros(N_t, N_t);
        for m = 1:N_RF
            for n = 1:N_RF
                BZ = kron(B(n,m), Z);
                Pi_bold(((m-1)*M+1):(m*M), ((n-1)*M+1):(n*M)) = BZ(((m-1)*M+1):(m*M), ((n-1)*M+1):(n*M));
            end
        end
    %     d = ones(N_t, 1); d(2:2:end) = -1;
        for k = 1:K
            B_k{k} = F_BB(:,k)*F_BB(:,k)'; H_k{k} = H(:,k)*H(:,k)';
        end
        for k = 1:K
            for i = 1:K
                for m = 1:N_RF
                    for n = 1:N_RF
                        BH = kron(B_k{i}(n,m), H_k{k});
                        Pi_ki{k,i}(((m-1)*M+1):(m*M), ((n-1)*M+1):(n*M)) = BH(((m-1)*M+1):(m*M), ((n-1)*M+1):(n*M));
                    end
                end
            end
        end
        lambda = ones(K, 1);
        while(1)
            manifold = complexcirclefactory(N_t);
            problem.M = manifold;
            problem.cost = @mycost;
            problem.egrad = @myegrad;
        %     options.verbosity = 1;
            [x, xcost, info] = steepestdescent(problem, x);
            for k = 1:K
                t = 0;
                for i = 1:K
                    if i ~= k
                        t = t + x'*Pi_ki{k,i}*x;
                    end
                end
                phi(k) = real(t*gamma_th+gamma_th*1 - x'*Pi_ki{k,k}*x);
                if phi(k) > epsilon2
                    lambda(k) = lambda(k)*10;
                end
            end
            if max(phi) < epsilon2
                break;
            end
        end
        objective_a = x'*Pi_bold*x;
        % k = 3; i = 3;
        % trace(F_RF'*H_k{k}*F_RF*B_k{i})
        % d'*Pi_ki{k,i}*d
        F_RF = zeros(N_t, N_RF);
        for m = 1:N_RF
            F_RF((M*m-M+1):(M*m), m) = x((M*m-M+1):(M*m));
        end
        delta = (objective_a - objective_old) / objective_old;
        objective_old = objective_a;
        if N_t == N_RF
            break;
        end
        convergence = [convergence, 1/real(1/objective_a/ xi0^2 *sigma2_radar)];
    end
    F_RF = zeros(N_t, N_RF);
    for m = 1:N_RF
        F_RF((M*m-M+1):(M*m), m) = x((M*m-M+1):(M*m));
    end
    Phi = F_RF'*A_dot'*inv(R_N)*A_dot*F_RF;
    CRB(n_gamma) = real(1/trace(F_BB'*Phi*F_BB) / xi0^2 * sigma2_radar);
    F_HAD{n_gamma}= F_RF*F_BB;
    complexity(n_gamma) = toc;
end
plot(1:length(convergence), convergence, '-bo', 'LineWidth', 1.2);
save('results_complexity_Nrf32_SDP.mat', 'gamma_th_dB_set', 'complexity', 'P_dB', 'N_t', 'N_RF', 'opt_FBB', 'F_HAD')

function [f_value] = mycost(x)
global Pi_bold Pi_ki lambda gamma_th K
f0 = -x'*Pi_bold*x;
for k = 1:K
    t = 0;
    for i = 1:K
        if i ~= k
            t = t + real(x'*Pi_ki{k,i}*x);
        end
    end
    if (t*gamma_th+gamma_th*1 - real(x'*Pi_ki{k,k}*x)) >0
        f0 = f0+ lambda(k)*(t*gamma_th+gamma_th*1 - x'*Pi_ki{k,k}*x)^2;
    end
end
f_value = real(f0)/1e8;
end

function [g] = myegrad(x)
global Pi_bold Pi_ki lambda gamma_th K
g0 = -2*Pi_bold*x;
for k = 1:K
    t = 0; Pi_sum = zeros(size(Pi_bold));
    for i = 1:K
        if i ~= k
            t = t + x'*Pi_ki{k,i}*x;
            Pi_sum = Pi_sum + Pi_ki{k,i};
        end
    end
    if (t*gamma_th+gamma_th*1 - x'*Pi_ki{k,k}*x) >0
        g0 = g0+ lambda(k)*4*(-x'*Pi_ki{k,k}*x + t*gamma_th + gamma_th*1) * (gamma_th*Pi_sum - Pi_ki{k,k})*x;
    end
end
g = g0/1e8;
end