
% num_Nt = 1:1:3;
% for num = num_Nt
%     N_t = 2^(num+3); %16/32/64
%     N_RF = 8;
%     M = N_t/N_RF;
%     noise_dB_user = -100;
%     sigma2_user = db2pow(noise_dB_user);
%     psi_user = [-75, -50, 30, 45];
%     K = length(psi_user);   % number of users
%     gamma_th_dB = 20;
%     gamma_th = db2pow(gamma_th_dB);
%     gamma_th_dB_set = 10:1:19;
%     CRB = zeros(size(gamma_th_dB_set));
%     n_NL = 4;   % mmWave NLOS path number
%     psi_nlos = zeros(n_NL, K);
%     alpha_nlos = zeros(n_NL, K);
%     h_user = zeros(N_t, K);
%     for k = 1: K
%         psi_nlos(:,k) = psi_user(k) + rand(n_NL,1)*4-2;
%         alpha_nlos(1:n_NL,k) = 1./db2pow(5*rand(n_NL,1)+5); % nlos path lossï¼?5-10 dB
%         h_user(:,k) = af(psi_user(k),N_t);
%         for i = 1:n_NL
%             h_user(:,k) = h_user(:,k) + alpha_nlos(i,k)*af(psi_nlos(i,k), N_t);
%         end
%     end
%     rho0_dB = -60; % channel gain at the reference distance
%     rho0 = db2pow(rho0_dB);
%     dist_comm = 5e1; % distance between BS and users
%     H_total{num} = h_user*sqrt(rho0*dist_comm^(-2))/sqrt(sigma2_user); %normalized by noise_power
% end
% save('channel_v3.mat', 'H_total');

load channel_v1.mat;
H{1} = H_total{1};
load channel_v2.mat;
H{2} = H_total{2};
load channel_v3.mat;
H{3} = H_total{3};
H_total = H;
save('channel.mat', 'H_total');