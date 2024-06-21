clear all;
load('results_Nt64_Nrf8_SDP.mat');
load('Rd.mat');
F = chol(Rd);
% dist_target_set = 500:500:2500;
% P_dB = 30;
% P_BS = db2pow(P_dB);
N_ite = 5e3;
M = 50;
SNRdB_set = -20:2:4;
SNR_set = db2pow(SNRdB_set);
MSE_sum = zeros(1, length(SNR_set));
MSE = zeros(1, length(SNR_set));
CRB = zeros(1, length(SNR_set));
MSE_sum_baseline = zeros(1, length(SNR_set));
MSE_baseline = zeros(1, length(SNR_set));
CRB_baseline = zeros(1, length(SNR_set));
N = 64; K = 1;
theta0 = 0;
a_psi = af(theta0, N_t);
a_dot = [0:(N_t-1)].'.*1j*pi*cos(theta0) .* a_psi;
A_dot = a_dot*a_psi.' + a_psi*a_dot.';
for n_ite = 1:N_ite
    noise_norm = randn(N,1);
    s = 2*randi(2,4,M) - 3;
    H_bu = (randn(4,N) + 1j* randn(4,N))/sqrt(2);
    for i_SNR = 1:length(SNR_set)
        SNR = SNR_set(i_SNR);
        N0 = trace(F_HAD{1}*F_HAD{1}')/SNR;
        noise = sqrt(N0)*noise_norm;
%         theta0 = 1;
        A = af(theta0,N)*af(theta0,N).';
        Y = A*F_HAD{1}*s + noise;
        y = vec(Y);
        angle_detection = [-2:0.002:2];
        P = zeros(1, length(angle_detection));
        for i_angle = 1:length(angle_detection)
            A_detec = af(angle_detection(i_angle),N)*af(angle_detection(i_angle),N).';
            vk = vec(A_detec*F_HAD{1}*s);
            P(i_angle) = abs(y'*vk + vk'*y) / (2*M*N);
        end
        [~, ind_max] = max(P);
        DOA = angle_detection(ind_max);
        MSE(i_SNR) = (DOA - theta0)^2;
        
        FHS0 = F'*H_bu'*s;
        [U0, S0, V0] = svd(FHS0);
        X0 = sqrt(M) * F* U0* eye(N,M)* V0';
        Y_baseline = A*X0 + noise;
        y_baseline = vec(Y_baseline);
        angle_detection = [-2:0.002:2];
        P = zeros(1, length(angle_detection));
        for i_angle = 1:length(angle_detection)
            A_detec = af(angle_detection(i_angle),N)*af(angle_detection(i_angle),N).';
            vk = vec(A_detec*X0);
            P(i_angle) = abs(y_baseline'*vk + vk'*y_baseline) / (2*M*N);
        end
        [~, ind_max] = max(P);
        DOA = angle_detection(ind_max);
        MSE_baseline(i_SNR) = (DOA - theta0)^2;
    end
    MSE_sum = MSE_sum + MSE;
    MSE_sum_baseline = MSE_sum_baseline + MSE_baseline;
end
for i_SNR = 1:length(SNR_set)
    SNR = SNR_set(i_SNR);
    N0 = trace(F_HAD{1}*F_HAD{1}')/SNR;
    R_N = N0*eye(N);
    CRB(i_SNR) = 1/2 *inv(real(trace(F_HAD{1}'*A_dot'*inv(R_N)*A_dot*F_HAD{1})));
    CRB_baseline(i_SNR) = 1/2 *inv(real(trace(A_dot'*inv(R_N)*A_dot*Rd)));
end
RMSE = sqrt(MSE_sum./N_ite);
RMSE_baseline = sqrt(MSE_sum_baseline./N_ite);
save('results_mle_v2.mat', 'SNRdB_set', 'RMSE', 'CRB', 'RMSE_baseline', 'CRB_baseline');