clear all;
load('results_Nt64_Nrf8_SDP.mat');
load('Rd.mat');
% dist_target_set = 500:500:2500;
% P_dB = 30;
% P_BS = db2pow(P_dB);
N_ite = 1e4;
SNRdB_set = -20:2:-4;
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
    for i_SNR = 1:length(SNR_set)
        SNR = SNR_set(i_SNR);
        N0 = trace(F_HAD{1}*F_HAD{1}')/SNR;
        noise = sqrt(N0)*noise_norm;
%         theta0 = 1;
        A = af(theta0,N)*af(theta0,N).';
        R_Z = A*F_HAD{1}*F_HAD{1}'*A' + noise*noise';
        R_Z1 = A*Rd*A' + noise*noise';
%         R_Z = af(theta0/180*pi,N)*real(af(theta0/180*pi,N).'*RX*conj(af(theta0/180*pi,N)))*af(theta0/180*pi,N)'/trace(RX);
        [V_R, D_R] = eig(R_Z);
        noiseSub = V_R(:, K+1:N);
        theta_est = -5:0.01:5;
        a_est = zeros(N, length(theta_est));
        res = zeros(length(theta_est), 1);
        for i_est = 1:length(theta_est)
            a_est(:,i_est) = af(theta_est(i_est), N);
            res(i_est) = 1/(norm(a_est(:, i_est)'*noiseSub).^2);
        end
        [resSorted, orgInd] = sort(res, 'descend');
        DOA = theta_est(orgInd(1:K));
        aa = diff(res);
        aa = sign(aa);
        aa = diff(aa);
        bb = find(aa==-2)+1;
        k = 0; DOA = [];
        for mm = 1:length(orgInd)
            if ismember(orgInd(mm),bb)
                DOA= [DOA, theta_est(orgInd(mm))];
                k = k + 1;
            end
            if k == K
                break;
            end
        end
        MSE(i_SNR) = (DOA - theta0)^2;
        
        [V_R, D_R] = eig(R_Z1);
        noiseSub = V_R(:, K+1:N);
        theta_est = -5:0.01:5;
        a_est = zeros(N, length(theta_est));
        res = zeros(length(theta_est), 1);
        for i_est = 1:length(theta_est)
            a_est(:,i_est) = af(theta_est(i_est), N);
            res(i_est) = 1/(norm(a_est(:, i_est)'*noiseSub).^2);
        end
        [resSorted, orgInd] = sort(res, 'descend');
        DOA = theta_est(orgInd(1:K));
        aa = diff(res);
        aa = sign(aa);
        aa = diff(aa);
        bb = find(aa==-2)+1;
        k = 0; DOA = [];
        for mm = 1:length(orgInd)
            if ismember(orgInd(mm),bb)
                DOA= [DOA, theta_est(orgInd(mm))];
                k = k + 1;
            end
            if k == K
                break;
            end
        end
        if isempty(DOA)
            DOA = 0;
        end
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
save('results_rmse_v2.mat', 'SNRdB_set', 'RMSE', 'CRB', 'RMSE_baseline', 'CRB_baseline');