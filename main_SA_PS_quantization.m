% 2017/08/21
% Impact of DAC quantization noise in fully-connected hybrid array.
% distortion in intended and nonintended angles are considered

% More than one beam is considered
% 3D mmW channel is used (azimuth and elevation)

% Conclusion: ???

clear;clc;
rng(2)
useQN = 1;
useEQN = 0;
noQN = 1;
MCtimes = 20;
RF_aware = 0;

% ------------ System parameters ----------------------

% opt power scaling with various parameters (from DA_required_pow_theo.m)
% case1 (U=8,16,32): -0.2984 -3.0350 -9.8275
% case2 (U=8,16,32): 3.5939 4.56 6.0829
% case3 (U=18): xx
power_scale = 10^(-3.17/10);

% Receiver planar antenna dimension
Nr_azim = 4;
Nr_elev = 2;

M = 3; % number of stream (mmW UE)
case_index = 2; % which use case in the paper
isNLOS = 1;
[SNR_origin,SINR_target] = get_target_SINR(case_index,M);

awgn_pow = 10^(-SNR_origin/10);


% Arrat geometry used in all cases (256 antenna with various output power)
Nt_range =         256;
Nt_azim_range =    32;
Nt_elev_range =    8;

if case_index == 1
    if M == 8
        Nt_azim_SA_range = 32;
        Nt_elev_SA_range = 1;
    elseif M == 16
        Nt_azim_SA_range = 16;
        Nt_elev_SA_range = 1;
    elseif M == 32
        Nt_azim_SA_range = 8;
        Nt_elev_SA_range = 1;
    end
elseif case_index == 2
    if M == 2
        Nt_azim_SA_range = 32;
        Nt_elev_SA_range = 4;
    elseif M == 4
        Nt_azim_SA_range = 32;
        Nt_elev_SA_range = 2;
    elseif M == 8
        Nt_azim_SA_range = 32;
        Nt_elev_SA_range = 1;
    end
elseif case_index == 3
    Nt_azim_SA_range = 32;
    Nt_elev_SA_range = 8;
end
    
% The optimal value can be evaluated in DA_required_pow_theo.m
% case1 (U=8,16,32): -6, -4, 4
% case2 (U=2,4,8): 30 24 22
%,case3 (U=1) 30 
BF_sigma_dB = 4;
BF_sigma = 10^(BF_sigma_dB/10);

% DSP/DAC finite number of bits
bits_range = 1:1:7;

% Make sure scheduled UEs has well seperated LOS paths
ang_grid_num = 32;
ang_grid = linspace(-pi/3,pi/3,ang_grid_num);

% -------- Channel Statistics Parameters (constant ones) -------------
cluster_num = 3; 
% LOS channel has 3 clusters and with K = 10dB
% NLOS channel has 2 clusters (LOS path has zero gain)

ray_num = 10; % number of rays in a cluster
sigma_delay_spread = 0; % No delay spread in NB model

sigma_AOA_az_spread = 10/180*pi; % 10 deg RMS spread
sigma_AOD_az_spread = 10/180*pi; % 10 deg RMS spread
sigma_AOA_el_spread = 10/180*pi; % 10 deg RMS spread
sigma_AOD_el_spread = 10/180*pi; % 10 deg RMS spread

RacianK = 13; % Ratio between LOS to the rest [dB]

% -------- zero initialization of matrices -------------
NMSE = zeros(length(Nt_azim_range), length(bits_range), MCtimes);

for MCindex = 1:MCtimes
    
    clc;fprintf('Monte Carlo Iteration %d/%d\n',MCindex,MCtimes)
    
    % -------- test signal baseband waveform  ----------
    L = 1e3;
    upsam = 5;
    for mm=1:M
        symbols = fix(L*2/upsam);   
        hmod = modem.qammod('M', 16, 'InputType', 'integer');
        hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
        hpulse = design(hdesign);
        data = randi(16,symbols,1)-1;
        data = modulate(hmod, data);
        data = upsample(data,upsam);
        temp_data = conv(data,hpulse.Numerator);
        sig = temp_data(end-L+1-1e3:end-1e3)./sqrt(temp_data(end-L+1-1e3:end-1e3)'...
            *temp_data(end-L+1-1e3:end-1e3)/L);
        sig_length = length(sig);
        sig_pow_dBm = -5;
        R = 50;
        sig_pow_V = sqrt(10^(sig_pow_dBm/10)*1e-3*R);
%         sig_pow_scale(:,mm) = sig * sig_pow_V;
        sig_pow_scale(:,mm) = sig./norm(sig)*sqrt(L);
    end
    
   % -------- Channel parameter statistics (dynamic ones) -------------
    ang_perm = randperm(ang_grid_num); % select LOS angle for UEs
    LOS_cluster_sel = ang_perm(1:M);

    % -------- Chan. Parameter Zero Initialization --------
    ray_gain = zeros(cluster_num*M,ray_num);
    ray_AOD_azim = zeros(cluster_num*M,ray_num);
    ray_AOD_elev = zeros(cluster_num*M,ray_num);
    ray_AOA_azim = zeros(cluster_num*M,ray_num);
    ray_AOA_elev = zeros(cluster_num*M,ray_num);

    % ------- Get mmW Chan. Parameters --------------
    for mm=1:M
        
        print_stat = (mm==1); % print channel statistics
        % printing channel statistic of first UE

        global_cluster_range = (mm-1)*cluster_num+1:mm*cluster_num;
        % indices for clusters that belong to the mm-th user 

        centroid_AOD_az = ang_grid(LOS_cluster_sel(mm));
        % Make sure LOS path is not overlapped

        % Randomly generate chan. parameters
        [ ray_gain(global_cluster_range,:),...
          raydelay,...
          ray_AOA_azim(global_cluster_range,:),...
          ray_AOD_azim(global_cluster_range,:),...
          ray_AOA_elev(global_cluster_range,:),...
          ray_AOD_elev(global_cluster_range,:)] =...
          get_chan_parameter_nogeo(   print_stat,...
                                      cluster_num,...
                                      ray_num,...
                                      sigma_delay_spread,...
                                      sigma_AOA_az_spread,...
                                      centroid_AOD_az,...
                                      sigma_AOD_az_spread,...
                                      sigma_AOA_el_spread,...
                                      sigma_AOD_el_spread,...
                                      isNLOS,...
                                      RacianK);
    end
    
    % -------------  for loop over Tx array size ---------------
    for tt=1:length(Nt_range)
        
        
        Nt_azim_SA = Nt_azim_SA_range(tt);
        Nt_elev_SA = Nt_elev_SA_range(tt);
        Nt_azim = Nt_azim_range(tt);
        Nt_elev = Nt_elev_range(tt);
        Nt = Nt_range(tt);
        Nt_SA = Nt_azim_SA * Nt_elev_SA;
        
        chan_raw = zeros(M, Nt);
        chan = zeros(M, Nt);
        
        % -------------  for loop for each UE ---------------
        for mm=1:M
            % Find range of parameter to generate channel
            temp_range = zeros(cluster_num,1);
            temp_range = (mm-1)*cluster_num+1:mm*cluster_num;

            MIMO_chan = get_H_MIMO_3d(      ray_gain(temp_range,:),...
                                            ray_AOD_azim(temp_range,:),...
                                            ray_AOD_elev(temp_range,:),...
                                            ray_AOA_azim(temp_range,:),...
                                            ray_AOA_elev(temp_range,:),...
                                            cluster_num,...
                                            ray_num,...
                                            Nt_azim,...
                                            Nt_elev,...
                                            Nr_azim,...
                                            Nr_elev);

            % The combining vector at Rx is from SVD and unit-magnitude scaling
            [U_mtx, Sigma_mtx, V_mtx] = svd(MIMO_chan);
            combiner = U_mtx(:,1)./abs(U_mtx(:,1)); 

            % Post combining channel is an MISO channel (a row vector) for each UE
            chan_raw(mm,:) = (combiner' * MIMO_chan)';

            % Normalize channel gain; Because each UE has the same SNR
            % Rx gain is normalized since it incorperated in link budget
            chan(mm,:) = chan_raw(mm,:)./norm(chan_raw(mm,:))*sqrt(Nt);
        end
        
        % -------------    Hybrid BF    ---------------
        for bb=1:length(bits_range)
            PS_bits_num = bits_range(bb);
            PS_scale = 2^(PS_bits_num-1);
            W_RF_PS_raw = zeros(1, Nt);
            W_RF_PS = zeros(1, Nt);
            W_RF = zeros(Nt, M);
            W_RF_orig = zeros(Nt,M);

            for mm=1:M
                sa_index = (mm-1)*Nt_SA + 1:mm*Nt_SA;
                W_RF_PS_raw(mm,:) = angle(chan(mm,:));
                W_RF_PS = (round(W_RF_PS_raw(mm,:)/pi*PS_scale+0.5)-0.5)/PS_scale*pi;
                W_RF_orig(sa_index,mm) = exp(-1j*transpose(W_RF_PS_raw(mm, sa_index)));
                W_RF(sa_index,mm) = exp(-1j*transpose(W_RF_PS(sa_index)));
            end

            for m1=1:M
                for m2=1:M
                    if RF_aware
                        chan_effct(m1, m2) = chan(m1,:) * W_RF(:,m2);
                    else
                        chan_effct(m1, m2) = chan(m1,:) * W_RF_orig(:,m2);
                    end
                end
            end

            W_BB = zeros(M,M);
            W_BB = chan_effct'*inv(chan_effct*chan_effct'+eye(M)*BF_sigma*sqrt(Nt));
            
            % ------- Simulate actual DAC Quantization -----------

            array_sig_PSQN_raw = zeros(Nt,sig_length);
            array_sig_PSQN = zeros(Nt,sig_length);
            array_sig_raw = zeros(Nt, L);
            array_sig = zeros(Nt, L);

            array_sig_raw =  W_BB * sig_pow_scale.';
            array_sig = array_sig_raw./norm(array_sig_raw)*sqrt(L*M);


            array_sig_PSQN_raw = W_RF * array_sig;
%                 array_sig_QN_raw = W_RF * array_sig;
            array_sig_PSQN = array_sig_PSQN_raw./norm(array_sig_PSQN_raw,'fro')...
                *sqrt(L)*sqrt(power_scale);



            %% Signal Distortion at Specific Angle
            for mm = 1:M
            
            rx_sig_PSQN = zeros(sig_length, 1);
            rx_sig_EQN = zeros(sig_length, 1);
            noise_pow = norm(chan(mm,:))^2/Nt * awgn_pow;
            AWGN = (randn(L,1)+1j*randn(L,1))/sqrt(2)*sqrt(noise_pow);
            
            rx_sig_PSQN = (chan(mm,:) * array_sig_PSQN).' + AWGN;
            % rx_sig is L by 1 vector (m-th rx signal at RF-chain with L samples)

            % Evaluate signal and non-signal part
            %  ------- min E||alpha * sig_org - sig_rx||^2 -----------
            sig_pow_norm = sig_pow_scale(:,mm)./norm(sig_pow_scale(:,mm))*sqrt(L);
            gain_chan(mm,bb,MCindex) = norm(chan(mm,:))^2;
            alpha_hat_PSQN = pinv(sig_pow_norm) * rx_sig_PSQN;
            gain_sig_PSQN(mm,bb,MCindex) = abs(alpha_hat_PSQN)^2;
            NMSE_QN(tt,bb,MCindex) = norm(alpha_hat_PSQN * sig_pow_norm...
                - rx_sig_PSQN)^2/norm(alpha_hat_PSQN * sig_pow_norm)^2;
            noise_pow = norm(chan(mm,:))^2/Nt * awgn_pow;
            IpN_PSQN(mm,bb,MCindex) = norm(alpha_hat_PSQN * sig_pow_norm - rx_sig_PSQN)^2/L;
            SINR_PSQN(mm,bb,MCindex) = abs(alpha_hat_PSQN)^2/(IpN_PSQN(mm,bb,MCindex));

            
            end
        end
    end
end

%% Take mean among MC realizations

SINR_mean_PSQN = zeros(M,length(bits_range));
for bb=1:length(bits_range)
    for mm=1:M
        SINR_mean_PSQN(mm,bb) = mean(squeeze(SINR_PSQN(mm,bb,:)));
    end
end


%% SINR w/ DAC QN
if useQN
    % SINR
    figure
    plot(bits_range, 10*log10(SINR_mean_PSQN'),'-o','linewidth',2);
    grid on
    xlabel('DAC Quantization Number')
    ylabel('SINR (dB)')

    % SE
    figure
    plot(bits_range, sum(log2(1+10^(-0.15/10)*SINR_mean_PSQN))','-o','linewidth',2);
%     plot(bits_range, M*log2(1+10^(0/10)*SINR_mean_QN)','-o','linewidth',2);

    grid on
    xlabel('DAC Quantization Number')
    ylabel('Spectral Efficiency (bps/Hz)')

end

