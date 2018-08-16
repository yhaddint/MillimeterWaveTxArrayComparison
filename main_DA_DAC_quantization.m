% 2017/12/18
% Impact of DAC quantization noise in digital array.
% distortion in intended and nonintended angles are considered

% More than one beam is considered
% 3D mmW channel is used (azimuth and elevation)

% Adjust channel model to be consistent (and compare) with hybrid one

% 256 Tx array serving 4 cell-edge UE (-10 dB pre-tx-BF SNR) and target for
% 8 dB post-BF SINR each (supporting 1e-3 uncoded BER in QPSK)

clear;clc;warning off
rng(1)
useQN = 1;
useEQN = 0;

% ------------ System parameters ----------------------

% opt power scaling with various parameters (from DA_required_pow_theo.m)
% case1 (U=8,16,32): -11.4, -19.6343 -23.2542
% case2 (U=2,4,8): -0.1957 -2.2548 -3.20
% case3 (U=1): -4.17
power_scale = 10^(-4.17/10);

M = 1; % number of stream (mmW UE)
case_index = 3; % which use case in the paper
isNLOS = 0;
[SNR_origin,SINR_target] = get_target_SINR(case_index,M);

% Normalization: AWGN has unit power
awgn_pow = 10^(-SNR_origin/10);

% Receiver planar antenna dimension
Nr_azim = 4;
Nr_elev = 2;

% The optimal value can be evaluated in DA_required_pow_theo.m
% case1 (U=8,16,32): -10, -10, 10
% case2 (U=8,16,32): 30 30 30
%,case3 (U=1) 30 
BF_alpha_dB = 30;
BF_alpha = 10^(BF_alpha_dB/10);

% Transmitter antenna size range (fixed)
Nt_azim_range = [32];
Nt_elev_range = [8];

Nt_range = Nt_elev_range .* Nt_azim_range;
DAC_bits_range = 3:10;

% bits_range = 5;

MCtimes = 30;
NMSE = zeros(length(Nt_azim_range), length(DAC_bits_range), MCtimes);

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

NMSE_EQN = zeros(length(Nt_range), length(DAC_bits_range), MCtimes);
SINR_EQN = zeros(length(Nt_range), length(DAC_bits_range), MCtimes);
NMSE_EQN = zeros(length(Nt_range), length(DAC_bits_range), MCtimes);
SINR_EQN = zeros(length(Nt_range), length(DAC_bits_range), MCtimes);
EQN_pow_rx = zeros(length(Nt_range), length(DAC_bits_range), MCtimes);

%------------- for loop of Monte Carlo simulation over MCtimes run -------------------
for MCindex = 1:MCtimes
    
    clc;fprintf('Monte Carlo Iteration %d/%d\n',MCindex,MCtimes)

    % Baseband Waveform using 16-QAM and oversampling ratio of 5
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
        sig_pow_scale(:,mm) = sig * sig_pow_V;
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
        
        Nt_azim = Nt_azim_range(tt);
        Nt_elev = Nt_elev_range(tt);
        Nt = Nt_azim * Nt_elev;
        chan_raw = zeros(M, Nt);
        chan = zeros(M, Nt);
        
        % -------------  for loop for each UE ---------------
        for mm=1:M
            % Find range of parameter to generate channel
            temp_range = zeros(cluster_num,1);
            temp_range = (mm-1)*cluster_num+1:mm*cluster_num;

            % -------------    Get mmW 3D Channel    ---------------
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
            chan(mm,:) = chan_raw(mm,:)./norm(chan_raw(mm,:))*sqrt(Nt);
        end
        
        % Regularized Zero-Forcing in MU-MIMO channel
        precoding_mtx = zeros(Nt, M);
        precoding_mtx = chan'*inv(chan*chan'+eye(M)*BF_alpha*sqrt(Nt));
        
        % Generate symbol sequences using precoding matrix
        array_sig_raw = zeros(Nt,L);
        array_sig = zeros(Nt, L);
        
        % Floating point operation of BB precoding 
%         array_sig_raw =  precoding_mtx * sig_pow_scale.';
        
        for bb=1:length(DAC_bits_range)
            bits_num = DAC_bits_range(bb);
            
            % Fixed point operation of BB precoding
            max_precoder = 2*max(max(max(real(precoding_mtx))),max(max(imag(precoding_mtx))));
            precoding_mtx_fpt = DAC_quan(precoding_mtx, bits_num+3, max_precoder);

            max_symbol = 2*max(max(max(real(sig_pow_scale))),max(max(imag(sig_pow_scale))));
            sig_pow_scale_fpt = DAC_quan(sig_pow_scale, bits_num+3, max_symbol);

            array_sig_raw =  precoding_mtx_fpt * sig_pow_scale_fpt.';

            array_sig = array_sig_raw./norm(array_sig_raw,'fro')*sqrt(L*Nt);
        
%         pow_store(MCindex) = 20*log10((1/norm(array_sig_raw,'fro')));

            
            
            % ------- Simulate actual DAC Quantization -----------

            array_sig_QN_raw = zeros(Nt,sig_length);
            array_sig_QN = zeros(Nt,sig_length);
            max_amp(MCindex) = max(max(max(abs(real(array_sig(:,:))))),max(max(abs(imag(array_sig(:,:))))));
            maxamp = 5;%sqrt(10);
            for ll=1:sig_length
                for aa = 1:Nt
                    array_sig_QN_raw(aa,ll) = DAC_quan(array_sig(aa,ll), bits_num, maxamp);
                end
            end
            array_sig_QN = array_sig_QN_raw./norm(array_sig_QN_raw,'fro')*sqrt(L)*sqrt(power_scale);

%             
            % Signal Distortion at Specific Angle (stream 1 & Angle 1)
            for mm=1:M
    %             mm = 1;
                rx_sig_QN = zeros(sig_length, 1);
                rx_sig_EQN = zeros(sig_length, 1);
                noise_pow = norm(chan(mm,:))^2/Nt * awgn_pow;
                AWGN = (randn(L,1)+1j*randn(L,1))/sqrt(2)*sqrt(noise_pow);


                rx_sig_QN = (chan(mm,:) * array_sig_QN).' + AWGN;

                %  Normalize the Beamforming Gain (real gain)
                %  ------- min E||sig_org - sig_rx * alpha||^2 -----------
                sig_pow_norm = sig_pow_scale(:,mm)./norm(sig_pow_scale(:,mm))*sqrt(L);

                alpha_hat_QN = pinv(sig_pow_norm) * rx_sig_QN;
                gain_sig_QN(tt,bb,MCindex) = abs(alpha_hat_QN)^2;
                NMSE_QN(tt,bb,MCindex) = norm(alpha_hat_QN * sig_pow_norm...
                    - rx_sig_QN)^2/norm(alpha_hat_QN * sig_pow_norm)^2;
                noise_pow = norm(chan(mm,:))^2/Nt * awgn_pow;
                IpN_QN(mm,bb,MCindex) = norm(alpha_hat_QN * sig_pow_norm - rx_sig_QN)^2/L;
                SINR_QN(mm,bb,MCindex) = abs(alpha_hat_QN)^2/(IpN_QN(mm,bb,MCindex));
            end
        end
    end
end

%% Take mean among MC realizations
for mm=1:M
    for bb=1:length(DAC_bits_range)
        if useQN
%             NMSE_mean_QN(tt,bb) = mean(squeeze(NMSE_QN(tt,bb,:)));
            SINR_mean_QN(mm,bb) = mean(squeeze(SINR_QN(mm,bb,:)));
        end

    end
end

%% SINR
figure
plot(DAC_bits_range, 10*log10(mean(SINR_mean_QN,1)'),'-o','linewidth',2);
grid on
xlabel('DAC Quantization Number')
ylabel('SINR (dB)')

%% SE
figure
% plot(DAC_bits_range, sum(log2(1+(SINR_mean_QN)),1)','-o','linewidth',2);
plot(DAC_bits_range, M*log2(1+mean(SINR_mean_QN,1))','-o','linewidth',2);

grid on
xlabel('DAC Quantization Number')
ylabel('Spectral Efficiency (bps/Hz)')

%% debug
% temp = (array_sig+EQN)./norm(array_sig_EQN_raw,'fro');
% 
% figure
% plot(abs(chan(1,:)*temp))
% hold on
% plot(abs(chan(1,:)*array_sig_EQN))
%%
% figure
% subplot(211)
% plot(10*log10(squeeze(IpN_QN(1,1,:))));hold on
% ylabel('Gain (dB)')
% title('Inteference Power Gain')
% 
% 
% subplot(212)
% plot(10*log10(squeeze(gain_sig_QN(1,1,:))));hold on
% plot(10*log10(ones(MCtimes,1)*mean(squeeze(gain_sig_QN(1,1,:)))));hold on
% ylabel('Gain (dB)')
% title('Signal Power Gain')