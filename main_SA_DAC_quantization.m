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

% ------------ System parameters ----------------------

% opt power scaling with various parameters (from DA_required_pow_theo.m)
% case1 (U=8,16,32): 0.1984 -2.0350 -9.8275
% case2 (U=8,16,32): 3.5939 4.56 6.0829
% case3 (U=18): xx
power_scale = 10^(6.05/10);

% Receiver planar antenna dimension
Nr_azim = 4;
Nr_elev = 2;

M = 8; % number of stream (mmW UE)
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
BF_sigma_dB = 22;
BF_sigma = 10^(BF_sigma_dB/10);

% DSP/DAC finite number of bits
bits_range = 3:1:10;

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
        PS_bits_num = 6;
        PS_scale = 2^(PS_bits_num-1);
        W_RF_PS_raw = zeros(1, Nt);
        W_RF_PS = zeros(1, Nt);
        W_RF = zeros(Nt, M);
        
        for mm=1:M
            sa_index = (mm-1)*Nt_SA + 1:mm*Nt_SA;
            W_RF_PS_raw(mm,:) = angle(chan(mm,:));
            W_RF_PS = (round(W_RF_PS_raw(mm,:)/pi*PS_scale+0.5)-0.5)/PS_scale*pi;
            W_RF(sa_index,mm) = exp(-1j*transpose(W_RF_PS(sa_index)));
        end
        
        for m1=1:M
            for m2=1:M
                chan_effct(m1, m2) = chan(m1,:) * W_RF(:,m2);
            end
        end
        
        W_BB = zeros(M,M);
        W_BB = chan_effct'*inv(chan_effct*chan_effct'+eye(M)*BF_sigma*sqrt(Nt));
%         W_BB = inv(chan_effct);
        
        for bb=1:length(bits_range)
            bits_num = bits_range(bb);
            
%             % ------- Adding Effective Quantization Noise ---------
%             if useEQN
%                 array_sig_EQN = zeros(Nt,sig_length);
%                 array_sig_raw = zeros(Nt, L);
%                 array_sig = zeros(Nt, L);
% 
%                 array_sig_raw =  W_BB_norm * sig_pow_scale.';
%                 for aa = 1:Nt
%                     % Adjust offset because original signal amp is not [-1, 1]
%                     bits_offset = 2.3;
%                     QN_pow = 10 ^ (- (bits_num + bits_offset) * 6.02 / 10);
%                     EQN_raw = (rand(sig_length,1) - 0.5) ...
%                             + 1j * (rand(sig_length,1) - 0.5);
%                     EQN_normalized = EQN_raw./norm(EQN_raw) * sqrt(sig_length);
%                     EQN =  EQN_normalized * sqrt(QN_pow) ;
%                     array_sig_EQN(aa,:) = array_sig(aa,:) + EQN.';
%                 end
%             end
            
            % ------- Simulate actual DAC Quantization -----------
            if useQN
                array_sig_QN_raw = zeros(Nt,sig_length);
                array_sig_QN = zeros(Nt,sig_length);
                array_sig_raw = zeros(Nt, L);
                array_sig = zeros(Nt, L);
                
                % ------ Fixed point operation of BB precoding ---------
                max_precoder = 2*max(max(max(real(W_BB))),max(max(imag(W_BB))));
                W_BB_fpt = DAC_quan(W_BB, bits_num+3, max_precoder);
                max_symbol = 2*max(max(max(real(sig_pow_scale))),max(max(imag(sig_pow_scale))));
                sig_pow_scale_fpt = DAC_quan(sig_pow_scale, bits_num+3, max_symbol);
                array_sig_raw =  W_BB_fpt * sig_pow_scale_fpt.';
                array_sig = array_sig_raw./norm(array_sig_raw)*sqrt(L*M);

                % ------ Quantization on Digitally Precoded Signals ---------
                maxamp = 5;
                for ll=1:sig_length
                    for mm = 1:M
                        array_BB_QN(mm,ll) = DAC_quan(array_sig(mm,ll), bits_num, maxamp);
                    end
                end
                array_sig_QN_raw = W_RF * array_BB_QN;
%                 array_sig_QN_raw = W_RF * array_sig;
                array_sig_QN = array_sig_QN_raw./norm(array_sig_QN_raw,'fro')...
                    *sqrt(L)*sqrt(power_scale);
                % sqrt(4) because total power is 4 times of digital array!

%             % ------- Power Statistics of Each PA (Bits = 7) -----------
%                 if bits_num == 7
%                     ave_pow = zeros(Nt,1);
%                     for aa=1:Nt
%                         ave_pow(aa) = mean(abs(array_sig_QN(aa,:).^2));
%                     end
%                     PAPR_array(tt,MCindex) = max(ave_pow) / mean(ave_pow);
%                 end
% 
%             % --------- Evaluate Quantization Noise Power -----------
%                 if Nt == 16
%                     QEN_pow_evl(bb,MCindex) = mean(abs(array_sig_QN(1,:) - array_sig(1,:)).^2);
%                 end
% 
%             % --------- Evaluate Quantization Noise xcorr (Bits = 7)-----------
%                 if Nt == 16 && bits_num == 7
%                     for aa = 1:Nt
%                         QN(:,aa) = array_sig_QN(aa,:).' - array_sig(aa,:).';
%                     end
%                     QNcorr_mtx(:,:,MCindex) = QN'*QN / sig_length;
%                 end
            end
            
            % not simulating quantization 
            if noQN
                W_hybrid_raw = W_RF *W_BB;
                W_hybrid = W_hybrid_raw./norm(W_hybrid_raw,'fro')*sqrt(power_scale);
                array_sig_noQN =  W_hybrid * sig_pow_scale.';
%                 array_sig_noQN = array_sig_noQN_raw./norm(array_sig_noQN_raw,'fro')*sqrt(L)*sqrt(4);
            end

%             %% quick watch quantization noise (it is white!)
%             temp = 20*log10(abs(fft(hann(sig_length).*QN(:,1))));
%             figure
%             plot([temp(sig_length/2+1:end);temp(1:sig_length/2)])
%             
            %% Signal Distortion at Specific Angle
            for mm = 1:M
            
            rx_sig_QN = zeros(sig_length, 1);
            rx_sig_EQN = zeros(sig_length, 1);
            noise_pow = norm(chan(mm,:))^2/Nt * awgn_pow;
            AWGN = (randn(L,1)+1j*randn(L,1))/sqrt(2)*sqrt(noise_pow);
            
            if noQN
                rx_sig_noQN = (chan(mm,:) * array_sig_QN).' + AWGN;
                % rx_sig is L by 1 vector (m-th rx signal at RF-chain with L samples)
                
%                 for m1=1:M
%                     sa_index = (m1-1)*Nt+1:m1*Nt;
%                     rx_sig_noQN = rx_sig_QN + (chan(mm,:) * array_sig_noQN(sa_index,:).') + AWGN;
%                 end
                
                % Evaluate signal and non-signal part
                %  ------- min E||alpha * sig_org - sig_rx||^2 -----------
                sig_pow_norm = sig_pow_scale(:,mm)./norm(sig_pow_scale(:,mm))*sqrt(L);
                gain_chan(mm,bb,MCindex) = norm(chan(mm,:))^2;
                alpha_hat_noQN = pinv(sig_pow_norm) * rx_sig_noQN;
                gain_sig_noQN(mm,bb,MCindex) = abs(alpha_hat_noQN)^2;
                NMSE_noQN(tt,bb,MCindex) = norm(alpha_hat_noQN * sig_pow_norm...
                    - rx_sig_noQN)^2/norm(alpha_hat_noQN * sig_pow_norm)^2;
                noise_pow = norm(chan(mm,:))^2/Nt * awgn_pow;
                IpN_noQN(mm,bb,MCindex) = norm(alpha_hat_noQN * sig_pow_norm - rx_sig_noQN)^2/L;
                SINR_noQN(mm,bb,MCindex) = abs(alpha_hat_noQN)^2/(IpN_noQN(mm,bb,MCindex));
            end
            
            
            if useQN
                rx_sig_QN = (chan(mm,:) * array_sig_QN).' + AWGN;
                % rx_sig is L by 1 vector (m-th rx signal at RF-chain with L samples)
                
                % Evaluate signal and non-signal part
                %  ------- min E||alpha * sig_org - sig_rx||^2 -----------
                sig_pow_norm = sig_pow_scale(:,mm)./norm(sig_pow_scale(:,mm))*sqrt(L);
                gain_chan(mm,bb,MCindex) = norm(chan(mm,:))^2;
                alpha_hat_QN = pinv(sig_pow_norm) * rx_sig_QN;
                gain_sig_QN(mm,bb,MCindex) = abs(alpha_hat_QN)^2;
                NMSE_QN(tt,bb,MCindex) = norm(alpha_hat_QN * sig_pow_norm...
                    - rx_sig_QN)^2/norm(alpha_hat_QN * sig_pow_norm)^2;
                noise_pow = norm(chan(mm,:))^2/Nt * awgn_pow;
                IpN_QN(mm,bb,MCindex) = norm(alpha_hat_QN * sig_pow_norm - rx_sig_QN)^2/L;
                SINR_QN(mm,bb,MCindex) = abs(alpha_hat_QN)^2/(IpN_QN(mm,bb,MCindex));
            end
            
%             if useEQN
%                 for ll = 1:sig_length
%                     rx_sig_EQN(ll) = chan(mm,:) * array_sig_EQN(:,ll);
%                 end
%                 alpha_hat_EQN = pinv(rx_sig_EQN) * sig_pow_scale(:,mm);
%                 NMSE_EQN(tt,bb,MCindex) = norm(alpha_hat_EQN * rx_sig_EQN...
%                     - sig_pow_scale(:,mm))^2/norm(sig_pow_scale)^2;
%             end
            end
        end
    end
end

%% Take mean among MC realizations

SINR_mean_QN = zeros(M,length(bits_range));
for bb=1:length(bits_range)
    for mm=1:M
        if useQN
%                 NMSE_mean_QN(mm,bb) = mean(squeeze(NMSE_QN(mm,bb,:)));
            SINR_mean_QN(mm,bb) = mean(squeeze(SINR_QN(mm,bb,:)));
        end
%             if useEQN
%                 NMSE_mean_EQN(bb) = mean(squeeze(NMSE_EQN(tt,bb,:)));
%             end
    end
end


%% Take mean among MC realizations
% if useQN
%     QEN_pow_dB_mean = 10*log10(mean(QEN_pow_evl,2));
% end

%% QN corr matrix
% if useQN
%     QNcorr_mtx_mean = zeros(4,4);
%     for aa = 1:4
%         for bb=1:4
%             QNcorr_mtx_mean(aa,bb) = abs(mean(squeeze(QNcorr_mtx(aa,bb,:))));
%         end
%     end
% end
%% Array sense PAPR 
% if useQN
%     figure(99)
%     for tt = 1:length(Nt_range)
%         [a, b] = ecdf(10*log10(PAPR_array(tt,:)));
%         plot(b,1-a,'linewidth',2);hold on
%     end
%     xlabel('Array Aspect PAPR (dB)')
%     ylabel('Probability')
%     legend('Nt = 4  (4 \times 1)',...
%        'Nt = 8  (4 \times 2)',...
%        'Nt = 16 (8 \times 2)',...
%        'Nt = 32 (8 \times 4)',...
%        'Nt = 64 (16 \times 4)' );
% end

%% sig and EQN gain
% if useQN
% for tt=1:length(Nt_range)
%     for bb=1:length(bits_range)
% %         gain_sig_EQN_mean(tt,bb) = mean(squeeze(gain_sig_EQN(tt,bb,:)));
%         gain_chan_mean(tt,bb) = mean(squeeze(gain_chan(tt,bb,:)));
%     end
% end
% 
% % IpN_EQN_mean = zeros(M,length(bits_range));
% IpN_QN_mean = zeros(M,length(bits_range));
% 
% 
% for mm=1:M
%     for bb=1:length(bits_range)
% %         IpN_EQN_mean(mm,bb) = mean(squeeze(IpN_EQN(mm,bb,:)));
%         IpN_QN_mean(mm,bb) = mean(squeeze(IpN_QN(mm,bb,:)));
%         gain_sig_QN_mean(mm,bb) = mean(squeeze(gain_sig_QN(mm,bb,:)));
%     end
% end
% 
% figure
% % plot(bits_range, 10*log10(gain_sig_EQN_mean'),'--x','linewidth',2);hold on
% plot(bits_range, 10*log10(gain_sig_QN_mean'),'-o','linewidth',2)
% grid on
% ylabel('Gain (dB)')
% xlabel('DAC Quantization Number')
% legend('EQN','QN')
% title('Signal Gain')
% 
% figure
% % plot(bits_range, 10*log10(IpN_EQN_mean'),'--x','linewidth',2);hold on
% plot(bits_range, 10*log10(IpN_QN_mean'),'-o','linewidth',2)
% grid on
% xlabel('DAC Quantization Number')
% ylabel('Gain (dB)')
% title('MUint+Tx Noise Power w/ EQN and QN')
% legend('EQN','QN')
% 
% % 
% % figure
% % plot(bits_range, 10*log10(gain_chan_mean'),'-o','linewidth',2)
% % grid on
% % xlabel('DAC Quantization Number')
% % ylabel('Gain (dB)')
% % 
% % title('Channel Gain')
% end
%% SINR w/ DAC QN
if useQN
    % SINR
    figure
    plot(bits_range, 10*log10(SINR_mean_QN'),'-o','linewidth',2);
    grid on
    xlabel('DAC Quantization Number')
    ylabel('SINR (dB)')

    % SE
    figure
    plot(bits_range, sum(log2(1+10^(0.1/10)*SINR_mean_QN))','-o','linewidth',2);
%     plot(bits_range, M*log2(1+10^(0/10)*SINR_mean_QN)','-o','linewidth',2);

    grid on
    xlabel('DAC Quantization Number')
    ylabel('Spectral Efficiency (bps/Hz)')

end

