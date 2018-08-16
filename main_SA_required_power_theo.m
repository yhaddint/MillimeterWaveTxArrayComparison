% This script provides theoretical value of required power
% The output Pout_adjust gives values in dB.
% It simulates expected signal gain, G, and multiuser inteference, I
% To reach SINR target, we need
% (P * G)/(N + P * I) > TH
% where P is the power adjustment factor (46dBm as reference)
% and N is the AWGN power, i.e., N = 1/SNR_origin.
% Simple steps show as long as G - I * TH > 0, adjusting P gives feasible
% solution for SINR to reach TH
% We use P^{\star} = (N * TH) / (G - I * TH) as required power
% Hardware impirment is not evaluated in this script.

% Different from DA and FH, the array division is important and we only
% show one specific design in the paper but simulator reveals it affects
% the required transmission power 

clear;clc;warning off
rng(3);
MCtimes = 30;

% ---- System Parameters (SNR and arrays) ------------
M = 8; % number of stream (num of mmW UE)
case_index = 2; % which use case in the paper
isNLOS = 1;
[SNR_origin,SINR_target] = get_target_SINR(case_index,M);

% Receiver planar antenna dimension
Nr_azim = 4;
Nr_elev = 2;

% Array dimension to find alphs in impairmnet figure of paper
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


% Arrat geometry used in Case I and 8 multiplexing level
% Nt_range =         [64, 128, 256, 512, 768, 1024];
% Nt_azim_range =    [32,  32,  32,  32,  32,   32];
% Nt_elev_range =    [ 2,   4,   8,  16,  24,   32];
% Nt_azim_SA_range = [ 8,  16,  32,  32,  32,   32];
% Nt_elev_SA_range = [ 1,   1,   1,   2,   3,    4];

% % Arrat geometry used in Case I and 16 multiplexing level
% Nt_range =         [64, 128, 256, 512, 768, 1024];
% Nt_azim_range =    [32,  32,  32,  32,  48,   32];
% Nt_elev_range =    [ 2,   4,   8,  16,  16,   32];
% Nt_azim_SA_range = [ 4,   8,  16,  16,  24,   32];
% Nt_elev_SA_range = [ 1,   1,   1,   2,   2,    2];


% Arrat geometry used in Case I and 32 multiplexing level
% Nt_range =         [64, 128, 256, 512, 768, 1024];
% Nt_azim_range =    [32,  32,  32,  32,  48,   32];
% Nt_elev_range =    [ 2,   4,   8,  16,  16,   32];
% Nt_azim_SA_range = [ 2,   4,   8,  16,  24,   32];
% Nt_elev_SA_range = [ 1,   1,   1,   1,   1,    1];


% % Arrat geometry used in Case II and 2 multiplexint level
% Nt_range =         [64, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048];
% Nt_azim_range =    [32,  32,  32,  32,  32,  32,  32,   32,   48,   64];
% Nt_elev_range =    [ 2,   4,   6,   8,  12,  16,  24,   32,   32,   32];
% Nt_azim_SA_range = [32,  32,  32,  32,  32,  32,  32,   32,   48,   64];
% Nt_elev_SA_range = [ 1,   2,   3,   4,   6,   8,  12,   16,   16,   16];


% Arrat geometry used in Case II and 4 multiplexint level
% Nt_range =         [64, 128, 256, 384, 512, 768, 1024, 1536, 2048];
% Nt_azim_range =    [32,  32,  32,  32, 32,  32,   32,    48,   64];
% Nt_elev_range =    [ 2,   4,   8,  12, 16,  24,   32,    32,   32];
% Nt_azim_SA_range = [16,  32,  32,  32, 32,  32,   32,    48,   64];
% Nt_elev_SA_range = [ 1,   1,   2,   3,  4,   6,    8,     8,    8];


% % Arrat geometry used in Case II and 8 multiplexint level
% Nt_range =         [64, 128, 256, 512, 768, 1024, 1536, 2048];
% Nt_azim_range =    [32,  32,  32,  32,  32,   32,   48,   64];
% Nt_elev_range =    [ 2,   4,   8,  16,  24,   32,   32,   32];
% Nt_azim_SA_range = [ 8,  16,  32,  32,  32,   32,   48,   64];
% Nt_elev_SA_range = [ 1,   1,   1,   2,   3,    4,    4,    4];

% Regularized-ZF coefficient alpha
alpha_range = linspace(-30,30,31);

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
gain_sig = zeros(length(Nt_range),length(alpha_range),MCtimes);
gain_MUint = zeros(length(Nt_range),length(alpha_range),MCtimes);

% ------------- Monte Carlo Simulations -------------
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
            
            % The 3D MIMO channel has been re-shapred to the following
            % kron(arx_elev, arx_azim) * kron(atx_elev, atx_azim)';
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

        % -------------    Hybrid BF    ---------------
        W_RF_PS_raw = zeros(1,Nt);
        W_RF_PS = zeros(1,Nt);
        W_RF = zeros(Nt, M);

        for mm=1:M
            sa_index = (mm-1)*Nt_SA+1:mm*Nt_SA;
            W_RF_PS_raw(mm,:) = angle(chan(mm,:));
            W_RF_PS = W_RF_PS_raw(mm,:);
            W_RF(sa_index,mm) = exp(-1j*transpose(W_RF_PS(sa_index)));
        end

        for m1=1:M
            for m2=1:M
                chan_effct(m1, m2) = chan(m1,:) * W_RF(:,m2);
            end
        end


        for aa=1:length(alpha_range)

            alpha = 10^(alpha_range(aa)/10);

            W_BB = zeros(M,M);
            W_BB = chan_effct'*inv(chan_effct*chan_effct'+eye(M)*alpha*sqrt(Nt));

            W_hybrid_raw = W_RF *W_BB;
            W_hybrid = W_hybrid_raw./norm(W_hybrid_raw,'fro');
            array_sig_noQN =  W_hybrid * sig_pow_scale.';
            % array_sig is Nt by L matrix (tx signal at all antenna and L samples)
            
            array_sig = array_sig_noQN./norm(array_sig_noQN,'fro')*sqrt(L);

            % -------------  for loop for all UEs  -------------
            for mm = 1:M
                rx_sig = (chan(mm,:) * array_sig).';
                % rx_sig is L by 1 vector (m-th rx signal at RF-chain with L samples)
                
                %  Normalize the Beamforming Gain (real gain)
                %  ------- min E||alpha * sig_org - sig_rx||^2 -----------
                sig_pow_norm = sig_pow_scale(:,mm)./norm(sig_pow_scale(:,mm))*sqrt(L);
                gain_chan(mm, MCindex) = norm(chan(mm,:))^2;
                g_hat = pinv(sig_pow_norm) * rx_sig;
                gain_sig_temp(mm) = abs(g_hat)^2;
                gain_int_temp(mm) = norm(g_hat * sig_pow_norm - rx_sig)^2/L;
            end
            
            % Average over all UE for signal gain and interference gain
            gain_sig(tt, aa, MCindex) = mean(gain_sig_temp);
            gain_MUint(tt, aa, MCindex) = mean(gain_int_temp);
        end
    end
end

%% Signal Gain and MU-Interference Evaluation

Int_mean = zeros(M,1);
gain_sig_mean = zeros(M,1);

for tt=1:length(Nt_azim_range)
    for aa=1:length(alpha_range)
        gain_sig_mean(tt,aa) = mean(squeeze(gain_sig(tt,aa,:)));
        Int_mean(tt,aa) = mean(squeeze(gain_MUint(tt,aa,:)));
    end
end

for tt=1:length(Nt_azim_range)
% for tt=1
    sig_gain = 10*log10(gain_sig_mean(tt,:));
    int_gain = 10*log10(Int_mean(tt,:));
    SIR = sig_gain - int_gain;
    index = ((gain_sig_mean(tt,:)-10^(SINR_target/10)*Int_mean(tt,:))>0);
    temp = (10^((SINR_target-SNR_origin)/10))./(gain_sig_mean(tt,index)-10^(SINR_target/10)*Int_mean(tt,index));
%     figure;plot(alpha_range(index),10*log10(temp))
    Pout_adjust(tt) = min(10*log10(temp)); 
end
%%
% SINR_now = (min(temp)*gain_sig_mean(tt,index))./(10^(-SNR_origin/10) + min(temp)*Int_mean(tt,index));
% figure
% plot(alpha_range(index), 10*log10(SINR_now));grid on
% xlable('alpha')
% ylabel('SINR')

figure
plot(alpha_range(index),10*log10(temp))
