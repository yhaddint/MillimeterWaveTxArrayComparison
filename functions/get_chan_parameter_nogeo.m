function [ raygain,...
           raydelay,...
           ray_AOA_azim,...
           ray_AOD_azim,...
           ray_AOA_elev,...
           ray_AOD_elev ] =...
    get_chan_parameter_nogeo( print_stat,...
                              cluster_num,...
                              ray_num,...
                              sigma_delay_spread,...
                              sigma_AOA_az_spread,...
                              centroid_AOD_az,...
                              sigma_AOD_az_spread,...
                              sigma_AOA_el_spread,...
                              sigma_AOD_el_spread,...
                              isNLOS,...
                              RacianK)
%GET_CHAN_PARAMETER_NOGEO Summary of this function goes here
%   Generate channel parameter using given statistics
%   [ raygain, raydelay, ray_AOA_azim, ray_AOD_azim ] = ...
%           get_chan_parameter_nogeo( print_stat, cluster_num, ray_num, sigma_delay_spread,...
%                                      centroid_AOA, sigma_AOA_spread,...
%                                      centroid_AOD, sigma_AOD_spread )
%   IP: print_stat, indicator for printing channel statistics
%   IP: cluster_num, number of MPC cluster
%   IP: ray_num, number of rays in each clusters
%   IP: sigma_delay_spread, desired intra-cluster delay spread (std with
%       unit second)
%   IP: sigma_AOA_az_spread, desired intra-cluster AOA spread (std var with unit rad)
%   IP: centroid_AOD_az, desired centroid of AOD of clusters (with unit rad), can
%       be 'random' if not to be fixed
%   IP: sigma_AOD_az_spread, desired intra-cluster az AOD spread (std var with unit rad)
%   IP: sigma_AOA_el_spread, desired intra-cluster el AOA spread (std var with unit rad)
%   IP: sigma_AOD_el_spread, desired intra-cluster el AOD spread (std var with unit rad)
%   OP: raygain, cluster_num by ray_num matrix with raygain for each ray
%   OP: raydelay, cluster_num by ray_num matrix with delay for each ray
%   OP: ray_AOA_azim, cluster_num by ray_num matrix with azimuth AOA for each ray
%   OP: ray_AOD_azim, cluster_num by ray_num matrix with azimuth AOD for each ray
%   OP: ray_AOA_elev, cluster_num by ray_num matrix with elevation AOA for each ray
%   OP: ray_AOD_elev, cluster_num by ray_num matrix with elevation AOD for each ray


if cluster_num>3
    fprintf('This function only support up to 3 multipath clusters\n');
end

if max(sigma_delay_spread)>1e-6
    fprintf('Delay spread has beter not to exceed 1e-6 second\n');
end

cluster_delay = [300e-9,250e-9,200e-9]; % Mean delay of two multipath clusters

% Zero initializations
raygain = zeros(cluster_num, ray_num);
raydelay = zeros(cluster_num, ray_num);
ray_AOA_azim = zeros(cluster_num, ray_num);
ray_AOA_elev = zeros(cluster_num, ray_num);
ray_AOD_azim = zeros(cluster_num, ray_num);
ray_AOD_elev = zeros(cluster_num, ray_num);


for cluster_index = 1:cluster_num

    % -------- Randomly generate chan. parameters (AOD in Azimuth) --------
    if cluster_index == 1
        ray_AOD_azim_mean = centroid_AOD_az;
        % Only first cluster (LOS) has well-seperated angle due to
        % scheduler assumption
    else
        ray_AOD_azim_mean = rand * 2 * pi / 3 - pi/3;
    end
    angle_spread = laprnd(1, ray_num, 0, sigma_AOD_az_spread);
%     angle_spread = randn(1, ray_num) * sigma_AOD_spread;
    ray_AOD_azim(cluster_index,:) = ray_AOD_azim_mean + angle_spread;

    
    % ------ Randomly generate chan. parameters (AOD in Elevation) --------
    ray_AOD_elev_mean = rand * pi / 3 - pi/6;
    angle_spread = laprnd(1, ray_num, 0, sigma_AOD_el_spread);
    ray_AOD_elev(cluster_index,:) = ray_AOD_elev_mean + angle_spread;

    
    % -------- Randomly generate chan. parameters (AOA in Azimuth) --------
    ray_AOA_azim_mean = rand * 2 * pi / 3 - pi/3;
    angle_spread = laprnd(1, ray_num, 0, sigma_AOA_az_spread);
%     angle_spread = randn(1, ray_num) * sigma_AOA_spread;
    ray_AOA_azim(cluster_index,:) = ray_AOA_azim_mean + angle_spread;

    
    % -------- Randomly generate chan. parameters (AOA in Elevation) --------
    ray_AOA_elev_mean = rand * pi / 3 - pi/6;
    angle_spread = laprnd(1, ray_num, 0, sigma_AOA_el_spread);
    ray_AOA_elev(cluster_index,:) = ray_AOA_elev_mean + angle_spread;


    % -------- Unit gain for each ray --------
    if cluster_index==1
        raygain(cluster_index,:) = exp(1j*rand(1,ray_num)*2*pi)/sqrt(ray_num);  
    else
        raygain(cluster_index,:) = exp(1j*rand(1,ray_num)*2*pi)/sqrt(ray_num)/sqrt(10^(RacianK/10)); 
    end
    
    % -------- Delay of each ray, Gaussian distributed --------
    raydelay(cluster_index,:) = cluster_delay(cluster_index) + randn(1,ray_num)*sigma_delay_spread;
    
    if print_stat
        fprintf('Cluster %d:\n',cluster_index);
        fprintf('DS Mean:    %4.2f ns\n', mean(raydelay(cluster_index,:)*1e9));
        fprintf('DS Std Dev: %4.2f ns\n', sqrt(var(raydelay(cluster_index,:)*1e9)));
        fprintf('AZ AOAS Mean:    %4.2f deg\n',...
            mean(ray_AOA_azim(cluster_index,:)/pi*180));
        fprintf('AZ AOAS Std Dev: %4.2f deg\n',...
            sqrt(var(ray_AOA_azim(cluster_index,:)/pi*180)));
        fprintf('EL AOAS Mean:    %4.2f deg\n',...
            mean(ray_AOA_elev(cluster_index,:)/pi*180));
        fprintf('EL AOAS Std Dev: %4.2f deg\n',...
            sqrt(var(ray_AOA_elev(cluster_index,:)/pi*180)));
        fprintf('AZ AODS Mean:    %4.2f deg\n',...
            mean(ray_AOD_azim(cluster_index,:)/pi*180));
        fprintf('AZ AODS Std Dev: %4.2f deg\n',...
            sqrt(var(ray_AOD_azim(cluster_index,:)/pi*180)));
        fprintf('EL AODS Mean:    %4.2f deg\n',...
            mean(ray_AOD_elev(cluster_index,:)/pi*180));
        fprintf('EL AODS Std Dev: %4.2f deg\n',...
            sqrt(var(ray_AOD_elev(cluster_index,:)/pi*180)));
    end
end

% In NLOS scenario, delete the first cluster (LOS cluster)
if isNLOS
    raygain(1,:) = zeros(1,ray_num);
end 


end

