function [ H_Tx_3d_rowvec ] = get_H_Tx_3d(  raygain,...
                                            rayAOD_azim,...
                                            rayAOD_elev,...
                                            cluster_num,...
                                            ray_num,...
                                            Nt_azim,...
                                            Nt_elev)
%GET_H_FREQ Summary of this function goes here
%   Detailed explanation goes here
    
    H_Tx_3d = zeros(Nt_elev, Nt_azim);
    for cluster_index = 1:cluster_num
        for ray_index = 1:ray_num
            
            theta_azim = rayAOD_azim(cluster_index, ray_index);
            atx_azim = exp(1j*(0:Nt_azim-1)'*pi*sin(theta_azim))/sqrt(Nt_azim);
            
            theta_elev = rayAOD_elev(cluster_index, ray_index);
            atx_elev = exp(1j*(0:Nt_elev-1)'*pi*sin(theta_elev))/sqrt(Nt_elev);
            
            H_Tx_3d = H_Tx_3d ...
                + raygain(cluster_index, ray_index) * atx_elev * atx_azim';
            
            H_Tx_3d_rowvec = reshape(H_Tx_3d, 1, Nt_azim*Nt_elev);

        end
    end

end