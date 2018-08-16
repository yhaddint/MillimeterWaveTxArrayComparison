function [ H_MIMO ] = get_H_MIMO_3d(        raygain,...
                                            rayAOD_azim,...
                                            rayAOD_elev,...
                                            rayAOA_azim,...
                                            rayAOA_elev,...
                                            cluster_num,...
                                            ray_num,...
                                            Nt_azim,...
                                            Nt_elev,...
                                            Nr_azim,...
                                            Nr_elev)
%GET_H_FREQ Summary of this function goes here
%   Detailed explanation goes here
    
    H_MIMO = zeros(Nr_elev * Nr_azim, Nt_elev * Nt_azim);
    for cluster_index = 1:cluster_num
        for ray_index = 1:ray_num
            
            theta_azim = rayAOD_azim(cluster_index, ray_index);
            atx_azim = exp(1j*(0:Nt_azim-1)'*pi*sin(theta_azim))/sqrt(Nt_azim);
            
            theta_elev = rayAOD_elev(cluster_index, ray_index);
            atx_elev = exp(1j*(0:Nt_elev-1)'*pi*sin(theta_elev))/sqrt(Nt_elev);
            
            phi_azim = rayAOA_azim(cluster_index, ray_index);
            arx_azim = exp(1j*(0:Nr_azim-1)'*pi*sin(phi_azim))/sqrt(Nr_azim);
            
            phi_elev = rayAOA_elev(cluster_index, ray_index);
            arx_elev = exp(1j*(0:Nr_elev-1)'*pi*sin(phi_elev))/sqrt(Nr_elev);
            
            H_MIMO = H_MIMO + raygain(cluster_index, ray_index) * ...
                        kron(arx_elev, arx_azim) * kron(atx_elev, atx_azim)';
            
%             H_MIMO = reshape(H_MIMO, Nr_azim*Nr_elev, Nt_azim*Nt_elev);

        end
    end

end