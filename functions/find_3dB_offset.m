function [ offset ] = find_3dB_offset( steer_angle, Nt )
%FIND_3DB_OFFSET Summary of this function goes here
%   Detailed explanation goes here
steer_angle_num = 500;
steer_angle_range = linspace(steer_angle-pi/Nt,steer_angle,steer_angle_num);
steer_vec = exp(1j * pi * (0:Nt-1)' * sin(steer_angle));

for ss = 1:steer_angle_num
    search_angle = steer_angle_range(ss);
    spatial_vec = exp(1j * pi * (0:Nt-1)' * sin(search_angle));
    RSS(ss) = abs(steer_vec' * spatial_vec/sqrt(Nt));
    RSS_dB(ss) = 20*log10(RSS(ss));
end
[~,index] = min(abs(RSS_dB - (max(RSS_dB)-3)));
offset = steer_angle_range(index) - steer_angle;

end

