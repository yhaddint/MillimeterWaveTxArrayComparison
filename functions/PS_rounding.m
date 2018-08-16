function [ phase_out ] = PS_rounding( phase_in, PS_bits )
%PS_ROUNDING Summary of this function goes here
%   Detailed explanation goes here
    N = length(phase_in);
    phase_out = zeros(N,1);
    step = 2^(PS_bits-1);
    for nn=1:N
        in = phase_in(nn);
        temp = ((round((in/pi)*step))+0.5)/step*pi-pi/(2^PS_bits);
        phase_out(nn) = temp;
    end

end

