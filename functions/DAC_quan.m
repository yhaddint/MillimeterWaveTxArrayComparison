function [ s_out_complx ] = DAC_quan( input_complx, bits_num, maxamp )
%DAC_QUAN Summary of this function goes here
%   Detailed explanation goes here
    
    temp_scale = 2^(bits_num-1);
    
    s_in_real = real(input_complx);
    s_in_imag = imag(input_complx);
    
    s_out_real = (round(s_in_real/maxamp * temp_scale+0.5)-0.5) / temp_scale * maxamp;
    s_out_imag = (round(s_in_imag/maxamp * temp_scale+0.5)-0.5) / temp_scale * maxamp;
    
    s_out_complx = s_out_real + 1j * s_out_imag;


end

