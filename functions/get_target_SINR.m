function [raw_SNR,target_SINR] = get_target_SINR(case_index,U)
%GET_TARGET_SINR Summary of this function goes here
%   Detailed explanation goes here

switch case_index
    case 1
        SE_target = 58.8;
        raw_SNR = 18.69;
    case 2
        SE_target = 4.7;
        raw_SNR = -14.71;
    case 3
        SE_target = 10/0.85;
        U = 1;
        raw_SNR = 15.51;
    otherwise
end

target_SINR = 10*log10(2^(SE_target/U)-1);
        
end

