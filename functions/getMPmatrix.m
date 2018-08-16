function [ A ] = getMPmatrix( sig_in , P)
%GETMPMATRIX Summary of this function goes here
%   Detailed explanation goes here

    % s(N), s(N-1), s(N-2),...
    % s(N)|s(N)|^2, s(N-1)|s(N-1)|^2, s(N-2)|s(N-2)|^2,...
    % s(N)|s(N)|^4, s(N-1)|s(N-1)|^4, s(N-2)|s(N-2)|^4...

    sig_d0 = sig_in(3:end);
    sig_d1 = sig_in(2:end-1);
    sig_d2 = sig_in(1:end-2);
    %----------- First Order ----------------
    A_10 = sig_d0;
    A_11 = sig_d1;
    A_12 = sig_d2;

    %----------- Third Order ----------------
    A_30 = sig_d0 .* sig_d0 .* conj(sig_d0);
    A_31 = sig_d1 .* sig_d1 .* conj(sig_d1);
    A_32 = sig_d2 .* sig_d2 .* conj(sig_d2);

    %----------- Fifth Order ----------------
    A_50 = sig_d0 .* sig_d0.^2 .* conj(sig_d0).^2;
    A_51 = sig_d1 .* sig_d1.^2 .* conj(sig_d1).^2;
    A_52 = sig_d2 .* sig_d2.^2 .* conj(sig_d2).^2;

    %----------- Seventh Order ----------------
    A_70 = sig_d0 .* sig_d0.^3 .* conj(sig_d0).^3;
    A_71 = sig_d1 .* sig_d1.^3 .* conj(sig_d1).^3;
    A_72 = sig_d2 .* sig_d2.^3 .* conj(sig_d2).^3;
    
    %----------- Nineth Order ----------------
    A_90 = sig_d0 .* sig_d0.^4 .* conj(sig_d0).^4;
    A_91 = sig_d1 .* sig_d1.^4 .* conj(sig_d1).^4;
    A_92 = sig_d2 .* sig_d2.^4 .* conj(sig_d2).^4;
    %----------- A matrix -------------
    if P==3
        A = [A_10, A_11, A_12,...
             A_30, A_31, A_32]; 
    elseif P==5
        A = [A_10, A_11, A_12,...
             A_30, A_31, A_32,...
             A_50, A_51, A_52];
    elseif P==7
        A = [A_10, A_11, A_12,...
             A_30, A_31, A_32,...
             A_50, A_51, A_52,...
             A_70, A_71, A_72];
    elseif P==9
        A = [A_10, A_11, A_12,...
             A_30, A_31, A_32,...
             A_50, A_51, A_52,...
             A_70, A_71, A_72,...
             A_90, A_91, A_92];
    end
end

