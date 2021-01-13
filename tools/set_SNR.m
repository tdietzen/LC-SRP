function [v_scaled, scaling] = set_SNR(x, v, SNR)
% [v_scaled, scaling] = set_SNR(x, v, SNR)
% scales the noise signal v relative to speech signal x in order to obtain specified SNR.
%
% IN:
% x           speech signal
% v           noise signal
% SNR         desired SNR
%
% OUT:
% v_scaled    scaled noise signal
% scaling     scaling factor


% scale signals
power_x = norm(x(:,1));
power_v  = norm(v(:,1));
if isinf(SNR)
    scaling = 0;
else
    scaling = ((power_x/power_v)/db2mag(SNR));
end
v_scaled = scaling*v;
end