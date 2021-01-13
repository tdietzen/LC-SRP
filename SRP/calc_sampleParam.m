function [ T, N_mm ] = calc_sampleParam(micPos, w_0, c)
% [ T, N_mm ] = calc_sampleParam(micPos, w_0, c)
% computes sampling period and number of samples within TDOA interval.
%
% IN:
% micPos         microphone positions - channels x coordinates
% w_0            band limit
% c              speed of sound
%
% OUT:
% T              sampling period
% N_mm           TDOAs - microphone pairs


M = size(micPos,1);
% microphone pairs
P = M*(M-1)/2;
% sampling period
T = pi/w_0;

dist = zeros(P,1);
p = 0;
for mprime = 1:M
    for m = mprime+1:M
        p = p+1;
        % distance between microphones
        dist(p) = norm(micPos(m,:) - micPos(mprime,:));        
    end
end

% distance between samples
Delta_t_0 = dist/c;

% samples inside TDOA interval
N_mm = floor(Delta_t_0/T);

end

