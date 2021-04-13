function [ SRP_stack ] = calc_SRPconv(Psi_STFT, omega, Delta_t_i)
% [ SRP_stack ] = calc_SRPconv(Psi_STFT, omega, Delta_t_i)
% performs conventional SRP. 
%
% IN:
% Psi_STFT        FD GCCs - frequency x frames x microphone pairs
% omega           frequency vector
% Delta_t_i       TDOA vector
%
% OUT:
% SRP_stack       SRP map - frames x candidate locations


[J, P] = size(Delta_t_i);
[K, L, P] = size(Psi_STFT);

SRP_stack = zeros(L,J);
for l = 1:L
    SRP_FD = zeros(K,J);
    for k = 2:K
        psi = squeeze(Psi_STFT(k,l,:));
        SRP_FD(k,:) = exp(1j*Delta_t_i*omega(k))*psi; % (8)
    end
    SRP_stack(l,:) = 2*sum(real(SRP_FD), 1); % (9)
end

