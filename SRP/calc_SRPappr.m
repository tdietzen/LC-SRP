function [ SRP_stack ] = calc_SRPappr(Psi_STFT,  omega, T, N_mm, N_aux, Delta_t_i)
% [ SRP_stack, nT, xi_mm_samp ] = calc_SRPappr(Psi_STFT,  omega, T, N_mm, N_aux, Delta_t_i)
% computes SRP approximation. 
%
% IN:
% Psi_STFT        FD GCCs - frequency x frames x microphone pairs
% omega           frequency vector
% T               sampling period
% N_mm            number of samples inside TDOA interval
% N_aux           bumber of auxilary samples outside TDOA interval (if vector, results are computed for each element)
% Delta_t_i       TDOA vector
%
% OUT:
% SRP_stack       SRP map - frames x candidate locations

[J, P] = size(Delta_t_i);
[K, L, P] = size(Psi_STFT);

SRP_stack = zeros(L, J, length(N_aux));
nT = cell(P);
xi_mm_samp = cell(L,P);
for l = 1:L
    
    SRP = zeros(J, length(N_aux));
    xi_mm_int = zeros(J, length(N_aux));
    
    for p = 1:P
                
        psi = Psi_STFT(:,l,p);
        
        % sample points
        N_half = N_mm(p) + max(N_aux);
        n = -N_half:1:N_half;
        nT{p} = (n*T).';
        
        % samples
        xi_mm_samp{l,p} = real(exp(1j*nT{p}*omega.')*psi); % (10); nDT{p}*omega.' is N x K matrix, psi is K x 1 vector
        
        % interpolate for different values of N_aux
        for N_aux_ind = 1:length(N_aux)
            N_offset = max(N_aux)-N_aux(N_aux_ind);
            xi_mm_int(:,N_aux_ind) = sinc(Delta_t_i(:,p)/T - n(N_offset+1:end-N_offset))*xi_mm_samp{l,p}(N_offset+1:end-N_offset); % (14); Delta_t(:,p)/T - n is J x N matrix, xi_mm_samp is N x 1 vector
        end
        SRP = SRP + 2*xi_mm_int; % (11)
        
    end
    SRP_stack(l,:,:) = SRP;
end