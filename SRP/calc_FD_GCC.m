function [ Psi_STFT ] = calc_FD_GCC(y_STFT)
% [ Psi_STFT ] = calc_FD_GCC(y_STFT)
% computes frequency domain GCCs.
%
% IN:
% y_STFT         microphone signal - frequency x microphones
%
% OUT:
% % Psi_STFT     FD GCCs - frequency x frames x microphone pairs


[K,L,M] = size(y_STFT);
P = M*(M-1)/2;
Psi_STFT = zeros(K,L,P);

for k = 1:K
    
    for l = 1:L
        
        y = squeeze(y_STFT(k,l,:));

        psi=zeros(P,1); 
        p = 0;
        for mprime = 1:M
            for m = mprime+1:M
                p = p+1;
                psi(p) = y(m)*y(mprime)';
            end
        end
       
        psi = psi./(abs(psi)+ 1e-9);
        Psi_STFT(k,l,:) = psi;
        
    end
end

end