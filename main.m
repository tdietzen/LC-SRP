%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Thomas Dietzen
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt).
%
% A Python version of this code is available at
% https://github.com/bilgesu13/LC-SRP-python.
%
% If you find it useful, please cite:
%
% [1] T. Dietzen, E. De Sena, and T. van Waterschoot, "Low-Complexity
% Steered Response Power Mapping based on Nyquist-Shannon Sampling," in
% Proc. 2021 IEEE Workshop Appl. Signal Process. Audio, Acoust. (WASPAA
% 2021), New Paltz, NY, USA, Oct. 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Example of the low-complexity SRP approximation described in [1]. The
% code in main.m defines the configuration, loads microphone signals,
% performs conventional and low-complexity SRP and computes approximation
% and localization errors. This is done for eight different source loations
% (32 frames each).


%% PREAMBLE

clear;
cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));


%% CONFIGURATION

%%% ACOUSTIC SETUP
% speed of sound
c = 340;
% sample rate
fs = 16000;
% bandlimit
w_0 = pi*fs;
% SNR in dB
SNR = 6;


%%% MICROPHONE ARRAY
% circular array, 10cm radius, six microphones 
tmp = load('coord_mic_array.mat');
% array center
arrayCenterPos = tmp.arrayCenterPos;
% microphone positions
micPos = tmp.micPos;
% number of microphones
M = size(micPos,1);


%%% SOURCE LOCATIONS
% 8 different locations
tmp = load('coord_loc_1_8.mat');
true_loc = tmp.true_loc;
% compute ground truth DOA vectors for source locations
true_DOAvec = calc_DOA(true_loc, arrayCenterPos);
% number of processed frames per location
L = 32;


%%% STFT PARAMETERS
% window size
N_STFT = 2048;
% shift
R_STFT = N_STFT/2;
% window
win = sqrt(hann(N_STFT,'periodic'));
N_STFT_half = floor(N_STFT/2)+1;
% frequency vector
omega = 2*pi*linspace(0,fs/2,N_STFT_half).';


%%% CANDIDATE LOCATIONS
% polar angles of candidate locations
ang_pol = 90:2:180;
% azimuth angles of candidate locations 
ang_az = 0:2:358;
% compute candidate DOA vectors and TDOAs
[DOAvec_i, Delta_t_i] = gen_searchGrid(micPos, ang_pol, ang_az, 'spherical', c);


%%% SRP APPROXIMATION PARAMETERS
% compute sampling period and number of samples within TDOA interval
[ T, N_mm ] =calc_sampleParam(micPos, w_0, c);
% number of auxilary samples (approximation will be computed for all values in vector)
N_aux = 0:2;


%% PROCESSING

% init results (per source location, frame, number of auxilary samples)
% approximation error in dB
res.approxErr_dB = zeros(size(true_loc,1), L, length(N_aux));
% localization error 
res.locErr = zeros(size(true_loc,1), L, length(N_aux)+1); % res.locErr(:,:,1) refers to  conventional SRP

for true_loc_idx = 1:size(true_loc,1)
    
    disp(['PROCESSING SOURCE LOCATION ' num2str(true_loc_idx)])
    %%% GENERATE MICROPHONE SIGNALS
    
    % speech component for selected source     
    x_TD = audioread(['x_loc' num2str(true_loc_idx) '.wav']);
    % noise component
    v_TD = audioread('v.wav');
    % scale noise component
    v_TD = set_SNR(x_TD, v_TD, SNR);
    
    % transform to STFT domain
    x_STFT  = calc_STFT(x_TD, fs, win, N_STFT, R_STFT, 'onesided');
    v_STFT  = calc_STFT(v_TD, fs, win, N_STFT, R_STFT, 'onesided');
    
    % discard frames that do not contain speech energy (local SNR 15dB below average)
    l = 1;
    useframe_idx = [];
    while length(useframe_idx) < L
        SNR_local = pow2db(sum(abs(x_STFT(:,l,1)).^2)/sum(abs(v_STFT(:,l,1)).^2));
        if SNR_local > SNR - 15
            useframe_idx = [useframe_idx, l];
        end
        l = l + 1;
    end
        
    % final nicrophone signal in STFT domain
    y_STFT = x_STFT(:,useframe_idx,:) + v_STFT(:,useframe_idx,:);
    
        
    %%% PROCESSING
    
    psi_STFT =calc_FD_GCC(y_STFT);
    
    % conventional SRP
    disp('* compute conventional SRP (stay tuned, this may take a minute)...')
    tic;
    SRP_conv = calc_SRPconv(psi_STFT, omega, Delta_t_i);
    toc;
    
    % SRP approximation based on shannon nyquist sampes
    disp('* compute SRP approximation...')
    tic;
    SRP_appr = calc_SRPappr(psi_STFT, omega, T, N_mm, N_aux, Delta_t_i);
    toc;
    
    % init errors
    approxErr_dB = zeros(L, length(N_aux));
    locErr = zeros(L, length(N_aux)+1); % res.locErr(:,1) refers to conventional SRP
    
    % compute localization error for conventional SRP
    [~, maxIdx_conv] = max(SRP_conv, [], 2);
    estim_DOAvec = DOAvec_i(maxIdx_conv,:);
    locErr(:,1) = rad2deg(acos(...
        (estim_DOAvec*transpose(true_DOAvec(true_loc_idx,:)))./(sqrt(sum(estim_DOAvec.^2, 2))*norm(true_DOAvec(true_loc_idx,:)))...
        ));
    
    % compute approximation and localization errors for SRP approximation
    for N_aux_ind = 1:length(N_aux)
        %
        % approximation error
        approxErr = sum((SRP_conv - SRP_appr(:, :, N_aux_ind)).^2, 2)./sum(SRP_conv.^2, 2);
        approxErr_dB(:, N_aux_ind) = pow2db(approxErr);
        %
        % localization error
        [~, maxIdx_conv] = max(SRP_appr(:,:,N_aux_ind), [], 2);
        estim_DOAvec = DOAvec_i(maxIdx_conv,:);
        locErr(:,N_aux_ind+1) = rad2deg(acos(...
            (estim_DOAvec*transpose(true_DOAvec(true_loc_idx,:)))./(sqrt(sum(estim_DOAvec.^2, 2))*norm(true_DOAvec(true_loc_idx,:)))...
            ));
    end
    
    
    %% SAVE
    
    res.approxErr_dB(true_loc_idx, :, :) = approxErr_dB;
    res.locErr(true_loc_idx, :, :) = locErr;
        
    
end

disp('DONE.')
