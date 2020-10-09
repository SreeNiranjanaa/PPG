function [ pir ] = pir_value(ppg_sig,Fs)
%pir_value computes the photoplethysmographic intensity ratio
% PIR is the ratio of intensity at peak to intensity 
% in the foot point  (pkamp/footamp)

%% INPUT
% ppg_sig - input ppg signal
% Fs - Sampling Frequency
%% OUTPUT
% pir - PIR value

%% STEP 1 : PEAK AND FOOT POINT DETECTION
[footin,footamp,pin,pamp] = ppg_footdetect(ppg_sig,Fs);
footamp = footamp(:);
pamp = pamp(:);
%% STEP 3 : DETERMINATION OF PHOTOPLETHYSMOGRAM INTENSITY RATIO
pir = pamp./footamp;

end
