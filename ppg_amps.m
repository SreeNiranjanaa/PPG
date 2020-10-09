function [sysamp,diasamp] = ppg_amps(ppg_sig,Fs)
% ppg_amps computes the photoplethysmographic systolic and 
% diastolic amplitude
% sysamp is the systolic peak amp - foot point amp
% diasamp is the diastolic peak amp - foot point amp

%% INPUT
% ppg_sig - input ppg signal
% Fs - Sampling Frequency
%% OUTPUT
% sysamp - systolic amplitude
% diasamp - diastolic amplitude

%% STEP 1 : PEAK, FOOT and DIASTOLIC POINT DETECTION
if Fs ==125
    morph  = ppg_morp(ppg_sig,Fs);
else if Fs == 1000
        morph = ppg_morp1000(ppg_sig,Fs);
    end
end
%% STEP 2 : PPG SYSTOLIC  & DIASTOLIC AMPLITUDE
sysamp = morph.amp_peak-morph.amp_base;
diasamp = morph.amp_dias-morph.amp_base;
end
