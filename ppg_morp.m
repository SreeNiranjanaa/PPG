function [morph] = ppg_morp(ppg,Fs)
% ppg_morp computes the indices and the amplitudes of 3 points in PPG
%   1.) upslope point of PPG
%   2.) Dicrotic notch of PPG
%   3.) Diastolic Peak of PPG

%% INPUT
% ppg_sig - input ppg signal
% Fs - Sampling Frequency

%% OUTPUT
% morph  -struct output
%   pts_upslope - amp_upslope
%   pts_dicro   - amp_dicro
%   pts_dias    - amp_dias
%% STEP 1: PPG PEAK AND FOOT DETECTION
[footin,footamp,pin,pamp] = ppg_footdetect(ppg,Fs);

%% STEP 2 : 1st and 2nd Deivative of PPG
% First Derivative
ppg_d = derivative(ppg,1);

% Second Derivative
ppg_a = derivative(ppg,2);

% Filtering of Acceleration PPG
[z1,p1] = butter(2,0.09,'low');
ppg_acc1 = filtfilt(z1,p1,ppg_a);
ppg_acc = [ppg_acc1;zeros(30,1)];

%% STEP 3: UPslope point Indices and Amplitude
for i = 1:length(pin)
    rng1 = [footin(i):pin(i)];
    [Upslope] = ppg_d(rng1);
    [max1(i),indmax1(i)] = max(Upslope);
    
    ind_dmax(i) = rng1(indmax1(i));
    ind_ppgup(i) = ind_dmax(i);
    ppg_up(i)  = ppg(ind_dmax(i));
end

pts_upslope = ind_ppgup;
der_upslope = max1;
amp_upslope = ppg_up;

%% STEP 4: Dicrotic Notch indices and Amplitude

for j = 1:length(pin)
    rng4 = [pin(j):pin(j)+25];
    [Dicro] = ppg_acc(rng4);
    [max4(j),indmax4(j)] = max(Dicro);
    
    ind_adicro(j) = rng4(indmax4(j));
    ind_dicro(j) = ind_adicro(j);
    ppg_dicro(j) = ppg(ind_adicro(j));
end

pts_dicro = ind_dicro;
acc_dicro = max4;
amp_dicro = ppg_dicro;

%% STEP 5 :  Diastolic Peak indices and amplitudes
% Maximum of Inverse of Acceleration PPG
inv_ppgacc1 = 1.01*max(ppg_acc)-ppg_acc;
inv_ppgacc = [inv_ppgacc1;zeros(20,1)];

for j = 1:length(pts_dicro)
     rng5 = [pts_dicro(j):pts_dicro(j)+20];
    [Dias] = inv_ppgacc(rng5);
    [~,indmax5(j)] = max(Dias);
    
    ind_adias(j) = rng5(indmax5(j)); 
    ppg_adias(j) = ppg_acc(ind_adias(j));
    ind_dias(j)  = ind_adias(j);
    ppg_dias(j)  = ppg(ind_adias(j));
end

pts_dias = ind_dias;
acc_dias = ppg_adias;
amp_dias = ppg_dias;

morph = struct('pts_base',footin(:),'amp_base',footamp(:),...............
        'pts_up',pts_upslope(:),'amp_up',amp_upslope(:),.................
        'pts_peak',pin(:),'amp_peak',pamp(:),.................
        'pts_dicro',pts_dicro(:),'amp_dicro',amp_dicro(:), ..............
        'pts_dias',pts_dias(:),'amp_dias',amp_dias(:));
end

