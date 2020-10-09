function [ppg_sig,morph] = ppg_morp1000(ppg,Fs)
% ppg_morp computes the indices and the amplitudes of 3 points in PPG
%   1.) Base point of PPG
%   2.) upslope point of PPG
%   3.) PEak point of PPG
%   4.) Dicrotic notch of PPG
%   5.) Diastolic Peak of PPG

%% INPUT
% ppg_sig - input resampled ppg signal  (1kHz)
% Fs - Sampling Frequency

%% OUTPUT
% morph  -struct output
%   pts_base    - amp_base
%   pts_upslope - amp_upslope
%   pts_peak    - amp_peak
%   pts_dicro   - amp_dicro
%   pts_dias    - amp_dias

% Butterworth filter design
[b,a] = butter(4,0.009,'low');
% Filtering using zero-phase filtering
ppg_sig = filter(b,a,ppg);

%% STEP 1 :PEAK Detection
[pamp,pin] = ppg_pkdetect1000(ppg_sig,Fs);

% Removing the extra points
diff1 = diff(pin);       
mean1 = mean(diff1);  std1 = std(diff1);
I1 = find(diff1 < 0.8*(mean1-std1));
pin(I1+1)=[];   pamp(I1+1)=[];
x = find(diff(pin)<300);
pin(x+1)=[];
pamp(x+1)=[];


%% STEP 2 : Base point of PPG
thresh = 250;
  % To Find the PPG Base point with Interval range of Two ECG Peaks
if pin(1)> thresh
    for i = 1:length(pin)
        rng2 = [pin(i)-thresh:pin(i)];
        [Base] = ppg_sig(rng2);
        [max2(i),indmax1(i)] = min(Base);

        ind_abase(i) = rng2(indmax1(i));
        ppg_abase(i) = max2(i);
        ind_base(i)  = ind_abase(i);
        ppg_base(i)  = ppg_sig(ind_abase(i));
    end
else
    for i = 1:length(pin)-1
        rng2 = [pin(i+1)-thresh:pin(i+1)];
        [Base] = ppg_sig(rng2);
        [max2(i),indmax1(i)] = min(Base);

        ind_abase(i) = rng2(indmax1(i));
        ppg_abase(i) = max2(i);
        ind_base(i)  = ind_abase(i);
        ppg_base(i)  = ppg_sig(ind_abase(i));
    end 
end

    footin = [ind_base(:)];
    footacc = [ppg_abase(:)];
    footamp = [ppg_base(:)];
    
    pts_base = footin;
    amp_base = footamp;
    
    % Removing the last peaks values
    pin(end)=[]; pamp(end)=[]; 
    % Removing the peak points before the ist base point
    condition = find(pin < footin(1));
    pin(condition)=[]; pamp(condition) = []; 
    pin=pin(:);  pamp = pamp(:);
    pts_peak = pin;   amp_peak = pamp;
    
%% STEP 3 : UPSLOPE of PPG
ppg_d1 = derivative(ppg_sig,1);
for i = 1:length(pin)
    rng1 = [footin(i):pin(i)];
    [Upslope] = ppg_d1(rng1);
    [max1(i),indmax1(i)] = max(Upslope);
    
    ind_dmax(i) = rng1(indmax1(i));
    ind_ppgup(i) = ind_dmax(i);
    ppg_up(i)  = ppg_sig(ind_dmax(i));
end

pts_upslope = ind_ppgup(:);
der_upslope = max1(:);
amp_upslope = ppg_up(:);

%% STEP 4 : Dicrotic Notch of PPG

ppg_ac = derivative(ppg_sig,2);
ppg_d2 = [ppg_ac;zeros(70,1)];
for j = 1:length(pin)-1
    rng4 = [pin(j):pin(j)+250];
    [Dicro] = ppg_d2(rng4);
    [max4(j),indmax4(j)] = max(Dicro);
    
    ind_adicro(j) = rng4(indmax4(j));
    ind_dicro(j) = ind_adicro(j);
    ppg_dicro(j) = ppg_sig(ind_adicro(j));
end

pts_dicro = ind_dicro(:);
acc_dicro = max4(:);
amp_dicro = ppg_dicro(:);

%% STEP 5: Diastolic Peak of PPG
% Maximum of Inverse of Acceleration PPG
inv_ppgacc1 = 1.01*max(ppg_d2)-ppg_d2;
inv_ppgacc = [inv_ppgacc1;zeros(20,1)];

for j = 1:length(pts_dicro)
     rng5 = [pts_dicro(j):pts_dicro(j)+170];
    [Dias] = inv_ppgacc(rng5);
    [max5(j),indmax5(j)] = max(Dias);
    
    ind_adias(j) = rng5(indmax5(j)); 
    ppg_adias(j) = ppg_d2(ind_adias(j));
    ind_dias(j)  = ind_adias(j);
    ppg_dias(j)  = ppg_sig(ind_adias(j));
end

pts_dias = ind_dias(:);
acc_dias = ppg_adias(:);
amp_dias = ppg_dias(:);

morph = struct('pts_base',pts_base(1:end-1,:),'amp_base',amp_base(1:end-1,:),...............
        'pts_up',pts_upslope(1:end-1,:),'amp_up',amp_upslope(1:end-1,:),.................
        'pts_peak',pts_peak(1:end-1,:),'amp_peak',amp_peak(1:end-1,:),.................
        'pts_dicro',pts_dicro,'amp_dicro',amp_dicro, ..............
        'pts_dias',pts_dias,'amp_dias',amp_dias);

end

