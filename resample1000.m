function [y,Ny] = resample1000(x,Nx,Type)
% resample1000 resamples the given signal of sampling rate Fs to sampling 
% rate of 1000
% INPUT  : x     - input signal
%          Nx    - Length of Signal
%          Type  - Type of Database
%           1    - MIMIC-II Database
%           2    - Queensland University
%           3    - Capnobase Database

% OUTPUT : y - resampled signal
%          Ny - No. of Samples

if size(x,2) == 1
    x = x';
end
% Resampling Parameters
xpad = 15;                              %  padding length
x_zeropad = [repmat(x(1),1,xpad), x , repmat(x(end),1,xpad)];
x_zeropad = x_zeropad(:);
switch(Type)
    %  Factor of 8 - 125 ---> 1000
    case 1      
        factor = 8;
        y = resample(x_zeropad,factor,1);
        % Removing padding elements
        y(1:(xpad*factor))=[];  
        y((Nx*factor)+(1:(xpad*factor)))=[]; 
        Ny = length(y);
    %  Factor of 10 - 100 ---> 1000
    case 2
        factor = 10;
        y = resample(x_zeropad,factor,1);
        % Removing padding elements
        y(1:(xpad*factor))=[];  
        y((Nx*factor)+(1:(xpad*factor)))=[]; 
        Ny = length(y);
    %  Factor of 10/3 - 300 ---> 1000
    case 3
        y = resample(x_zeropad,10,3);
        % Removing padding elements  Geneates 50 values on each side
        y(1:50)=[];  
        y(end-49:end)=[]; 
        Ny = length(y);
end


end

