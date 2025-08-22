%% Grasp_Force_Testing: interpData
%% Author: Alexander B. Ambrose
%% Date: 8-12-2025

%% Description:
%  This script syncronizes the time series data from the finger controller
%  and the F/T sensor. The data stream from the finger has a much lower
%  sample rate than the F/T sensor, so the time series data are made the 
%  same length by interpolation.
%  
%  Outputs:
%  Matched_Data.mat -> Strucutre of trials with the finger data and the 
%  6-axis force and torque readings from the F/T sensor with time stamps.

clear all; close all; clc;

load('Data/Finger_Data.mat'); % load the unpacked teensy data. produces a struct called 'finger'

load('Data/Rokubi_Data.mat'); % load the unpacked rokubi. produces a struct called 'rokubi'

fieldNames = fieldnames(finger);
kk = 1;
for ii = 1:numel(fieldNames)
    tic;
    Name = fieldNames{ii};
    fin{ii} = finger.(Name);
    rok{ii} = rokubi.(Name);
    ff = 1/round(mean(diff(fin{ii}.time)), 3);
    fr = 1/round(mean(diff(rok{ii}.device_timestamp/1000000)), 3);
    num = ceil(ff/2); % 0.5 second moving average filter
    coeff = ones(1, num)/num;
    fin{ii}.F = filter(coeff, 1, fin{ii}.F);
    num = ceil(fr/2); % 0.5 second moving average filter
    coeff = ones(1, num)/num;
    rok{ii}.Fx = filter(coeff, 1, rok{ii}.Fx);
    rok{ii}.Fy = filter(coeff, 1, rok{ii}.Fy);
    rok{ii}.Fz = filter(coeff, 1, rok{ii}.Fz);
    rok{ii} = downsample(rok{ii}, fr/ff);
    offset = 105;
    stop = 415;
    idf = find(fin{ii}.F(offset:stop) == max(fin{ii}.F(offset:stop)), 1, 'first')+offset;
    fin{ii} = fin{ii}(idf:end, :); fin{ii}.time = fin{ii}.time - fin{ii}.time(1);
    stop = 543;
    idr = find(rok{ii}.Fx(1:stop) == max(rok{ii}.Fx(1:stop)), 1, 'first');
    rok{ii} = rok{ii}(idr:end, :); 
    rok{ii}.device_timestamp = (rok{ii}.device_timestamp - rok{ii}.device_timestamp(1))/1E6;
    rok{ii} = removevars(rok{ii}, ["host_time_s", "status"]);
    rok{ii}.Properties.VariableNames{1} = 'time';

    dfx = diff(rok{ii}.Fx);
    offset = 3600;
    [npks, nlocs] = findpeaks(-dfx(offset:end), ...
          'MinPeakDistance', 300, 'MinPeakHeight', 0.01);
    ix(ii) = nlocs(end)+offset;
    idx(ii) = find(abs(dfx(1:ix(ii))) < 0.001, 1, 'last');
    rok{ii} = rok{ii}(1:idx(ii), :);
    Temperature_C = NaN*ones(size(fin{ii}.time));
    tmp = fin{ii}; tmp = addvars(tmp, Temperature_C);
    for jj = 1:numel(tmp.Properties.VariableNames)
        xq = linspace(0, fin{ii}.time(end), size(fin{ii}, 1))';
        x = linspace(0, fin{ii}.time(end), size(rok{ii}, 1))';
        tmp{:, jj} = interp1(x, rok{ii}{:, jj}, xq);
    end
    tmp.Properties.VariableNames{1} = 'time';
    tmp.Properties.VariableNames{2} = 'Fx';
    tmp.Properties.VariableNames{3} = 'Fy';
    tmp.Properties.VariableNames{4} = 'Fz';
    tmp.Properties.VariableNames{5} = 'Mx';
    tmp.Properties.VariableNames{6} = 'My';
    tmp.Properties.VariableNames{7} = 'Mz';
    tmp.Properties.VariableNames{8} = 'Temperature_C';
    rok{ii} = tmp;
    rok{ii}.time = fin{ii}.time;
    clear tmp;
    % plot(fin{ii}.F)
    % hold on
    % plot(rok{ii}.Fx)
    % hold off
    fprintf('-> Done Unpacking Trial %d in %.2f s\n', ii, toc);
end

save('Data/Matched_Data.mat', 'fin', 'rok');