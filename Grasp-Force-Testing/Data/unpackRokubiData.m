%% Grasp_Force_Testing: unpackRokubiData
%% Author: Alexander B. Ambrose
%% Date: 8-12-2025

%% Description:
%  This script unpacks the data streams from the 6-axis BOTA Systems Rokubi 
%  Serial Force-Torque (F/T) sensor. The sensor communicates with Python on a
%  host PC and this data is stored in a csv file with time stamps
%  
%  Outputs:
%  Rokubi_Data.mat -> Strucutre of trials with the 6-axis force and torque
%  readings from the F/T sensor with time stamps.

clear all; close all; clc;

%% Import Rokubi Data
RData.t1 = readtable('T1.csv', 'FileType', 'text');
RData.t2 = readtable('T2.csv', 'FileType', 'text');
RData.t3 = readtable('T3.csv', 'FileType', 'text');
RData.t4 = readtable('T4.csv', 'FileType', 'text');
RData.t5 = readtable('T5.csv', 'FileType', 'text');
RData.t6 = readtable('T6.csv', 'FileType', 'text');
RData.t7 = readtable('T7.csv', 'FileType', 'text');
RData.t8 = readtable('T8.csv', 'FileType', 'text');
RData.t9 = readtable('T9.csv', 'FileType', 'text');
RData.t10 = readtable('T10.csv', 'FileType', 'text');

fieldNames = fieldnames(RData);
kk = 1;
for ii = 1:numel(fieldNames)
    tic;
    Name = fieldNames{ii};
    data = RData.(Name);
    % determine the average sample rate of the data
    fr = 1/round(mean(diff(data.device_timestamp/1000000)), 3);
    num = ceil(fr/4); % 0.26 second moving average filter
    coeff = ones(1, num)/num;
    df = diff(filter(coeff, 1, data.Fx)); % filter the data
    offset = 400; % find when the sensor makes contact with the finger
    [pks, idx] = findpeaks(df(offset:end), 'NPeaks', 1, 'SortStr', 'descend');
    idx=idx-fr/2;
    % Zero the sensor out
    Zfx = mean(data.Fx(1:idx));
    Zfy = mean(data.Fy(1:idx));
    Zfz = mean(data.Fz(1:idx));
    Zmx = mean(data.Mx(1:idx));
    Zmy = mean(data.My(1:idx));
    Zmz = mean(data.Mz(1:idx));
    % assemble the structure
    data.Fx = abs(data.Fx - Zfx);
    data.Fy = abs(data.Fy - Zfy);
    data.Fz = abs(data.Fz - Zfz);
    data.Mx = data.Mx - Zmx;
    data.My = data.My - Zmy;
    data.Mz = data.Mz - Zmz;

    rokubi.(Name) = data;
    clear x data
    fprintf('-> Done Unpacking Trial %d in %.2f s\n', ii, toc);
end
% save the assembled structure to a mat file
save('Rokubi_Data.mat', 'rokubi')