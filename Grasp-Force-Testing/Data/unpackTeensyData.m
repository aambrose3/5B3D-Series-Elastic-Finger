%% Grasp_Force_Testing: unpackTeensyData
%% Author: Alexander B. Ambrose
%% Date: 8-12-2025

%% Description:
%  This script unpacks the data streams from the mid-level motor
%  controllers for a single finger of the gripper.
%
%  Outputs:
%  Finger_Data.mat -> Strucutre of trials with reports on the finger states
%  from a single finger of the gripper

clear all; close all; clc;

%% Import Teensy Data (.log files from PuTTY)
TData.t1 = readtable('T1', 'FileType', 'text');
TData.t2 = readtable('T2', 'FileType', 'text');
TData.t3 = readtable('T3', 'FileType', 'text');
TData.t4 = readtable('T4', 'FileType', 'text');
TData.t5 = readtable('T5', 'FileType', 'text');
TData.t6 = readtable('T6', 'FileType', 'text');
TData.t7 = readtable('T7', 'FileType', 'text');
TData.t8 = readtable('T8', 'FileType', 'text');
TData.t9 = readtable('T9', 'FileType', 'text');
TData.t10 = readtable('T10', 'FileType', 'text');

fieldNames = fieldnames(TData);
kk = 1;
for ii = 1:numel(fieldNames)
    tic;
    Name = fieldNames{ii};
    data = TData.(Name);
    for jj = 1:size(data, 1)
        if jj == 1
            x.time(jj, 1) = 0.03; % Relative time in seconds
            x.t1(jj, 1) = data{1, 2}; % \theta_1 of the finger (rad)
            x.t2(jj, 1) = data{2, 2}; % \theta_2 of the finger (rad)
            x.t3(jj, 1) = data{3, 2}; % \theta_A/actuator angle (rad) 
            x.p1(jj, 1) = data{4, 2}; % approximated p_1 (m)
            x.p1(jj, 1) = data{5, 2}; % approximated p_2 (m)
            x.F(jj, 1) = data{6, 2};  % approximated total enveloping force (N)
            kk = kk+1;
        else
            try
                str = cell2mat(data{jj, 1});
                % search for the string 'Reporting' to collect the next
                % frame
                if strcmp(str(1:9), 'Reporting') == 1
                    % the significant digits of time changes
                    try 
                        x.time(kk, 1) = str2num(str(end-9:end-3));
                    catch
                        try 
                            x.time(kk, 1) = str2num(str(end-8:end-3));
                        catch
                            x.time(kk, 1) = str2num(str(end-7:end-3));
                        end
                    end
                end
                x.t1(kk, 1) = data{jj+1, 2};
                x.t2(kk, 1) = data{jj+2, 2};
                x.t3(kk, 1) = data{jj+3, 2};
                x.p1(kk, 1) = data{jj+4, 2};
                x.p2(kk, 1) = data{jj+5, 2};
                x.F(kk, 1) = data{jj+6, 2};
                kk = kk+1;
            end
        end
    end
    % fix if there is an additional NaN in theta 1 data...
    if any(isnan(x.t1)) == 1
        idx = find(isnan(x.t1) == 1, 1, 'first');
        x.t1 = x.t1(1:idx-1);
    end
    x = struct2table(x);
    % assemble the strucutre
    finger.(Name) = x;
    clear x data
    fprintf('-> Done Unpacking Trial %d in %.2f s\n', ii, toc);
    kk = 1;
end
% save the assembled structure to a mat file
save('Finger_Data.mat', 'finger')