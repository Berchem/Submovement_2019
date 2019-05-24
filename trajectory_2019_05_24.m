close all; clear all; clc; tic;
% We would use the following scripts to parse the experimental data:
%     parsing_path 
%     parsing_data 
%     parsing_trial
%
% And some defined functions:
%     func_filtmat_class
%     func_start_stop_filter
%     func_elliott
%
% Adjustment:
% measures{sub,1}{dist,cond}(trial,4)
%     original: if under-shoot existed
%                   the point at vel < 0.3
%               else
%                   the end of movement
%     adjusted: the end of movement
%
% parsing_data
% parameter:
%            index of subject,   integer / double
%            index of distance,  integer / double
%            index of condition, integer / double
%            index of trial,     integer / 'all'
%            filter method,      1 for low-pass
%                                2 for kalman filter
%                                3 for weighted moving average
%            initial_path,       set the initial path to data folder
%
% return:
%         Kinematic{dist, cond}{trial, sub} = [Pos, Vel, Acc, Jerk]
%         Submoves{dist, cond}{trial, sub} = infomation of peaks
%         Counter{sub, 1}{dist, cond}(trial, :) = [I, II, III, IV]
%         Measures{sub, 1}{dist, cond}(trial, :) = [peak acc, peak vel, peak negative acc, the end of the movement]
%         Amount{dist,cond}(sub,:) = [cumulative amount of submovement s0:s20]

index_sub = 1:12; % the indices of subjects
index_dist = 1:3; % the distances
index_cond = 1:5; % condition
index_trial = 'all';
algorithm = 'elliott';
initial_path = 'E:\Data\exp1\12 subjects'

[Kinematic, Submoves, Counter, Measures, Amount] = parsing_data(...
    index_sub, index_dist, index_cond, index_trial, algorithm, initial_path);

for sub = index_sub
    for dist = index_dist
        for cond = index_cond
            Type_avg{dist,cond}(sub,:) = mean(Counter{sub, 1}{dist, cond});
            Various{dist, cond}(sub, :) = std(Measures{sub, 1}{dist, cond});
            for trial = 1:100
                peak = Submoves{dist, cond}{trial, sub};
                if ~isempty(peak)
                    kinematic = Kinematic{dist, cond}{trial, sub};
                    Position{dist, cond}{trial, sub} = kinematic(peak(:, 1), 1);
                end
            end
        end
    end
end

toc

% checking_graph;