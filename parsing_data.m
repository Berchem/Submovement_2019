function [Kinematic, Submoves, Counter, Measures, Amount] = parsing_data(index_sub, index_dist, index_cond, index_trial, algorithm, initial_path)
    condi = {'Fast','Fast-Mid','Middle','Mid-Accurate','Accurate'};

    path_data = uigetdir(initial_path);
    
    for sub = index_sub
        for dist = index_dist
            for cond = index_cond
                % parsing the data path
                [~, ~, path_cond] = parsing_path(path_data, sub, dist * 10, cond);
                path_current = [path_cond, '\*.dat'];
                listdir = dir(path_current);
                submoves_ct = zeros(1, 21);
                if strcmp(index_trial, 'all'); index_trial = 1:size(listdir, 1); end
                for trial = index_trial
                    %% load data
                    fname = listdir(trial).name;
                    fname = [path_cond, '\', fname];
                    data = load(fname);
                    %% parsing each trial
                    [kinematic, peaks, typecounter, measures] = parsing_trial(data, algorithm);
                    Kinematic{dist, cond}{trial, sub} = kinematic;
                    Submoves{dist, cond}{trial, sub} = peaks;
                    Counter{sub, 1}{dist, cond}(trial, :) = typecounter;
                    Measures{sub, 1}{dist, cond}(trial, :) = measures;
                    %% count the number of submoves
                    n = size(peaks, 1) + 1;
                    if n > 20; n = 21; end
                    submoves_ct(1, n) = submoves_ct(1, n) + 1;
                end
                Amount{dist,cond}(sub,:) = submoves_ct;
            end
        end
    end
end

function [path_subj, path_dist, path_cond] = parsing_path(path_data, subj, dist, cond)
    path_subj = [path_data, '\S', num2str(subj)];
    path_dist = [path_subj, '\', num2str(dist)];
    path_cond = [path_dist, '\', num2str(cond)];
end

function [kinematic, peaks, typecounter, measures] = parsing_trial(data, algorithm)
    % adjust pixel to cm, 0.00201 is a right trzcform value for pixel to cm
    data(:, 7) = data(:, 7) * 0.00201;
    
    [start, stop] = func_started_stopped_filter(data(:,7));
    data(:, 7) = func_filtmat_class((1/130),5,data(:,7), 1); % 5 = cutoff frequency
    
    % get Position, Velocity, Acceleration, Jerk
    Pos = data(start:stop, 7) - data(start, 7);
    Vel = diff(data(start:stop+1, 7)) * 130;
    Acc = diff(Vel) * 130;
    Jerk = diff(Acc) * 130;
    
    % get kinematic
    kinematic = [Pos-Pos(1), Vel, [0; Acc], [0; 0; Jerk]];
    
    % get measures
    measures(1,1) = Pos(Acc==max(Acc)); % peak acceleration,
    measures(1,2) = Pos(Vel==max(Vel)); % peak velocity,
    measures(1,3) = Pos(Acc==min(Acc)); % peak negative acceleration
    measures(1,4) = Pos(end);           % the end of the movement
    
    % detect sub-movements in each trial
    if strcmp(algorithm, '2019.05')
      [peaks, typecounter] = parsing_submovement2(Vel, Acc);
    else
      [peaks, typecounter] = parsing_submovement(Vel, Acc);
    end
end

function [start, stop] = func_started_stopped_filter(pos)
    start = 0;
    stop = length(pos) - 1;
    ct_index = 0;
    for i = 1:length(pos)-1
        vel = (pos(i+1) - pos(i)) * 130;
        if vel >= 0.3  % velocity > 3 mm/sec
            ct_index = ct_index + 1;
        else
            ct_index = 0;
        end

        ct_time = ct_index / 130;
        if ct_time >= 0.030  % duration > 30 msec
            start = i - ct_index + 1;
            break
        end
    end

    ct_index = 0;
    for j = i:length(pos)-1
        vel = (pos(i + 1) - pos(i)) * 130;
        if abs(vel) <= 0.3
            ct_index = ct_index + 1;
        else
            ct_index = 0;
        end

        ct_time = ct_index / 130;
        if ct_time >= 0.030  % duration > 30 msec
            stop = i - ct_index;
            break
        end
    end
end

function [ fdata ] = func_filtmat_class( dt, cutoff, data, ftype, forder)

if nargin == 3
   ftype = 1;
   forder = 2;
end

if nargin == 4
   forder = 2;
end

cutoff = cutoff / (sqrt(2)-1)^(0.5/forder);

if ftype == 1
   [ b, a] = butter( forder, 2*cutoff*dt);					%  low-pass
elseif ftype == 2
   [ b, a] = butter( forder, 2*cutoff*dt, 'high');	%  high-pass
else
   [ b, a] = butter( forder, 2*cutoff*dt);					%  band-pass
end

[ n_rows, n_cols] = size( data );
fdata = zeros(n_rows, n_cols);

for i=1:n_cols
   fdata( :, i) = filtfilt( b, a, data(:,i) );
end
