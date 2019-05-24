function [peak, typecounter] = parsing_submovement(Vel, Acc)
    Vmax = find(Vel==max(Vel));
    % the location of the max velocity
    under = find(Vel(Vmax:end)<=0.3,1) + (Vmax-1);
    
    % after the max velocity, find the first point of the velocity < 0.3cm/s
    zc = find(Vel(Vmax:end-1).*Vel(Vmax+1:end)<=0) + (Vmax-1);
    % after the max velocity, find the zero crossing of the velocity
    peak1=[]; % type I counter
    peak2=[]; % type II counter
    peak3=[]; % type III counter
    peak4=[]; % type IV counter
    %% calculating
    
    peak1 = func_elliott(Acc, 1, Vmax, 1); % Type I
    if length(zc) < 1
        % in this case, there was no zero crossing occured
        % it would exist type II or III
        peak3 = func_elliott(Vel,under,length(Vel),3); % find type III
        if ~isempty(peak3) % Type III exists
            peak2 = func_elliott(Acc,Vmax,under,2);
        else % Type II only
            peak2 = func_elliott(Acc,Vmax,length(Acc),2);
        end
    else
        % zero crossing was observed in this trial
        % 3 sub-case may occur
        if sum(Vel(under:zc(1))>=5)>0 % case 1 : Type II >> III >> IV
            peak2 = func_elliott(Acc,Vmax,under,2);
            peak3 = func_elliott(Vel,under,zc(1),3);
            peak4 = func_elliott(Vel,zc(1),length(Vel),4);
            %                         subplot(3,1,1); hold on; plot(peak4,Vel(peak4),'^')
        else
            peak2 = func_elliott(Acc,Vmax,zc(1),2);
            if length(zc)<2 % case 3 : Type II >> IV
                peak4 = func_elliott(Vel,zc(1),length(Vel),4);
            else % case 2 : Type II >> IV >> III. There must be 2 zero crossing
                peak4 = func_elliott(Vel,zc(1),zc(2),4);
                peak3 = func_elliott(Vel,zc(2),length(Vel),3);
            end
        end
    end
    %% concatenate different peaks
    peak = [peak1;peak2;peak3;peak4];
    if isempty(peak)~=1
        % sort the peak
        [~,index] = sort(peak(:,1));
        peak = peak(index,:);
    end
    %% which types of peaks were detected
    if size(peak,1) == 0;
        typecounter = [1 0 0 0 0];
    else
        typecounter = [0 size(peak1,1) size(peak2,1) size(peak3,1) size(peak4,1)];
        typecounter = ~isnan(typecounter./typecounter);
    end
end

function [peaks] = func_elliott(x,initial,final,type)
% it's a function of elliott's method
% [loc,pks,category] = elliott(x,initial,final,type)
% where x is a time series acc/vel
% give the initial and final point to define the boundary of x
% most of all, we need the 'type' input to go different algorithm
temporal = 72 * 130 / 1000;
% temporal = 0;
peaks = [];
if type==1
    x = [0; x];
    [primpeak, primloc] = findpeaks(x, 'sortstr', 'descend', 'npeak', 1); % peak acc of primary movement
    [p1, locp1] = findpeaks(x(initial:final)); % Type I : acceleration
    dx = diff(x);
    p1(locp1 == find(dx(1:end-1).*dx(2:end)<=0,1)+1)=[]; % delete primary peak
    locp1(locp1 == find(dx(1:end-1).*dx(2:end)<=0,1)+1)=[]; % delete primary peak
    
    for i =1:length(locp1)
        if length(x(locp1(i):-1:initial))>=3
            [v1, locv1] = findpeaks(-x(locp1(i):-1:initial),'npeaks',1);
            if isempty(v1)
                v1 = x(initial);
                locv1 = initial;
            else
                v1 = -v1;
                locv1 = locp1(i) - locv1 + 1;
            end
            % get duration
            dur_1 = find(x(locv1-1:-1:initial, 1) >= p1(i), 1);
            dur_2 = find(dx(locv1-1:-1:2) .* dx(locv1-2:-1:1) <= 0, 1) ;
            if isempty([dur_1, dur_2])
                dur = 0;
            else
                dur = min([dur_1, dur_2]) + locp1(i) - locv1;
            end
            % get amplitude
            p1(i) = p1(i) - v1;

            if p1(i)>=0.1*primpeak && dur>temporal  
                peaks = [peaks; locp1(i), locv1, p1(i), dur, 1];
            end
        end
    end
elseif type == 2
    [primpeak, primloc] = findpeaks(x, 'sortstr', 'descend', 'npeak', 1); % peak acc of primary movement
    [p2, locp2] = findpeaks(-x(initial:final)); % Type II : deceleration
    dx = diff(x);
    locp2 = initial + locp2; % find all valley
    
    for i = 1:length(locp2)
        if length(x(locp2(i):final))>=3
            [v2,locv2] = findpeaks(x(locp2(i):final),'npeaks',1);
            if isempty(v2)
                v2 = x(final);
                locv2 = final;
            else
                v2 = -v2;
                locv2 = locp2(i) + locv2;
            end
            % get duration
            dur_1 = find(-x(locv2:final, 1) >= p2(i), 1);
            dur_2 = find(dx(locv2:final-2) .* dx(locv2+1:final-1) <= 0, 1) + 1;
            if isempty([dur_1, dur_2])
                dur = 0;
            else
                dur = min([dur_1, dur_2]) + locv2 - locp2(i);
            end
            % get amplitude
            p2(i) = p2(i) - v2;
            % elliott's condition for Type II
            if p2(i)>=0.1*primpeak && dur>temporal  
                peaks = [peaks; locp2(i), locv2, p2(i), dur, 2];
            end
        end
    end
elseif type==3
    loc3 = find(x(initial:final)>=5,1) + initial - 1;
    if isempty(loc3)~=1
        p3 = max(x(initial:final));
        locp3 = find(x(initial:final) == p3) + initial - 1;
        dur = final - initial;
        peaks = [initial, locp3, p3, dur, 3];
    end
elseif type==4
    loc4 = find(x(initial:final)<=-1,1) + initial-1;
    if isempty(loc4)~=1
        p4 = min(x(initial:final));
        locp4 = find(x(initial:final) == p4) + initial - 1;
        dur = final - initial;
        peaks = [initial, locp4, p4, dur, 4];
    end
else
    return
end
end
