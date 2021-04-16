function [processed,removed] = removeSpikes(raw)
% based on removeremoved.m, removes scanner-looking large spikes in the data
spike_thresh = 3;   % in ydiffstd, what we are removing

% Extract x and y components of raw
x = raw(:,2);
y = raw(:,3);

% Look at how y data behaves to see when there are large jumps before
% and after blink
ydiff     = diff(y);
ydiffmean = nanmean(abs(ydiff));
ydiffstd  = nanstd(abs(ydiff));

spikes = (abs(ydiff) > (ydiffmean + spike_thresh*ydiffstd));
spikes  = find(spikes==1);

removed = [];
% Delete data before start of blink
for i = 1:length(spikes)
  s = spikes(i); % Timepoint right before the start of the blink
  go = 20; % Will continue deleting until there are this many consecutive 
           % non-blink datapoints. Since the RecMem experiment records at
           % 1000 Hz
  while s>0 & go>0
    % Check to see if y value shows a significant jump in the data
    if (abs(ydiff(s)) > (ydiffmean + 3*ydiffstd))
      raw(s,2:4)=nan; % Get rid of bad data before blink
      % Add timepoint to removed
      removed = [removed s];
      go=20; % Reset go
    else
      raw(s,2:4)=nan; % Get rid of a few extra timepoints before the blink
                      % to make sure the data is clean
      removed = [removed s];
      go=go-1;
    end
    s=s-1; % Look at previous timepoint
  end
end


% Delete data after end of spike
for i = 1:length(spikes)
  e = spikes(i) + 1; % Timepoint right after the end of the blink
  go = 20; % Will continue deleting until there are this many consecutive 
           % non-blink datapoints
  while e<length(x) & go>0
    % Check to see if y value shows a significant jump in the data
    if (abs(ydiff(e)) > (ydiffmean + 3*ydiffstd))
      raw(e,2:4)=nan; % Get rid of bad data before blink
      removed = [removed s];
      go=20; % Reset go 
    else
      raw(e,2:4)=nan; % Get rid of a few extra timepoints before the blink
                      % to make sure the data is clean
      removed = [removed s];
      go=go-1;
    end
    e=e+1;
  end
end

% Return processed
processed = raw;

end