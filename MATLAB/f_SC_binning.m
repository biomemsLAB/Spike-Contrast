%% bin a sequence of events (e.g. time stamps of spikes)
% Input:    y: signal (e.g. time stamps)
%           time_start:         time where first bin starts (in seconds)
%           time_end:           time where last bin ends (in seconds)
%           bin: binsize in seconds to discretize spiketrain
%           step: bin-step, if step==bin -> no overlap
%           flag_binary=1: frequencies can be only 0 or 1, flag_binary=0: normal mode
% Output:   y_binned:       binned signal
%           x_step:         corresponding x values (=bin center) of y_binned
%           edges_step:     edges of bins
% needed functions: [y_binned,x_step,edges_step]=binning_halfOverlap(y,rec_dur,binsize,step,flag_binary) 
%
% written by Manuel Ciba, 2016/2017

function [y_binned,x_step,edges_step]=f_SC_binning(y,time_start,time_end,binsize,step,flag_binary) 

    if binsize/step == 2 % if half overlap use faster function
        [y_binned,x_step,edges_step]=f_SC_binning_halfOverlap(y,time_start,time_end,binsize,flag_binary);
    else
        %% Init
        edges_bin=time_start:binsize:time_end;      % edges of bins without step
        x_bin=edges_bin(1:end-1) + binsize/2;       % center of bins without step

        numSteps=binsize/step;                      % number of steps per bin
        numBins=length(x_bin);                      % number of bins without step
        numEdges=length(edges_bin);
        temp_binned=zeros(numBins,numSteps);        % each row contains binned signal for another step
        temp_x=zeros(numBins,numSteps);             % temp_x containts corresponding x values for temp_binned
        temp_edges=zeros(numEdges,numSteps);        % contains corresponding edges for temp_binned

        %% 1) calc one histograms for each overlapping step. Bin-start is shiftet by "step" respectively.       
        for i=1:numSteps
            %edges=0+((i-1)*step):binsize:rec_dur+((i-1)*step);
            edges_nu = edges_bin+((i-1)*step);
            temp_x(:,i)= edges_nu(1:end-1) + binsize/2;
            temp_edges(:,i)= edges_nu;
            temp_binned(:,i)=histcounts(y,edges_nu);
        end

        %% 2) merge all histograms to one histogram
        y_binned=zeros(numSteps*numBins,1);         % final histogram
        x_step=zeros(numSteps*numBins,1);           % corresponding x values of final histogram
        edges_step=zeros(numSteps*(numBins+1),1);   % corresponding edges of final histogram
        for k=1:numSteps
            y_binned(k:numSteps:end)=temp_binned(:,k);
            x_step(k:numSteps:end)=temp_x(:,k);
            edges_step(k:numSteps:end)=temp_edges(:,k);
        end

        %% 3) if flag_binary is true set all values >0 to 1
        if flag_binary
           y_binned(y_binned>0)=1; 
        end
    
    end

end