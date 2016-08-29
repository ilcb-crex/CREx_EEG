function [NoisyDetect, F1]=NoisyChanDetector(DataIn, NoisyDetect)
% Programmed by: D. Bolger                   June 2015
% Method to isolate noisy channels (adapted from the Prep Pipeline (Bigdely-Shamlo N, Mullen T, Kothe C, Su K-M & Robbins KA; 2015))
% The robust z-score is calculated using the median and the absolute deviant from the median; it is therefore less sensitive to outliers than the z-score.  
% Calculates:
%               - the robust standard deviation (z-score)
%               - the correlation between electrodes for low frequency band
% Define:
% NoisyConfig.ChanNum= [];
% NoisyConfig.robustChanDevThresh=5;
% Output:
% A structure, NoisyDetect, containing signal statistics.
%
%*********************************************************************************
chanlim= NoisyDetect.chanNum
ep=size(DataIn.data,3);
chans_inds=1:chanlim;
NoisyDetect.highFrequencyNoiseThreshold=0.4;   %z-score threshold

chan_noms={DataIn.chanlocs(1:chanlim).labels};
%% Calculate robust standard deviation
display('Calculating robust standard deviation...');

ChanDev=zeros(chanlim,ep);
robustChanDev=zeros(chanlim,ep);

fs      = DataIn.srate;
winpnts = floor(DataIn.times*fs/1000); % to samples
stepnts = floor((DataIn.times(2)-DataIn.times(1))*fs/1000);% to samples
dursam1 = DataIn.pnts;

for counter=1:ep
    
    ChanDev(:,counter)=0.7413*iqr(DataIn.data(1:chanlim,:,counter),2);
    
    ChanDevSTD=0.17413*iqr(ChanDev(:,counter));
    ChanDevMed=nanmedian(ChanDev(:,counter));
    robustChanDev(:,counter)=abs((ChanDev(:,counter)-ChanDevMed)./ChanDevSTD);
    
end

NoisyDetect.robustChanDev=robustChanDev;

% Identify the channels with unusally high SD according to this calculation
badChanDev=abs(NoisyDetect.robustChanDev)>NoisyDetect.robustChanDevThresh | isnan(NoisyDetect.robustChanDev);
NoisyDetect.NoisyChanDev=badChanDev;
NoisyDetect.ChanDevSTD=ChanDevSTD;
NoisyDetect.ChanDevMed=ChanDevMed;

%% Calculate correlation (best for continuous data)
% Determine level of low and high frequency noise for each electrode
% Filter design as it is based on the correlation of low frequency portion
% of signals
display('************************Calculating correlation of low frequency portion of EEG*******************************');

noisiness=zeros(chanlim, ep);
noisinessMedian=zeros(1,ep);
noisinessSD=zeros(1,ep);
badChannelsFromHFNoise=zeros(ep,chanlim);
zscoreHFNoise=zeros(ep,chanlim);

F=[2*[0 30 35]/DataIn.srate 1];   % to remove freqencies above 35 and below 1Hz
A= [1 1 0 0];   %amplitude to speficy low pass
N=100;   %filter order
nfft = max(DataIn.srate,2^ceil(log(N)/log(2)));   %no of fft bins
W = 0.54 - 0.46*cos(2*pi*(0:N)/N);                    %hamming window

% calculate interpolated frequency response
F = interp1(round(F*nfft),A,(0:nfft),'pchip');

% set phase & transform into time domain
F = F .* exp(-(0.5*N)*sqrt(-1)*pi*(0:nfft)./nfft);
B = real(ifft([F conj(F(end-1:-1:2))]));

% apply hamming window to kernel
B = B(1:N+1).*W(:)';
find(isnan(B))
X = zeros(size(DataIn.data,2),chanlim);

for k1=1:ep  %for each epoch
    d=DataIn.data(:,:,k1);
    d=d';
    size(d)
    for k = 1:chanlim
        
        X(:,k)= filtfilt(double(B), 1, double(d(:,k)));   %
    end

    noisi = mad(d(:,1:chanlim)- X, 1)./mad(X, 1);   % median absolute deviation
    noisi=squeeze(noisi);
    noisiMedian = nanmedian(noisi);
    noisiSD = mad(noisi, 1)*1.4826;
    
    noisiMedian=repmat(noisiMedian,1,size(noisi,1));
    noisiSD=repmat(noisiSD,1,size(noisi,1));
    zscoreHFNoiseTemp =abs (noisi - noisiMedian) ./ noisiSD;

    % Determine z-scored level of EM noise-to-signal ratio for each channel
    noiseMask = (zscoreHFNoiseTemp < NoisyDetect.highFrequencyNoiseThreshold) | ...
        isnan(zscoreHFNoiseTemp);
    
    badChannelsFromHFNoise(k1,:)  = noiseMask;
    
    noisiness(:,k1)=noisi;
    noisinessMedian(k1)=noisiMedian;
    noisinessSD(k1)=noisiSD;
    zscoreHFNoise(k1,:)=zscoreHFNoiseTemp;
end

badChan_T=cell(1,chanlim);

if ep>1
    
    
for chan_cnt=1:chanlim

 badChan_T{1,chan_cnt}=find(badChannelsFromHFNoise(:,chan_cnt));
 
end

elseif ep==1
    
   badChan_cont=find(badChannelsFromHFNoise);
   NoisyDetect.noisyChannels.badChan_cont=badChan_cont; 
   
end

NoisyDetect.noisyChannels.badChannelsFromHFNoise = badChannelsFromHFNoise;
NoisyDetect.noisyChannels.noiseMask = noiseMask;
NoisyDetect.noisyChannels.badChan_T=badChan_T;

NoisyDetect.zscoreHFNoise = zscoreHFNoise;
NoisyDetect.noisinessMedian = noisinessMedian;
NoisyDetect.noisinessSD = noisinessSD;

%% Calculate the Signal Statistics for each Electrode

if size(DataIn.data,3)==1
    
    display('****************************Calculating Signal Statistics**************************************************');

    M=cell(chanlim,ep);          %for mean
    SD=cell(chanlim,ep);         %for standard deviation
    med=cell(chanlim,ep);
    sk=cell(chanlim,ep);           %for skewness data
    k=cell(chanlim,ep);             %for kurtosis data
    zlow=cell(chanlim,ep);
    zhi=cell(chanlim,ep);
    tM=cell(chanlim,ep);          %trimmed mean
    tSD=cell(chanlim,ep);         %trimmed SD
    tndx=cell(chanlim,ep);        %index of the data retained after trimming
    ks_test=cell(chanlim,ep);     %Kolmogorov Smirnov test of gaussianity


  for ep_cnt=1:ep  
for chan_cnt=1:chanlim
    
    [M{chan_cnt,ep_cnt},SD{chan_cnt,ep_cnt},sk{chan_cnt,ep_cnt},k{chan_cnt,ep_cnt},med{chan_cnt,ep_cnt},zlow{chan_cnt,ep_cnt},zhi{chan_cnt,ep_cnt},tM{chan_cnt,ep_cnt},tSD{chan_cnt,ep_cnt},tndx{chan_cnt,ep_cnt},ks_test{chan_cnt,ep_cnt}] = signalstat(DataIn.data(chan_cnt,:), 0, [], [], strcat(DataIn.setname,'_channel statistics'),...
        mean(DataIn.data,2), DataIn.chanlocs);
    %ks_test(chan_cnt)=kstest(DataIn.data(chan_cnt,:));
end
  end
    
NoisyDetect.ChannelStats.Mean={{DataIn.chanlocs(1:chanlim).labels}',M};
NoisyDetect.ChannelStats.SD={{DataIn.chanlocs(1:chanlim).labels}',SD };
NoisyDetect.ChannelStats.Skew={{DataIn.chanlocs(1:chanlim).labels}',sk };
NoisyDetect.ChannelStats.Median={{DataIn.chanlocs(1:chanlim).labels}',med };
NoisyDetect.ChannelStats.Kurtosis={{DataIn.chanlocs(1:chanlim).labels}',k };
NoisyDetect.ChannelStats.zlow_qq={{DataIn.chanlocs(1:chanlim).labels}',zlow };
NoisyDetect.ChannelStats.zhi_qq={{DataIn.chanlocs(1:chanlim).labels}',zhi };
NoisyDetect.ChannelStats.tMean={{DataIn.chanlocs(1:chanlim).labels}',tM };
NoisyDetect.ChannelStats.tSD={{DataIn.chanlocs(1:chanlim).labels}',tSD };
NoisyDetect.ChannelStats.trim_indx={{DataIn.chanlocs(1:chanlim).labels}',tndx };
NoisyDetect.ChannelStats.KolSmirnov={{DataIn.chanlocs(1:chanlim).labels}',ks_test };
elseif size(DataIn.data,3)>1
    
    display('Epoched data! Could take too long to calculate all the statistics for all channels and epochs.');
    display('************************************Calculation of kurtosis**********************************************')
    
    [kurtosis,rej]=rejkurt(DataIn.data(1:chanlim,:,:),0,[],0);
    
    NoisyDetect.ChannelStats.Kurtosis={{DataIn.chanlocs(1:chanlim).labels}',kurtosis };
    NoisyDetect.ChannelStats.MeanKurtosis={{DataIn.chanlocs(1:chanlim).labels}',mean(kurtosis,2) };
    NoisyDetect.ChannelStats.RejKurtosis=rej;
    
end



%% Visualisation

F1=figure;

if ep==1
    
    ih=imagesc(NoisyDetect.robustChanDev);
    colormap(jet);
    colorbar()
    set(gca,'YTick',1:1:chanlim,'YTickLabel',{DataIn.chanlocs(1:chanlim).labels})
    set(get(gca,'XLabel'),'String','Trial Number','FontSize',12);
    set(get(gca,'YLabel'),'String','Electrode Labels','FontSize',12);
    set(get(gca,'Title'),'String','Noisy Electrodes: robust z-score: ','FontSize',12);
    set(gca,'FontSize',10);
    set(ih,'HitTest','on','SelectionHighlight','on','UserData',{chan_noms 1:chanlim},'XData',1:size(length(DataIn.times)),'YData',1:size(NoisyDetect.robustChanDev,1));
    
    dcm_obj=datacursormode(F1);
    set(dcm_obj,'UpdateFcn',@myupdateFcn)
    set(dcm_obj,'enable','on')
    
elseif ep>1
    
    ih=imagesc(NoisyDetect.robustChanDev);
    colormap(jet);
    colorbar()
    set(gca,'YTick',1:1:chanlim,'YTickLabel',{DataIn.chanlocs(1:chanlim).labels})
    set(get(gca,'XLabel'),'String','Trial Number','FontSize',12);
    set(get(gca,'YLabel'),'String','Electrode Labels','FontSize',12);
    set(get(gca,'Title'),'String','Noisy Electrodes: robust z-score: ','FontSize',12);
    set(gca,'FontSize',10);
    set(ih,'HitTest','on','SelectionHighlight','on','UserData',{chan_noms 1:chanlim},'XData',1:ep,'YData',1:size(NoisyDetect.robustChanDev,1));
    
    dcm_obj=datacursormode(F1);
    set(dcm_obj,'UpdateFcn',@myupdateFcn1)
    set(dcm_obj,'enable','on')
end

    function txt=myupdateFcn1(~, event_obj)
        % Sub-function to display the current time (seconds) and electrode on the
        % data-tips
        
        pos1=get(event_obj,'Position');
        txt={['Trial Number:  ',num2str(pos1(1))], ['Electrode: ',chan_noms{pos1(2)}]};
        
    end

    function txt=myupdateFcn(~, event_obj)
        % Sub-function to display the current time (seconds) and electrode on the
        % data-tips
        
        pos=get(event_obj,'Position');
        txt={['Robust z-score:  ',num2str(NoisyDetect.robustChanDev(pos(2)))], ['Electrode: ',chan_noms{pos(2)}]};
        
    end

end

