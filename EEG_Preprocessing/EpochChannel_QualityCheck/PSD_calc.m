%% CALCULATE THE POWER SPECTRUM
% Function to calculate the PSD of input eeg data channel using a pwelch.
% It uses a smoothing function smooth out small changes.
% Allows user to click on electrodes to highlight electrode and display its
% label, thus facilitating the detection of bad electrodes.
% DataIn : input data matrix
% fs : sampling frequency (Hz)
% chan: cell array of channel names
% chan_ind : vector of channel indices to include in analysis
% freqs : vector of the range of frequencies of interest
% datanom: title of the current dataset (EEG.setname)
% scale: log ou lin
% Exemple :
%PSD_calc(EEG.data, EEG.srate, {EEG.chanlocs.labels}, 1:64, [1:0.1:60],EEG.setname);
function [PSDOut,f1]=PSD_calc(DataIn, fs, chans,chan_ind, freqs,datanom,scale)
% Date : 01 2015   Programmed by: D. Bolger
%******************************************************
f1=figure;
wbh=waitbar (0,'Please wait...');  %Initialise the waitbar
h=zeros(length(chan_ind),1);       % Initialise the handles vector
if freqs(end)>(fs/2)
    error('Highest frequency of interest is higher than Nyquist!!')
    return;
end

%% If the data is segmented ...
if size(DataIn,3)>1
    seg=size(DataIn,3);  %number of trials
    pcent=20;   %take into account 20% of the data
    pcent_seg=ceil(seg*(pcent/100));
    segind=randi(seg,[pcent_seg,1]);
    DataIn=DataIn(chan_ind,:,segind);
end

%% Carry out PSD analysis
data_seg=0.5;
seg=ceil(size(DataIn,2)*data_seg);
DataIn_seg=DataIn(chan_ind,1:seg,size(DataIn,3));   %Just taking a segment of the data to calculate the PSD
assignin('base','DataIn_seg',DataIn_seg); 
Nfft=2.^nextpow2(size(DataIn_seg,2));                %Nfft length
assignin('base','Nfft',Nfft); 
PSDOut=cell(1,length(chan_ind));

for counter =1: length(chan_ind) %tsake each electrode in turn
   
    for scounter=1:size(DataIn,3)
           
        na=20;  %the averaging factor
        nx=length(DataIn_seg(chan_ind(counter),:));  %the window function
        fenetre=blackmanharris(floor(nx/na));  %define the windowing function
        olap=length(fenetre)/2;     %define the overlap function
        DIn=DataIn_seg(chan_ind(counter),:);
        [Pxx, Freqx]=pwelch(DIn,fenetre,olap,Nfft,fs);   % calculates psd using pwelch method
        assignin('base','Pxx',Pxx); 
        assignin('base','Freqx',Freqx); 
        fbin=Freqx(2)-Freqx(1);    %the frequency bin
        CG = sum(fenetre)/(nx/na);  % calculate the coherent gain
        NG = sum(fenetre.^2)/(nx/na);  %calculate the noise gain
        ifreqs=round(freqs./fbin)+1;
        
        %% SMOOTHING OF THE RESULTING ENVELOPE USING A MOVING AVERAGER
        intval=(1/fs)*2;  %specifying the interval over which to carry out the smoothing
        span=intval/(1/fs); %number of points over which to apply MA
        Pxx_smooth=smooth(Pxx,span,'moving'); %call of smoothing function
        assignin('base','Pxx_smooth',Pxx_smooth); 
        
        
    end %end of scounter loop
    
    Pxx_smooth2=mean(Pxx_smooth,2);
    PSDOut{1,counter}=Pxx_smooth2(ifreqs);
    
    %% PLOT THE PSD FOR EACH ELECTRODE SPECIFIED
    if strcmp(scale,'lin')==1
    
    h(counter)=plot(Freqx(ifreqs),Pxx_smooth2(ifreqs));   %
    set(h(counter),'HitTest','on','SelectionHighlight','on','UserData',chans{counter});
    set(h(counter),'ButtonDownFcn',@dispElectrode);   % displays the channel label upon mouse click
    
    hold all
    xlabel('Frequency (Hz)')
    ylabel('Power Spectral Density (V^2/Hz)')
    title(strcat('PSD: ',datanom))
    xlim([freqs(1) freqs(end)]);
    
    waitbar(counter/length(chan_ind));
    
    elseif strcmp(scale,'log')==1
        
        h(counter)=plot(Freqx(ifreqs),10*log10(Pxx_smooth2(ifreqs)));   %
        set(h(counter),'HitTest','on','SelectionHighlight','on','UserData',chans{counter});
        set(h(counter),'ButtonDownFcn',@dispElectrode);   % displays the channel label upon mouse click
        
        hold all
        xlabel('Frequency (Hz)')
        ylabel('Log Power Spectral Density (V^2/Hz)')
        title(strcat('PSD: ',datanom))
        xlim([freqs(1) freqs(end)]);
        
        waitbar(counter/length(chan_ind));
        
    end
    
end  %end of for loop



delete(wbh)  %close the waitbar




end

%% CALL OF BUTTONDOWNFCN FUNCTION
function dispElectrode(hdl,~)

disp(get(hdl,'UserData'));
set(hdl,'LineWidth',2.5);


end

%end of function