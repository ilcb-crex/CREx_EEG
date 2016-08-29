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
% Requires the smooth.m function from the Curve Fitting Toolbox.
% Exemple :
%PSD_calc(EEG.data, EEG.srate, {EEG.chanlocs.labels}, 1:64, [1:0.1:60],EEG.setname);
function [PSDOut,f1]=PSD_calc(DataIn, fs, chans,chan_ind, freqs,datanom)
% Date : 01 2015   Programmed by: D. Bolger
%******************************************************
f1=figure;
wbh=waitbar (0,'Please wait...');  %Initialise the waitbar
h=zeros(length(chan_ind),1);       % Initialise the handles vector
if freqs(end)>(fs/2)
    error('Highest frequency of interest is higher than Nyquist!!')
    return;
end

PSDOut=cell(1,length(chan_ind));

    for counter =1: length(chan_ind) %take each electrode in turn
      
            
        Nfft=2.^nextpow2(size(DataIn,2));  %Nfft length
      
        na=20;  %the averaging factor
        nx=length(DataIn(chan_ind(counter),:));  %the window function
        fenetre=blackmanharris(floor(nx/na));  %define the windowing function
        olap=length(fenetre)/2;     %define the overlap function
        DIn=DataIn(chan_ind(counter),:);
        [Pxx, Freqx]=pwelch(DIn,fenetre,olap,Nfft,fs);   % calculates psd using pwelch method
        fbin=Freqx(2)-Freqx(1);    %the frequency bin
        CG = sum(fenetre)/(nx/na);  % calculate the coherent gain
        NG = sum(fenetre.^2)/(nx/na);  %calculate the noise gain
        ifreqs=round(freqs./fbin)+1;
        
        %% SMOOTHING OF THE RESULTING ENVELOPE USING A MOVING AVERAGER
        intval=0.1;  %specifying the interval over which to carry out the smoothing
        span=intval/(1/fs); %number of points over which to apply MA
        Pxx_smooth=smooth(Pxx,span,'moving'); %call of smoothing function
        PSDOut{1,counter}=Pxx_smooth(ifreqs);
        
        %% PLOT THE PSD FOR EACH ELECTRODE SPECIFIED
        h(counter)=plot(Freqx(ifreqs),10*log10(Pxx_smooth(ifreqs)));   %
        set(h(counter),'HitTest','on','SelectionHighlight','on','UserData',chans{counter});
        set(h(counter),'ButtonDownFcn',@dispElectrode);   % displays the channel label upon mouse click
        
        hold all
        xlabel('Frequency (Hz)')
        ylabel('Power Spectral Density (V^2/Hz)')
        title(strcat('PSD: ',datanom))
        xlim([freqs(1) freqs(end)]);
        
        waitbar(counter/length(chan_ind));
        
    end  %end of for loop



delete(wbh)  %close the waitbar




end

%% CALL OF BUTTONDOWNFCN FUNCTION
function dispElectrode(hdl,~)

disp(get(hdl,'UserData'));
set(hdl,'LineWidth',2.5);


end

%end of function
