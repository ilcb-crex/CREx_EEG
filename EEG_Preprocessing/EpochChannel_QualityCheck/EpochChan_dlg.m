
function [NoisyOutput, DataIn]= EpochChan_dlg(DataIn)
% Programmed by : D. Bolger                Date: June 2015
% A function to enable the user to identify bad electrodes for continuous
% data and both bad electrodes and bad trials for epoched EEG data.
% It incorporates the idea that choice may have to be made, in epoched data in particular, between
% rejecting noisy electrodes that are "bad" for a certain number of epochs
% or rejecting epoches. 
% The function opens a dialogue box in which the user can specify:
% 1. A threshold (mV) for the identification of noisy data (default: 80mV).
% 2. A maximum number of above-threshold electrodes per epoch before
% considering the epoch as "bad" (default: 3).
% 3. The number of electrodes over which to carry out the analysis
% (default: length(EEG.chanlocs)).
% 4. A threshold z-score value for the identification of "bad" electrodes.
% 5 & 6: The upper and lower frequency bounds for the PSD calculation (only
% applicable for the moment on continuous data)(defaults: 30Hz & 1Hz, respectively).
% 6. If the user wishes to reject ica components (default: No). 
% Calls the EpochChan_Reject() function
% Exemple:
% EpochChan_dlg(EEG)
% It takes as input the EEGLAB dataset structure.
% 
%
%*********************************************************************

repeat='Yes';
while strcmp(repeat,'Yes')==1
    
    prompt={'Enter the threshold (mV): ', 'Enter max. no. above-threshold electrodes: ', 'Number of Electrodes: ','Enter z-score threshold value','Lower Frequency Limit (Hz) for PSD: ','Upper Frequency Limit(Hz) for PSD: ','Frequency Step: ','Remove ica components (Yes/No)'};
    dlg_title='Input Parameters';
    num_lignes=1;
    chan_num=length({DataIn.chanlocs.labels});
    deflts={'80','3', num2str(chan_num),'5','1','30','0.1','No'};
    thresh_ans=inputdlg(prompt,dlg_title,num_lignes,deflts);
    
    thresh_mV=str2double(thresh_ans{1,1});
    thresh_Elec=str2double(thresh_ans{2,1});
    numElec=str2double(thresh_ans{3,1});
    
    NoisyOutput=[]; %Initialise structure in which to record data concerning analysis of quality of data
    NoisyOutput.chanNum=numElec;         %number of electrodes
    NoisyOutput.thresh_mV=thresh_mV;     %the threshold(mV) for identifying bad data
    NoisyOutput.thresh_Elec=thresh_Elec;   %define the number of above threshold electrodes per trial
    NoisyOutput.robustChanDevThresh=str2double(thresh_ans{4,1});
    NoisyOutput.icaans=thresh_ans{8,1};
    
    freqs=[str2double(thresh_ans{5,1}):str2double(thresh_ans{7,1}):str2double(thresh_ans{6,1})];  %specify the frequency range over which to calculate the PSD
    repeat=[];
    
    %% REMOVE (IF SPECIFIED) ICA COMPONENTS
    
    if strcmp(NoisyOutput.icaans,'Yes')==1
        
        if isempty(DataIn.reject.icarejmanual)==0
            warn_comps(DataIn.reject.icarejemanual);
        end
        
        prompt2={'Which components would you like to remove? '};
        dlg_title2='Components to remove';
        num_lignes2=[1 20];
        deflts2={'1'};
        comp_ans=inputdlg(prompt2,dlg_title2,num_lignes2,deflts2);
        
        comps_remov=str2double(comp_ans{1,:});
        
        DataIn=pop_subcomp(DataIn,comps_remov,1);  %call of function pop_subcomp() to remove ICA components. 
        DataIn.reject.icarejmanual=comps_remov;          %record the indices of the removed components in the EEG structure. 
    else
        
        display('Not removing any ica components for the moment!');
        
    end
    %% 
    
    [BadT,BadChanAmp]=EpochChan_Reject( DataIn, numElec, thresh_mV, thresh_Elec);  
    DataIn.BadChanAmp=BadChanAmp;     % record the indices of the bad channels identified on the basis of their amplitudes
    
    NoisyOutput.badChanAmp=BadChanAmp; %number of time points in which each electrode > threshold
    NoisyOutput.badTAmp=BadT;                   %number of trials with an above-threshold number of above-threshold electrodes
    
    [NoisyOutput, F1]=NoisyChanDetector(DataIn, NoisyOutput);
    
    if size(DataIn.data,3)==1
        display('Continuous data. Calculate the PSD per Electrode.');
        f1=PSD_calc(DataIn.data, DataIn.srate, {DataIn.chanlocs.labels},1:chan_num, freqs,DataIn.setname);
        uiwait(f1)
    elseif size(DataIn.data,3)>1
        display('No PSD calculation of segmented data.');
        uiwait(F1)
    end
    
    
    repeat=repdlg;
    
end

end

function rep_ans=repdlg

prompt={'Do you wish to change parameters and analyse again? (Yes/No): '};
dlg_title='Repeat Analysis Query';
num_lignes=1;
deflts={'No'};
ansr=inputdlg(prompt,dlg_title,num_lignes,deflts);
rep_ans=ansr{1,1};

end

function warn_comps(ica_reject)

d=dialog('Position',[300 300 250 150],'Name','ICA components removed!');

txt=uicontrol('Parent',d,'Style','text','Position',[20 80 210 40],'String',strcat('The following components were already removed: ',double2str(ica_reject)));
btn=uicontrol('Parent',d,'Position',[85 20 70 25],'String','Close','Callback','delete(gcf)');

end
