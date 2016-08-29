%Date: 21-07-2016                        Programmed by: D. Bolger
%function to adjust the T0 of EEGLAB datasets.
%***********************************************************************

prompt={'Enter Condition name to process: ', 'Group to Process (NOL/DYS): ', 'Enter new T0 (ms after current T0):  '};
dlg_title='Topo Plot Parameters';
num_lignes=[1, 1, 1];
deflts={'MotsCC','NOL','100'};
answr=inputdlg(prompt,dlg_title,num_lignes,deflts);

cond=cellstr(answr{1,1});
group=cellstr(answr{2,1});
newT0=str2double(cellstr(answr{3,1}));
currdir='F:/BLRI/EEG/Projets_en_cours/Projet_MotInter/ExpEEG_Phase1/Data_Biosemi/P1AUD_Results/';                %basic file path-need to change this. 
Condnom={cond,group};

%% LOAD SET FILES INTO EELAB SESSION

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; 
currdirectory=strcat(currdir,cond,'/',group,'/');

files_cond= dir(currdirectory{1,1});
currfiles=files_cond;
fileIndex= find(~[currfiles.isdir]);
FInd=fileIndex;
filenum=dir(strcat(currdirectory{1,1},'*.set'));
filenom={filenum.name};                           %titles of the all the .set files to be loaded
EEG = pop_loadset('filename',filenom,'filepath',currdirectory{1,1});
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
EEG=eeg_checkset(EEG); eeglab redraw;

%%

for Count=1:length(ALLEEG)
    
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',Count,'study',0);
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
    Time=EEG.times;                                                                     %define the new T0 in ms
    diffs=(1/EEG.srate)*1000;                                                   %time interval in ms
    i0=find(Time==0);                                                             %find the index of current T0
    inew=find(Time>=(newT0-diffs) & Time<=(newT0+diffs));
    [~,i]=min(abs(newT0-Time(inew))); indxnew=inew(i);             %index of the new T0
    
    bline_new=Time(1)-Time(indxnew):diffs:(diffs*-1);
    poststim_new=0:diffs:(Time(end)-Time(indxnew));
    TimeNew=cat(2,bline_new,poststim_new);
    
    EEG.times=[]; ALLEEG(Count).times=TimeNew;
    EEG.xmin=TimeNew(1)/1000;
    EEG.xmax=ALLEEG(Count).times(end)/1000;
    lat=[EEG.event.latency];
    latnew=lat+(indxnew-i0);
    
    for counter=1:length(ALLEEG(Count).event)
        
        EEG.event(counter).latency=latnew(counter);
        
    end
    
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    
    title_new=strcat(EEG.setname,'-T0new');                                                  %rename the new file with new T0
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(title_new),'gui','off'); % current set = xx;
    EEG = eeg_checkset( EEG );
    savepath=strcat(currdir,cond,'/',group,'-T0/');
    EEG = pop_saveset( EEG, 'filename',char(title_new),'filepath',savepath{1,1});  % Saves a copy of the current resampled dataset to the current directory
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw
end




