function [BadTrials badtrial_total]=EpochChan_Reject( DataStruct, nE, thresh, badE)
%*******************************************************************
%Date: 3-06-2015      Programmed by: Deirdre Bolger
%Function to isolate electrodes and channels for rejection
%Input:
%        DataStruct ==> the EEG structure
%        nE ==> actual number of electrodes
%        thresh ==> threshold for detection of artifactual activity
%        badE ==> number of bad electrodes per trial before considering a
%        trial as bad.
%
%The function can take electrode X TimePoint data or
%                                electrode X TimePoint X Trial data
%It includes the possibility of specifying the time interval over which to
%analyse the data. 
%
%This function is called by the function EpochChan_dlg() which allows the
%user to specify the number of electrodes (nE) and the threshold data
%(thresh and badE). It also allows the user to redo the analysis,
%redefining the parameters. 
%
%*******************************************************************
%% DEFINE CHANNEL NUMBERS AND

chans=1:nE;
DataIn=DataStruct.data;
DataIn=DataIn(1:nE,:,:);
[~,~,p]=size(DataIn);
fs=DataStruct.srate;
fpath=DataStruct.filepath;
fnom=DataStruct.setname;
global chan_noms
chan_noms={DataStruct.chanlocs.labels};
global myData
myData=[];
global Artcnt_sec;


%% CHECK TO SEE IF DATA IS CONTINUOUS OR EPOCHED

if p==1
    disp('Continuous data')
    
    prmpt_t={'Specify a different start time interval: ','Specify a different end time interval'};
    dlg_title_t='Time Interval Query';
    num_lignes_t=1;
    deflts_t={num2str(DataStruct.times(1)), num2str(DataStruct.times(end))};
    ans_t=inputdlg(prmpt_t, dlg_title_t, num_lignes_t, deflts_t);
    rept_start=str2double(ans_t{1,1});
    rept_end=str2double(ans_t{2,1});
    
elseif p>1
    disp('Epoched data')
    
    prmpt_t={'Specify a different start time interval: ','Specify a different end time interval'};
    dlg_title_t='Time Interval Query';
    num_lignes_t=1;
    deflts_t={num2str(DataStruct.times(1)), num2str(DataStruct.times(end))};
    ans_t=inputdlg(prmpt_t, dlg_title_t, num_lignes_t, deflts_t);
    rept_start=str2double(ans_t{1,1});
    rept_end=str2double(ans_t{2,1});
    
end

%% FOR EACH TRIAL FIND NUMBER OF ELECTRODES EXCEEDING THRESHOLD VOLTAGE
t_int=find(DataStruct.times>=rept_start & DataStruct.times<=rept_end);
Artcnt=zeros(nE,p);  % initialise variable
Zerocnt=zeros(nE,p);

for pcount=1:p   %for each trial
    
    
    for ecount=1:nE  % for each electrode
        
        i1=find(abs(DataIn(chans(ecount),t_int,pcount))>=thresh);     %indices of timepoints for current electrode where voltage > 100 or < 100
        iz=find(abs(DataIn(chans(ecount),t_int,pcount))==0);
        
        if isempty(i1)==1
            nArt=0;
            disp('No data points exceeding threshold');
            
        else
            nArt=length(i1);
           
        end
         
        Artcnt(chans(ecount),pcount)=nArt;   % No. of data-points exceeding threshold for the current electrode/current trial
        
        if isempty(iz)==1
            nzero=0;
            %disp('No data points equalling zero')
            
        else
            nzero=length(iz);
        end
        
        Zerocnt(chans(ecount),pcount)=nzero;  %No. of data-points equalling zero for the current electrode/trial
        assignin('base','ZeroCount',Zerocnt);
        
    end   %end of electrode count
    
end   %end of pcount (trial count)

Artcnt_total=sum(Artcnt,1);  %For each trial find the number of total above-threshold points
iArtcnt_total=(Artcnt_total>0);

%% FOR EACH ELECTRODE FIND NUMBER OF TRIALS EXCEEDING THRESHOLD VOLTAGE
ecount=[];  %initialise variables
numE=zeros(nE,1);
numT=zeros(1,p);
BadTrials = [];

for ecount=1:nE
    
    numE(chans(ecount))=length(find(Artcnt(chans(ecount),:)>0));           %find number of trials in which current electrode exceeds threshold
    
end

if isequal(unique(Artcnt),0)==0  % if electrodes with above-threshold noise have been isolated
    
    
    for pcount=1:p
        
        numT(pcount)=length(find(Artcnt(:,pcount)>0)); %find no. bad electrodes in each trial
        
    end
    
    BadTrials = [find(numT>=badE)' numT(numT>=badE)'];
    x=find(numT>=badE); % find those trials in which more than x electrodes are bad
    badtrial_total=unique(cat(2,x,find(iArtcnt_total)));
    
    winrej_mat=zeros(length(badtrial_total), 5+nE);
    
    if x>=1   
        
    for w_cnt=1:length(badtrial_total)
        
    u=(size(DataIn,2)*badtrial_total(w_cnt))-size(DataIn,2);
    winrej_mat(w_cnt,:)=[u u+size(DataIn,2)  0 1 0 zeros(1,nE)] ;
 
    
    end
    
    eegplot(DataIn(:,t_int,:), 'srate',fs,'eloc_file',DataStruct.chanlocs(1:nE),'dispchans',nE,'title','Epochs Marked for Rejection','xgrid','off','ygrid','off','ploteventdur','on','command','reject','winrej',winrej_mat);  %,
    else
       display(['No trial with more than ',char(num2str(badE)),' above-threshold electrodes.'])
    end
    
    
elseif isequal(unique(Artcnt),0)==1
    
    disp(strcat('No data points exceed defined threshold of ',num2str(thresh),'mV'));
    
    
end  %end of if

f1=figure;

if p>1
    
    ih=imagesc(Artcnt);
    colormap(jet);
    colorbar()
    set(gca,'YTick',1:1:nE,'YTickLabel',{DataStruct.chanlocs.labels})
    set(get(gca,'XLabel'),'String','Trial Number','FontSize',12);
    set(get(gca,'YLabel'),'String','Electrode Labels','FontSize',12);
    set(get(gca,'Title'),'String',strcat('Bad Electrodes and Trials Analysis Image: ',num2str(rept_start),'ms - ', num2str(rept_end), 'ms'),'FontSize',12);
    set(gca,'FontSize',10);
    set(ih,'HitTest','on','SelectionHighlight','on','UserData',{chan_noms 1:nE},'XData',1:size(Artcnt,2),'YData',1:size(Artcnt,1));
    
    dcm_obj=datacursormode(f1);
    set(dcm_obj,'UpdateFcn',@myupdateFcn)
    set(dcm_obj,'enable','on')
    
    f2=figure
    Artcnt_total_sec=Artcnt_total.*(1/fs);
    Trial_num=[1:size(DataIn,3)]
    bar(Artcnt_total_sec)
    set(gca,'XLim',[1 size(DataIn,3)],'XTick',1:1:size(DataIn,3),'XTickLabel',1:size(DataIn,3),'FontSize',8)
    set(gca,'YLim',[0 max(Artcnt_total_sec)],'YTick',1:5:max(Artcnt_total_sec),'FontSize',8)
    set(get(gca,'XLabel'),'String','Trial Number','FontSize',12)
    set(get(gca,'YLabel'),'String','Total above-threshold time (s)','FontSize',12)
    set(get(gca,'Title'),'String',strcat('Total Above-Threshold Time (seconds) for each trial' ),'FontSize',12);
    set(gca,'HitTest','on','SelectionHighlight','on','UserData',{Trial_num}); %#ok<NBRAK>
    
    dcm_obj3=datacursormode(f2);
    set(dcm_obj3,'UpdateFcn',@myupdateFcn3)
    set(dcm_obj3,'enable','on')
    
    Rdata=[{DataStruct.chanlocs(1:nE).labels}' num2cell(numE)];
    
elseif p==1  %if the data is continuous
    
    Artcnt_sec= Artcnt.*(1/fs);
    bar(1:nE, Artcnt_sec);
    set(gca,'XLim',[1 nE],'XTick',1:1:nE,'XTickLabel',{DataStruct.chanlocs(1:nE).labels},'FontSize',8)
    set(gca,'YLim',[0 max(Artcnt_sec)],'YTick',1:5:max(Artcnt_sec),'FontSize',8)
    set(get(gca,'XLabel'),'String','Electrode Labels','FontSize',12)
    set(get(gca,'YLabel'),'String','Outside Threshold Time Intervals (s)','FontSize',12)
    set(get(gca,'Title'),'String',strcat('Total Above-Threshold Time (seconds) for each Electrode: ',num2str(rept_start),'ms - ', num2str(rept_end), 'ms' ),'FontSize',12);
    set(gca,'HitTest','on','SelectionHighlight','on','UserData',{chan_noms Artcnt_sec});
    
    dcm_obj2=datacursormode(f1);
    set(dcm_obj2,'UpdateFcn',@myupdateFcn2)
    set(dcm_obj2,'enable','on')
    
    Rdata=[{DataStruct.chanlocs(1:nE).labels}' num2cell(numE)];
    
%% Plot eegplot of continuous data with possible bad electrodes marked. 

   
     ArtcntMean=repmat(mean(Artcnt_sec,1),[nE,1]);
     Adiff=Artcnt_sec-ArtcntMean;
     ArtcntStd=std(Artcnt_sec);
     ArtcntStd_mat=repmat(ArtcntStd',[nE,1]);
     Artcnt_zscor=abs(Adiff./ArtcntStd_mat)
     
     ibadE=(Artcnt_zscor>4);
     
     Erej_mat=[1 size(DataIn,2)  1 1 1 ibadE'];
     eegplot(DataIn, 'srate',fs,'eloc_file',DataStruct.chanlocs(1:nE),'dispchans',nE,'title','Potentially Bad Electrodes Marked','xgrid','off','ygrid','off','ploteventdur','on','command','reject','winrej',Erej_mat);  %,
     
end %% p if

    function txt=myupdateFcn(~, event_obj)
    % Sub-function to display the trial number and electrode on the
    % data-tips
        
        pos=get(event_obj,'Position');
        txt={['Trial: ',num2str(pos(1))], ['Electrode: ',chan_noms{pos(2)}]};
        
    end

    function txt=myupdateFcn2(~, event_obj)
    % Sub-function to display the current time (seconds) and electrode on the
    % data-tips
        
        pos=get(event_obj,'Position');
        txt={['Time Int (s):  ',num2str(Artcnt_sec(pos(1)))], ['Electrode: ',chan_noms{pos(1)}]};
        
    end

    function txt=myupdateFcn3(~, event_obj)
    % Sub-function to display the current time (seconds) and electrode on the
    % data-tips
        
        pos=get(event_obj,'Position');
        txt={'Trial Number: ', num2str(Trial_num(pos(1)))};
        
    end



   
end

