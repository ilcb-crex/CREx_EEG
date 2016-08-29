
function [AllData, gfp_cond,pvalue_corr,Cond_data]=CREx_GroupPlot_v2(cfg)
% Date: 7-12-2015             Programmed by : Deirdre Bolger
% For this script to run you have placed the datasets for each condition in a
% separate folder with the condition name.
% It loads all the subject-level datasets for each condition.
% For each condition it calculates the Grand Average over all subjects -
% for GFP calculation.
% For all subjects, for each electrode of interest (cfg.eoi) it carries out a
% permutation test and plots those time points that exceed critical p-value
% (0.05) against the ERP with 95% confidence intervals/
% The Global Field Power (GFP) of each condition is plotted against the erp plot.
% Very Important! It works for two conditions only!!
% This script also accepts datasets whose T0 do not fall on the same
% time-point. Upon visualisation, the T0 points of each are forced to
% correspond and th T0 point is marked. The GFPs and GMDs are also
% adjusted.
%To adjust the T0 of epoched data apply the "CREx_AdjustT0.m"
% script.
% This script calls the functions :
%           - plot_ConfInt2( ) : to plot the 95% confidence intervals.
%           - plot_perm( )     : to carry out permutation test and plot results if any.
% Input to the functions:
% The function accepts as input a configuration structure (cfg) in which
% the following parameters are defined :
%       1. currdir - current directory in which the condition folders are
% located.
%       2. Condnom - the names of the conditions to be compared.
%       3. Subnum - the number of subjects;
%       4. Channum - the number of channels (electrodes) to include in
%       analysis.
%       5. TimeWin - the upper and lower limits (in ms) of the time window
%       of interest.
%       6. eoi - the electrodes of interest (labels)
%***************************************************************


%%  INITIALISE THE VARIABLES

length(cfg.Condnom)
dir_cond=cell(1,length(cfg.Condnom));
AllData=cell(cfg.Subnum,length(cfg.Condnom));
files_cond=cell(1,length(cfg.Condnom));
fileIndex=cell(1,length(cfg.Condnom));
CondData=cell(1,length(cfg.Condnom));
CondGA=cell(1,length(cfg.Condnom));
Dirs=cfg.dirs;


%% FOR EACH CONDITION LOAD THE SUBJECT LEVEL FILES



for condcnt=1:length(cfg.Condnom)
    
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;  %open an EEGLAB session
    
    
    dir_cond{1,condcnt}=Dirs{1,condcnt};    %the directory of the current condition
    
    display(dir_cond{1,condcnt});                    %display the contents of the current folder
    
    files_cond{1,condcnt}= dir(dir_cond{1,condcnt});
    currfiles=files_cond{1,condcnt};
    fileIndex{1,condcnt} = find(~[currfiles.isdir]);
    FInd=fileIndex{1,condcnt};
    filenum=dir(strcat(dir_cond{1,condcnt},'*.set'));
    filenom={filenum.name};                           %titles of the all the .set files to be loaded
    EEG = pop_loadset('filename',filenom,'filepath',dir_cond{1,condcnt});
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
    EEG=eeg_checkset(EEG); eeglab redraw;
    
    for counter=1:length(ALLEEG)
        AllData{counter,condcnt}=ALLEEG(counter);
    end
    
end

size(AllData)
assignin('base','AllData',AllData);

%% DEFINE THE TIME INTERVAL AND ADJUST TIME VECTOR OF T0 DIFFER BETWEEN CONDITIONS
adder_start=[];
adder_end=[];
if strcmp(cfg.TimeWind,'all')==1
    
    if isequal(AllData{1,1}.xmin,AllData{1,2}.xmin)==1   %if datasets of both conditions use the same T0
        Time=AllData{1,1}.times;
        Tind=1:AllData{1,1}.pnts;
    elseif isequal(AllData{1,1}.xmin, AllData{1,2}.xmin)==0   %if datasets of both conditions apply different T0
        
        [T0ind,isort]=sort([find(AllData{1,1}.times==0) find(AllData{1,2}.times==0)],'descend');
        T0ind_diff=T0ind(1)-T0ind(2);                                                                            %this difference needs to be added to the end of dataset with earlier T0 (1) and to the start of dataset with later T0 (2).
        intval=(1/AllData{1,1}.srate)*1000;
        adder_start=AllData{1,isort(1)}.times(1):intval:(AllData{1,isort(2)}.times(1)-intval);    %time interval to be added at the start
        adder_end=(AllData{1,isort(1)}.times(end)+intval):intval:(AllData{1,isort(2)}.times(end));  %time interval to be added at the end
    end
else
    Time=AllData{1,1}.times;
    Tind=find(Time>=cfg.TimeWind(1) & Time<=cfg.TimeWind(2));
    Time=Time([Time>=cfg.TimeWind(1) & Time<=cfg.TimeWind(2)]);
end

%% EXTRACT DATA FOR EACH SUBJECT AND EACH CONDITION AND CALCULATE THE GRANDAVERAGE OVER ALL SUBJECTS (FOR GFP CALCULATION)

for condcnt2=1:size(AllData,2)
    currCond=zeros(cfg.Channum,size(AllData{1,1}.data,2),size(AllData,1));
    assignin('base','currCond',currCond);
    if isempty(adder_start)==1
        for sujcnt=1:size(AllData,1)
            currCond(:,:,sujcnt)=mean(AllData{sujcnt,condcnt2}.data(1:cfg.Channum,:,:),3);
        end
        
    elseif isempty(adder_start)==0 && (condcnt2==isort(1))==1   %if this is the condition with later T0 add adder_end to end
        disp('Later T0 dataset')
        Time=cat(2,AllData{1,condcnt2}.times,adder_end);
        Tind=1:length(Time)
        newpnts=zeros(cfg.Channum,length(adder_end),size(AllData,1));
        assignin('base','currCond',currCond)
        
        for sujcnt=1:size(AllData,1)
            display(sujcnt)
            currCond(:,:,sujcnt)=mean(AllData{sujcnt,condcnt2}.data(1:cfg.Channum,:,:),3);
        end
        currCond=cat(2,currCond,newpnts);
        
    elseif isempty(adder_start)==0 && (condcnt2==isort(2))==1   %if this is the condition with earlier T0 add adder_start to start
        disp('Earlier T0 dataset')
        
        Time=cat(2,adder_start,AllData{1,condcnt2}.times);
        Tind=1:length(Time);
        newpnts=zeros(cfg.Channum,length(adder_start),size(AllData,1));
        assignin('base','currCond',currCond)
        
        for sujcnt=1:size(AllData,1)
            display(sujcnt)
            currCond(:,:,sujcnt)=mean(AllData{sujcnt,condcnt2}.data(1:64,:,:),3);
        end
        
        currCond=cat(2,newpnts,currCond);
        
    else
        disp('LOST!!')
    end
    
    CondData{1,condcnt2}=currCond;
    CondGA{1,condcnt2}=mean(currCond,3);
    
end

assignin('base','CondGA',CondGA);
assignin('base','CondData',CondData);
assignin('base','Time',Time);


%% FOR EACH ELECTRODE OF INTEREST PLOT THE MEAN FOR EACH cfg.eoi


Cond_data=cell(length(cfg.eoi),length(cfg.Condnom));
Eint=zeros(length(cfg.eoi),1);

for ECnt=1: length(cfg.eoi)  %for each electrode of interest
    
    Eint(ECnt,1)=find(strcmp({AllData{1,1}(1).chanlocs.labels},cfg.eoi{1,ECnt})); % find the indice of the current electrode of interest
    
    for condcnt3=1:length(cfg.Condnom)
        
        Cond_data{ECnt,condcnt3}=squeeze(CondData{1,condcnt3}(Eint(ECnt,1),Tind,:))';
        
    end
    
end
%save('E:\BLRI\EEG\Projets_en_cours\projet_DECODAGE_LEXIQUE\cond_cell2.mat','Cond_data','CondData','CondGA','-v7.3')
%% PLOT THE GFP AND THE RESULTS OF THE PERMUTATION TEST

peplot=[1]   % 2 5 6 9 10 13 14];
gfpplot=[2]  % 4 7 8 11 12 15 16];
sh1=zeros(length(cfg.eoi),1);
sh2=zeros(length(cfg.eoi),1);
sh3=zeros(length(cfg.eoi),1);

if length(cfg.Condnom)==2   %if there are only two conditions
    
    colrs=eye(length(cfg.Condnom),3).*0.6;
    disp('*************Two experimental conditions. Will carry out permutation test.*******************');
    
    %% CALCULATE THE GLOBAL FIELD POWER (GFP)
    
    %     [gfp_cond1,tind1]=GFP_calc(CondGA{1,1},64,1:64,[Time(1) Time(end)],Time);   %Calculate the GFP for condition 1
    %     [gfp_cond2,tind2]=GFP_calc(CondGA{1,2},64,1:64,[Time(1) Time(end)],Time);   %Calculate the GFP for condition 2
    
    [SC,GMD,gfp_cond1,gfp_cond2]=GMD_calc(CondGA{1,1}, CondGA{1,2},64,1:64,[Time(1) Time(end)],[Time(1) Time(end)],Time); %Finds the Global Map Dissimilarity between the two conditions
    size(GMD)
    size(Time)
    
    
    %% PREPARE THE FIGURE TO PLOT THE ERPs WITH ELECTRODE CONFIGURATION
    
    hndl=figure; set(hndl,'PaperUnits','normalized','PaperPosition',[0.0235308 0.0272775 0.894169 0.909249],'Color',[1 1 1]);
    orient portrait; axis ('normal')
    
    xvals=zeros(length(EEG.chanlocs),1);
    yvals=zeros(length(EEG.chanlocs),1);
    pwidth    = 0.75;     % 0.75, width and height of plot array on figure
    pheight   = 0.75;
    axwidth=0.04;
    axheight=0.08;
    
    %Identify all the non-empty channel
    notmtchans=cellfun('isempty',{EEG.chanlocs.theta});
    notmtchans=find(~notmtchans);
    
    
    %Read in channel locations file
    [elocs,titres,theta,rads,inds]=readlocs(EEG.chanlocs(notmtchans));
    channoms=strvcat(EEG.chanlocs.labels);
    Th=pi/180*theta;           %convert degrees to radians
    
    %Convert from polar to cartesian
    [ycart,xcart]=pol2cart(Th,rads);
    xvals(notmtchans)=ycart;
    yvals(notmtchans)=xcart;
    
    %Find the positions of all the channels
    allchans=length(EEG.chanlocs);
    mtchans=setdiff(1:allchans,notmtchans);        %find the channels indices common to both
    allchans_sqrt=floor(sqrt(allchans))+1;
    
    for i=1:length(mtchans)
        
        xvals(mtchans(i))=0.7+0.2*floor((i-1)/allchans);
        yvals(mtchans(i))=-0.4+mod(i-1,allchans)/allchans;
        
    end
    
    channoms2=channoms(1:64,:);
    xvals=xvals(1:64);
    yvals=yvals(1:64);
    
    if length(xvals) > 1
        if length(unique(xvals)) > 1
            gcapos = get(gca,'Position'); axis off;
            xvals = (xvals-mean([max(xvals) min(xvals)]))/(max(xvals)-min(xvals)); % this recenters
            xvals = gcapos(1)+gcapos(3)/2+pwidth*xvals;                                 %this controls width of plot
        end;
    end;
    gcapos = get(gca,'Position'); axis off;
    yvals = gcapos(2)+gcapos(4)/2+pheight*yvals;  % controls height of plot
    ho=zeros(length(cfg.eoi),1);
    
    Axes = [];
    for ecounter=1:length(cfg.eoi)
        
        xcenter=xvals(ecounter);
        ycenter=yvals(ecounter);
        Axes=[Axes axes('Units','normalized','Position', [ycenter-axheight/2 xcenter-axwidth/2 axheight axwidth])];
        hold on;
        
        X1=Cond_data{ecounter,1}';
        X1mean=mean(X1,2);
        
        X2=Cond_data{ecounter,2}';
        X2mean=mean(X2,2);
        
        [axs1,ploth1]=plot_ConfInt2(X1,Time',X1mean,[colrs(1,:)], 0.25,0.05);
        set(axs1,'YDir','reverse');
        hold on
        [axs2,ploth2]=plot_ConfInt2(X2,Time',X2mean,[colrs(2,:)], 0.25,0.05);
        set(axs2,'YDir','reverse','XTickLabel',[],'XGrid','off','YGrid','off','Box','off')
        hold on
        ploth3=plot(Time,X2mean-X1mean,'r');
        pvalue_corr=plot_perm({X1 X2},Time',(cfg.ylims(2)-0.4),[Time(1) Time(end)]);
        title(cfg.eoi{1,ecounter});
        ylim(cfg.ylims);
        h1=vline(0,'k','0');
        
        
        set(hndl,'CurrentAxes',Axes(ecounter));
        set(Axes(ecounter),'HitTest','on','SelectionHighlight','on','UserData',{Time,X1,X2,X1mean,X2mean,colrs,cfg.eoi{1,ecounter},gfp_cond1,gfp_cond2,GMD,cfg.Condnom,cfg.ylims,CondGA,EEG.chanlocs},'NextPlot','new');
        set(Axes(ecounter),'ButtonDownFcn',@plotsingle_pe)
        
    end %%
    assignin('base','Axes',Axes);
    assignin('base','hndl',hndl);
    
    gfp_cond=cat(3,gfp_cond1,gfp_cond2);
    
elseif length(cfg.Condnom)>2
    
    disp('*************More than two experimental conditions. Cannot carry out permutation test.*******************');
    colrs1=eye(3,3).*0.5;
    colrs2=fliplr(eye(3,3).*0.2)+colrs1;
    Clrs=cat(1,colrs1,colrs2);
    
    pvalue_corr=[];
    gfp_cond=zeros(length(Time),length(cfg.Condnom));
    tind1=zeros(length(cfg.Condnom),length(Time));
    
    for ccounter=1:length(cfg.Condnom)
        
        [gfp_cond(:,ccounter),tind(ccounter,:)]=GFP_calc(CondGA{1,ccounter},64,1:64,[Time(1) Time(end)],Time);   %Calculate the GFP for condition 1(ccounter,:)(ccounter,1)
        
    end
    
    X1mean=zeros(length(Time),length(cfg.eoi),length(cfg.Condnom));
    axs1=zeros(length(cfg.Condnom),1);
    ploth1=zeros(length(cfg.Condnom),1);
    p1=zeros(length(cfg.Condnom),1);
    colors={'b' 'g' 'm' 'k' 'r' 'c'};
    
    
    
    for ecounter=1:length(cfg.eoi)    %for each electrode
        figure;
        set(gcf,'Color',[1 1 1])
        sh1(ecounter)=subplot(2,1,1);
        pos1=get(sh1(ecounter),'Position');
        pos1=[pos1(1)-0.05 pos1(2)-0.02 pos1(3)+0.05 pos1(4)+0.025];
        set(sh1(ecounter),'Position',pos1);
        pos1b=get(sh1(ecounter),'Position');
        
        for ccnter=1:length(cfg.Condnom)  %for each condition
            
            X1=Cond_data{ecounter,ccnter}';       %
            X1mean(:,ccnter,ecounter)=mean(X1,2);
            
            [axs1(ccnter,1),ploth1(ccnter,1)]=plot_ConfInt2(X1,Time',X1mean(:,ccnter,ecounter),[Clrs(ccnter,:)], 0.25,0.05);
            set(axs1(ccnter,1),'YDir','reverse');
            ylim([-8 8])
            hold on
        end
        h=vline(0,'k','T0');
        set(h,'LineWidth',2);
        
        ploth1
        title(horzcat('Global Field Power (GFP) of Conditions : Electrode :', char(cfg.eoi{1,ecounter})));
        L2=legend(ploth1,cfg.Condnom);
        set(L2,'Box','off');
        
        sh2(ecounter)=subplot(2,1,2,'Position',[pos1b(1) (pos1b(2)-pos1b(4)/1.5) pos1b(3) pos1b(4)/1.5]);
        for ccnter1=1:length(cfg.Condnom)  %for each condition
            
            
            p1(ccnter1)=plot(Time,gfp_cond(:,ccnter1),colors{1,ccnter1});  %plot the difference if necessary
            set(p1(ccnter1),'LineWidth',2);
            
            ylabel('GFP');
            set(gca,'Box','off','Color',[0.95 0.95 0.95]);
            hold on
            grid on
        end
        
        p1
        L3=legend(p1,cfg.Condnom);
        set(L3,'Box','off','Location','NorthWest');
    end
    
end

f3=figure; set(f3,'Color',[1 1 1]);
GMD=zeros(length(Time),size(CondGA,2));
AX=cell(size(CondGA,2),1); H1=zeros(size(CondGA,2),1); H2=zeros(size(CondGA,2),1);

for fcnt=1:size(CondGA,2)
    
    [SCorr,GMD(:,fcnt),~,~]=GMD_calc(CondGA{1,fcnt}, CondGA{1,fcnt},64,1:64,[Time(1) Time(end)],[Time(1) Time(end)],Time);
    
    subplot(2,1,fcnt)
    [AX{fcnt,1},H1(fcnt),H2(fcnt)]=plotyy(Time,GMD(:,fcnt),Time,gfp_cond(:,fcnt));
    set(get(AX{fcnt,1}(1),'YLabel'),'String','Global Map Dissimilarity (GMD)');
    set(get(AX{fcnt,1}(2),'YLabel'),'String','Global Field Power (GFP)');
    h=vline(0,'k','T0');
    title(horzcat('GMD and GFP : Condition ',cfg.Condnom{1,fcnt}));
    
    
end





end

function plotsingle_pe(hdl,~)

D=get(hdl,'UserDat');
time=D{1,1};
x1=D{1,2};
x2=D{1,3};
x1mean=D{1,4};
x2mean=D{1,5};
colors=D{1,6};
eoi=D{1,7};
gfp1=D{1,8};
gfp2=D{1,9};
gmd=D{1,10};
condnames=D{1,11};
ylims=D{1,12};
alldata=D{1,13};
chaninfo=D{1,14};

f1=figure; set(f1,'PaperUnits','normalized','PaperPosition',[0.0235308 0.0272775 0.894169 0.909249],'Color',[1 1 1]);
orient portrait; axis ('normal');

sh1=subplot(3,1,1);
pos1=get(sh1,'Position');
pos1=[pos1(1)-0.05 pos1(2)-0.02 pos1(3)+0.05 pos1(4)+0.025];
set(sh1,'Position',pos1);
pos1b=get(sh1,'Position');

[axs1,ploth1]=plot_ConfInt2(x1,time',x1mean,[colors(1,:)], 0.25,0.05);
set(axs1,'YDir','reverse');
hold on
[axs2,ploth2]=plot_ConfInt2(x2,time',x2mean,[colors(2,:)], 0.25,0.05);
set(axs2,'YDir','reverse','XTickLabel',[])
hold on
ploth3=plot(time,x2mean-x1mean,'r');
pvalue_corr=plot_perm({x1 x2},time',(ylims(2)-1),[time(1) time(end)]);
title(strcat('ERP signal with 95% CI:  ',char(condnames{1,1}),' &  ', char(condnames{1,2}),' : Electrode :', char(eoi)));
L2=legend([ploth1 ploth2 ploth3],condnames{1,1},condnames{1,2},strcat(condnames{1,2},'-',condnames{1,1}));
set(L2,'Box','off');
h=vline(0,'k','T0');


sh2=subplot(3,1,2,'Position',[pos1b(1) (pos1b(2)-pos1b(4)/1.5) pos1b(3) pos1b(4)/1.5]);
p1=plot(time,gfp1,'b');   %plot the difference if necessary
set(p1,'LineWidth',2);
hold on
p2=plot(time,gfp2,'g');
set(p2,'LineWidth',2);
grid on
ylabel('GFP');
set(gca,'Box','off','Color',[0.95 0.95 0.95],'XTickLabel',[])
L3=legend([p1 p2],condnames{1,1},condnames{1,2});
set(L3,'Box','off','Location','NorthWest');
pos1c=get(sh2,'Position');

sh3=subplot(3,1,3,'Position',[pos1c(1) (pos1c(2)-pos1c(4)/1) pos1c(3) pos1c(4)/1]);
p3=plot(time,gmd(1:length(time)),'k');   %plot the difference if necessary
set(p3,'LineWidth',2);
grid on
ylabel('GMD');
set(gca,'Box','off','Color',[0.95 0.95 0.95])
xlabel('Time (ms)');

%% OPEN DIALOGUE BOX TO DEFINE THE TIME INTERVALS TO SHADE

name='Define ERP & Topo Figure Parameters';
prompt={'Enter time interval for plotting topographies:','Enter lower limits of shaded areas (ms):', 'Enter upper limits of shaded areas (in ms):'};
numlines=[1 6 6]';
answr=inputdlg(prompt,name,numlines);
intval=str2double(answr{1});
lowlims=str2double(answr{2});
hilims=str2double(answr{3}); 
assignin('base','answr',answr); 

%% PLOT ERP WITH CONFIDENCE INTERVALS AGAINST TOPOGRAPHIES OF CONDITION DIFFERENCE

numtops=floor((time(end)+abs(time(1)))/intval);  %number of topographies
maplims=[-4 4];
lim_low=(ceil((abs(time(1))/100))*sign(time(1)))*100;
lim_hi=(ceil((abs(time(end))/100))*sign(time(end)))*100; 
ylimits=[-15 15];

f2=figure; set(f2,'PaperUnits','normalized','PaperPosition',[0.0235308 0.0272775 0.994169 0.909249],'Position',[586 512 1156 586],'Color',[1 1 1]);
orient portrait; axis ('normal'); 

sh1=subplot(3,numtops,1:numtops);
pos1=get(sh1,'Position');
% pos1=[pos1(1)-0.05 pos1(2)-0.02 pos1(3)+0.05 pos1(4)+0.055];
% set(sh1,'Position',pos1);
% pos1b=get(sh1,'Position');

[axs1,ploth1]=plot_ConfInt2(x1,time',x1mean,[colors(1,:)], 0.25,0.05);
set(axs1,'YDir','reverse');
set(axs1,'YDir','reverse','XLim',[lim_low lim_hi],'YLim',ylimits,'XTick',[lim_low:intval:lim_hi]);
hold on
[axs2,ploth2]=plot_ConfInt2(x2,time',x2mean,[colors(2,:)], 0.25,0.05);
set(axs2,'YDir','reverse','XLim',[lim_low lim_hi],'YLim',ylimits,'XTick',[lim_low:intval:lim_hi]);
ylabel('Potential(\muV)'); 
hold on

% Add shading over the intervals defined in the dialogue box
shadcols={[1 1 0.7],[0.7 1 1],[1 0.8 0.8],[0.8 0.8 1],[0.8 1 0.8],[10.8 1],[1 1 0.7],[0.7 1 1],[1 0.8 0.8],[0.8 0.8 1],[0.8 1 0.8],[10.8 1]}; 
for icnt=1:length(lowlims)
yintval=ones(1,length( lowlims(icnt):hilims(icnt)))*ylimits(1);
ha=fill([lowlims(icnt) lowlims(icnt):hilims(icnt) hilims(icnt)],[ylimits(2) yintval ylimits(2)],shadcols{1,icnt},'FaceAlpha',0.5,'EdgeColor','none'); 
assignin('base','ha',ha); 
hold on
end

ploth3=plot(time,x2mean-x1mean,'r');
title(strcat('ERP signal with 95% CI:  ',char(condnames{1,1}),' &  ', char(condnames{1,2}),' : Electrode :', char(eoi)));
L2=legend([ploth1 ploth2 ploth3],condnames{1,1},condnames{1,2},strcat(condnames{1,2},'-',condnames{1,1}));
set(L2,'Box','off');
h=vline(0,'k','T0');

%% PLOT THE TOPOGRAPHIES 
xdiff=alldata{1,2}-alldata{1,1};
tcurr=time(1);
sh2=zeros(numtops,1);
for cnt=1:numtops
    sh2(cnt)=subplot(3,numtops,cnt+numtops);   %subplot for topographies
    i=find(time>=tcurr & time<=tcurr+intval);
    display(time(i)); 
    if cnt==numtops
        topoplot(mean(xdiff(:,i),2),chaninfo,'maplimits',maplims,'electrodes','on');
        colorbar; 
    else
        topoplot(mean(xdiff(:,i),2),chaninfo,'maplimits',maplims,'electrodes','on');
    end
    tcurr=tcurr+intval;
end

set(f2,'Color',[1 1 1]); 







end
















