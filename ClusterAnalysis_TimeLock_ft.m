% Date: 01-02-1976   Programmed by: Deirdre Bolger
% Script to carry out a non-parametric statistic test on ERP data. It uses
% a cluster-based test statistic to deal with the multiple comparisons
% problem (MCP).
% The main steps of the test are as follows :
% 1. The experimental conditions are compared for every channel-time pair by
% means of a t-test.
% 2. All those samples whos t-value exceeds a predefined
% threshold are selected.
% 3.The selected samples are clustered in connected sets on the basis of temporal and spatial adjacency.
% 4. The sume of the t-values within each cluster are taken which yields
% the cluster-level statistics.
% 5. The maximum of the cluster-level statistics are found, which gives the
% test statistic used to evaluate the effect of the experimental
% conditions.
% As the MCP is resolved by the cluster-based test statistic, the
% significance probability has to be calculated using the Monte Carlo
% method.
% Very important: Need to know if the structure of your experiment falls
% under the heading of a "between-" or a "within" unit of observation
% design. The Monte Carlo significance probability is calculated
% differently depending on the experimental structure (within or between UO
% design).

clear all;
close all;


%% Define the condition names and groups and root directory *************************
Conds=  {'Control' 'Relie'};                             
Groups ={'sanschange' 'sanschange'};
dir_root='F:\BLRI\EEG\Projets_en_cours\projet_VariAcoustique\Conditions\MotHF\'; 
Time_Interval=[0 1]; %time interval over which to carry out the analysis
SLOption='No';      %choose to calculate the surface laplacian
channum=64;

%% **************************************Load in the data files and convert to Fieldtrip form for spatial permutation analysis*****************************************

if strcmp(Groups{1,1},Groups{1,2})==1
    lenGroup=1;
else
    lenGroup=length(Groups);
end

AllFiles=cell(lenGroup,size(Conds,2));
AllDirs=cell(lenGroup,size(Conds,2));
DataIn_eeglab=cell(lenGroup,size(Conds,2)); 

for Gcounter=1:lenGroup
    for Condcounter=1:size(Conds,2)
        
        [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;      
        
  
        AllDirs{Gcounter,Condcounter}=strcat(dir_root,char(Groups{1,Gcounter}),'\',char(Conds{1,Condcounter}),'\');  %Current directory
        AllFiles{Gcounter,Condcounter}=dir(strcat(AllDirs{Gcounter,Condcounter},'*.set'));   %declare empty variable for files
        currfiles=AllFiles{Gcounter,Condcounter};
        filenom={AllFiles{Gcounter,Condcounter}.name};
        
        EEG=pop_loadset('filename',filenom,'filepath',AllDirs{Gcounter,Condcounter});
        [ALLEEG, EEG, CURRENTSET]=pop_newset(ALLEEG,EEG,0,'study',0);
        EEG=eeg_checkset(EEG); eeglab redraw; 
        
        DataIn_eeglab{Gcounter,Condcounter}=ALLEEG; 
        numsuj=length(ALLEEG); 
    end
end

DataIn=cell(numsuj,size(DataIn_eeglab,1),size(DataIn_eeglab,2));
DataIn_pp=cell(numsuj,size(DataIn_eeglab,1),size(DataIn_eeglab,2));

for ccnt=1:size(DataIn_eeglab,2)
    for gcnt=1:size(DataIn_eeglab,1)
        for fcnt=1:numsuj
            DataIn_pp{fcnt,gcnt,ccnt}= eeglab2fieldtrip(DataIn_eeglab{gcnt,ccnt}(fcnt), 'preprocessing', []);
            
            alleegs=DataIn_eeglab{gcnt,ccnt}(fcnt); 
            cfg=[];
            cfg.channel={alleegs.chanlocs(1:64).labels};
            cfg.trials='all';
            cfg.keeptrials='no';
            cfg.covariance='no';
            
            DataIn{fcnt,gcnt,ccnt}=ft_timelockanalysis(cfg, DataIn_pp{fcnt,gcnt,ccnt});
        end
    end
end


%% OPTION TO CALCULATE THE SURFACE LAPLACIAN

if strcmp(SLOption,'Yes')==1
    
    disp('Calculating the Surface Laplacian')
    
    csdFileName='C:\Users\bolger\Documents\MATLAB_tools\CSDtoolbox\CSDtoolbox\resource\10-5-System_Mastoids_EGI129.csd';
    Labels=textread('C:\Users\bolger\Documents\MATLAB_tools\CSDtoolbox\CSDtoolbox\resource\E64_10_20.txt','%s');
    
    montage=CREx_ExtractMontage(csdFileName,Labels);    %Extract the Montage
    [Gout,Hout]=CREx_CalcGH(montage,5);                          %Calculate the G and H transformation matrices
    
    SL=cell(size(DataIn));
    CSD=cell(size(DataIn));
    
    for GCount=1:size(DataIn,2)
        for SCount=1:size(DataIn,1)
            
            [CSD{SCount,GCount},SL{SCount,GCount}]=CREx_CalcCSD(DataIn{SCount,GCount}.avg(1:64,:),Gout,Hout,64);
            DataIn{SCount,GCount}.avg=SL{SCount,GCount};
            %DataIn{SCount,GCount}.SL='yes';
            
        end
    end
end

%% Compute non-parametric statistical test


if lenGroup>1
    DEV={DataIn{:,:,1}};
    STD={DataIn{:,:,2}};
    DEV_char=Conds{1,1};
    STD_char=Conds{1,2};
elseif lenGroup==1
    DEV={DataIn{:,2}};
    STD={DataIn{:,1}};
    DEV_char=Conds{1,2};
    STD_char=Conds{1,1};
end


%% CALCULATE THE CLUSTER-BASED NON-PARAMETRIC TEST

% Calculate the neighbours for current electrode configuration
cfg=[];
cfg.layout ='C:\Users\bolger\Documents\MATLAB_tools\fieldtrip-master\template\layout\biosemi64.lay' ;
lcfg.rotate=90;
layout1 = ft_prepare_layout(cfg,DataIn{1,1});
figure
ft_plot_lay(layout1,'label','yes','point','no','outline','yes','box','no');

% Define the neighbours
cfg=[];
cfg.method= 'triangulation';
%cfg.neighbourdist = 5;   %define the maximum distance between neighbouring sensors
cfg.layout = layout1;
cfg.channel={DataIn{1,1}.label{1:channum}};
neighbs=ft_prepare_neighbours(cfg,DataIn{1,1});

cfg=[];
cfg.neighbours=neighbs;
cfg.enableedit='no';
cfg.layout=layout1;
ft_neighbourplot(cfg)

cfg=[];
cfg.channel = {EEG.chanlocs(1:channum).labels};
cfg.latency=Time_Interval;
cfg.avgoverchan = 'no';
cfg.avgovertime= 'no';
cfg.parameter='avg';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT'; % independent samples T-statistic to measure effect size at sample level.
cfg.correctm = 'cluster';                          % method to deal with multiple comparisons
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail = 0;
cfg.clustertail=0;
cfg.correcttail='prob';
cfg.alpha=0.025;
cfg.minnbchan = 2;                  %the minimum number of neighbourhood channels required for selected channel to be included in the clustering algorithm
cfg.numrandomization = 1000;  %number of draws from the permutation distribution
cfg.neighbours=neighbs;

% Create the design matrix to store information about the UOs --for
% MotInter

design=zeros(2,size(DataIn,1)*size(DataIn,2)*lenGroup);  %initialise the design

for gcount=1: lenGroup
    starter=(size(DataIn,1)*length(Conds))*(gcount-1);
    design(1,starter+1:size(DataIn,1)+starter)=1:size(DataIn,1);
    for ccounter=1:length(Conds)-1
        starter1=(size(DataIn,1)*length(Conds))*(gcount-1);
        design(1,((size(DataIn,1)*ccounter)+1)+starter1:(size(DataIn,1)*(ccounter+1))+starter1)=1:size(DataIn,1);
    end
end
clear starter starter1 gcount

% second line of design
hdesign=size(design,2)/length(Groups);
for gcount=1:length(Groups)
    starter=gcount+((hdesign-1)*(gcount-1))
    design(2,starter:hdesign*gcount)=gcount;
end

display(design);
cfg.design=design;


cfg.ivar=2;                   %specifies the row of the design matrix corresponding to the independent variable
cfg.uvar=1;                  %specifies the row of the design matrix corresponding to the unit of observation ...here the subjects
tstats_all=ft_timelockstatistics(cfg,DEV{:},STD{:})   %DEV{:},STD{:}


%% Calculate the difference between the two experimental conditions

cfg=[];
cfg.latency=[tstats_all.time(1) tstats_all.time(end)]
avgDEV=ft_timelockgrandaverage(cfg,DEV{:})   %, DataIn{:,2});   %Calculate the grand average over all subjects
avgSTD=ft_timelockgrandaverage(cfg,STD{:})  %, DataIn{:,4});

cfg=[];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
diff_DEVvsSTD_all=ft_math(cfg,avgDEV,avgSTD);

%% Construct a boolean matrix indicating which sensors are members of the significant clusters
% Construct two matrices: one for positive and one for negative clusters
% Begin with the positive clusters

if isempty(tstats_all.posclusters)==0
    pos_cluster_tstatsall= [tstats_all.posclusters(:).prob];
    pos_sig=find(pos_cluster_tstatsall<=tstats_all.cfg.alpha);
    
    if isempty(pos_sig)==0
        pos=ismember(tstats_all.posclusterslabelmat,pos_sig);      %construct a boolean matrix of which channel-time pairs are part of significant clusters
    elseif isempty(pos_sig)==1
        
        display('No significant positive clusters!')
    end
elseif isempty(tstats_all.posclusters)==1
    
    display('No positive clusters!');
    pos_sig=[];
end


%Now do the negative clusters
if isempty(tstats_all.negclusters)==0
    neg_cluster_tstatsall=[tstats_all.negclusters(:).prob];
    neg_sig=find(neg_cluster_tstatsall<=tstats_all.cfg.alpha);
    if isempty(neg_sig)==0
        neg=ismember(tstats_all.negclusterslabelmat,neg_sig)
    elseif isempty(neg_sig)==1
        
        display('No significant negative clusters!')
    end
elseif isempty(tstats_all.negclusters)==1
    
    display('No negative clusters!');
    neg_sig=[];
end


%% Plot the Results of the Non-Parametric Cluster-based Test
figure;
set(gcf,'Color',[1 1 1])
D=diff(tstats_all.time);     %resolution is generally in the order of 2ms.
timestep=D(1)*25;
sample_count=length(tstats_all.time);
SR=DataIn_pp{1,1,1}.fsample;
j=[tstats_all.time(1):timestep:tstats_all.time(end)];    %[DataIn{1,1}.time(1):timestep:tstats_all.time(end)]
m=[1:floor(timestep*SR):sample_count]

for k=1:length(m)-1

    subplot(ceil(length(m)/4),4,k)
    cfg=[];
    cfg.xlim=[j(k) j(k+1)];
    cfg.zlim=[-2 2];
    
    if isempty(pos_sig)==0 && isempty(neg_sig)==0
        disp('**************Some significant stuff: positive and negative*********************************');
        pos_int=all(pos(:,m(k):m(k+1)),2);
        neg_int=all(neg(:,m(k):m(k+1)),2);
        cfg.highlightchannel=find(pos_int | neg_int);
    elseif isempty(pos_sig)==0 && isempty(neg_sig)==1
        display('No significant negative clusters!');
        pos_int=all(pos(:,m(k):m(k+1)),2);
        cfg.highlightchannel=find(pos_int);
    elseif isempty(pos_sig)==1 && isempty(neg_sig)==0
        display('No significant positive clusters!');
        neg_int=all(neg(:,m(k):m(k+1)),2);
        cfg.highlightchannel=find(neg_int);
    elseif isempty(pos_sig)==1 && isempty(neg_sig)==1
        cfg.highlightchannel=[];
        display('*****************No significant anything!*************************');
    end
    
    
    
    cfg.highlight='labels';
    cfg.comment='xlim';
    cfg.commentpos='title';
    cfg.layout = layout1;
    cfg.zlim=[-3 3]
    ft_topoplotER(cfg,diff_DEVvsSTD_all)
    
    if k==1
        colorbar('location','westoutside')
    end
    
end

%% PLOT THE RESULTS OF CLUSTER-BASED NON-PARAMETRIC TEST ON INDIVIDUAL ELECTRODES

%reuse the design defined above.

cfg1=[];
cfg1.method='montecarlo';
cfg1.statistic= 'ft_statfun_depsamplesT';
cfg1.correctm='cluster';
cfg1.clusteralpha=0.05;
cfg1.clustertail=0;
cfg1.tail=0;  % two-sided test
cfg1.correcttail='prob';
cfg1.alpha=0.025;
cfg1.numrandomization=1000;
cfg1.design=design;
cfg1.ivar=2;    %independent variable
cfg1.uvar=1;   %unit of observation
cfg1.channel={DataIn{1,1}.label{1:64}};
cfg1.neighbours=[];

tstats_single=ft_timelockstatistics(cfg1,DEV{:},STD{:});

f1=figure; set(f1,'Color',[1 1 1]);
s=zeros(length(cfg1.channel),1);
st=zeros(length(cfg1.channel),1);

for chan_cnt=1:length(cfg1.channel)
    
    pcurr=tstats_single.prob(chan_cnt,:);
    i=find([pcurr>cfg1.alpha]);
    i2=find([pcurr<=cfg1.alpha]);
    pcurr(i)=0;
    pcurr(i2)=10;
    s(chan_cnt)=subplot(8,8,chan_cnt);
    st(chan_cnt)=stem(s(chan_cnt),tstats_single.time,pcurr);
    set(st(chan_cnt),'MarkerSize',2,'Color','g');
    
    grid on
    xlabel('time (seconds)');
    ylim([0 1])
    xlim([-0.2 1])
    title(s(chan_cnt),DataIn{1,1}.label{chan_cnt});
    
    set(s(chan_cnt),'HitTest','on','SelectionHighlight','on','UserData',{avgSTD.avg(chan_cnt,:),avgDEV.avg(chan_cnt,:),avgSTD.time,tstats_single.time,pcurr,STD_char,DEV_char,DataIn{1,1}.label{chan_cnt},i2},'NextPlot','add');
    set(s(chan_cnt),'ButtonDownFcn',@plotsingle1)
    disp(get(s(chan_cnt),'UserData'));
    
end











