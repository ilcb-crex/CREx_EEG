function CREx_plottopo_config()

h=figure; set(h,'PaperUnits','normalized','PaperPosition',[0.0235308 0.0272775 0.894169 0.909249]);
orient portrait; axis ('normal')

xvals=zeros(length(EEG.chanlocs),1);
yvals=zeros(length(EEG.chanlocs),1);
pwidth    = 0.95;     % 0.75, width and height of plot array on figure
pheight   = 0.88;
axwidth=0.05;
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

 Axes = [];
for c=1:64

    xcenter=xvals(c);
    ycenter=yvals(c);
    Axes=[Axes axes('Units','normalized','Position', [xcentre-axwidth/2 ycenter-axheight/2 axwidth axheight])]; 
    axis('off')     
    hold on;
    
    plot(
    
end 













end