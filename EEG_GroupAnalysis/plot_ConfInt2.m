function[haxe,h1]= plot_ConfInt2(varargin)
% Programmer: Deirdre Bolger
% Input Variables:
%                         DataIn = the input matrix of data;
%                                      no. data points X no. subjects
%                         points
%                         T= time vector (vertical)
%                         Xin_mean = mean data vector over all subjects
%                         ci = optional; confidence interval. If the confidence
%                         interval is not specified, the default is 0.05
% A subfunction SEM_calc is included to calculate the standard error of the
% mean
%***************************************************************************************
if nargin ==5
    Xin=varargin{1,1};
    t=varargin{1,2};
    Xin_mean=varargin{1,3};
    dath=varargin{1,4};
    alf=varargin{1,5};
    ci = 0.05;       %specify confidence interval
    
elseif nargin ==7
    Xin=varargin{1,1};
    t=varargin{1,2};
    Xin_mean=varargin{1,3};
    dath=varargin{1,4};
    alf=varargin{1,5};
    ci=varargin{1,6};
    cond_title=varargin{1,7};  %title for current figure
    
elseif nargin ==6
    Xin=varargin{1,1};
    t=varargin{1,2};
    Xin_mean=varargin{1,3};
    dath=varargin{1,4};
    alf=varargin{1,5};
    ci=varargin{1,6};
    
    
    
elseif nargin ==4
    Xin=varargin{1,1};
    t=varargin{1,2};
    Xin_mean=varargin{1,3};
    dath=varargin{1,4};
    alf=0.25;
    ci = 0.05;
    
elseif nargin < 3
    disp('Not enough arguments!');
    return
elseif nargin > 4
    disp('Too many arguments!');
    return
end

disp(dath)

sem =SEM_calc(Xin',ci);  % call of function to calculate the standard error of mean
assignin('base','SEM',sem);
assignin('base','Xin_mean',Xin_mean);

X=[t', fliplr(t')];

Y=[(Xin_mean'-sem),fliplr((Xin_mean'+sem))];


f1=fill(X,Y,dath,'FaceAlpha',alf,'LineStyle','none');
set(f1,'edgecolor','none');
alpha(0.25);
hold on;
h1=plot(t,Xin_mean,'Color',dath,'LineStyle','-','LineWidth',2.5);
set(gca,'Xlim',[t(1) t(end)]);
set(gca,'Ylim',[-15 15]);
set(gca,'XGrid','on','YGrid','on','Box','off');


if length(varargin)==7
    set(get(gca,'Title'),'String',cond_title,'FontSize',10);
elseif length(varargin)<7
    set(get(gca,'Title'),'String','Condition: ','FontSize',10);
end

haxe=gca;
end

function sem=SEM_calc(vect, CI)
% SEM_calc - standard error of the mean, confidence interval
%
% function sem=SEM_calc(vect, CI)
%
%***********************************************************************

error(nargchk(1,2,nargin))

if isvector(vect)
    vect=vect(:);
end


if nargin==1
    stdCI = 1.96 ; % if the interval specified is 5%
elseif nargin==2  %if a different interval is specified
    CI = CI/2 ; %Convert to 2-tail
    stdCI = abs(norminv(CI,0,1)) ;  %uses inverse of the normal cdf.
end

sem = ( (std(vect)) ./ sqrt(sum(~isnan(vect))) ) * stdCI ;    % ignores nans
end
