function [GFP,i]=GFP_calc(varargin)
%************************************************************************************************************
% Date: 16-02-2015         Programmer: D. BOLGER
%
% Function to calculate the Global Field Power for a given map defined by a
% vector of time limits. 
% The Global Field Power (GFP) is the standard deviation of the potentials
% at all electrodes of an average reference map.
% It is used in spatial analysis of EEG topographic data and is adapted
% from the GFP calculation implemented in the Cartool program (Brunet et
% al, 2010). 
% Plotting GFP over time allows one to detect moments of high
% signal-to-noise ration which correspond with moments of high global
% neuronal synchronisation (Brunet et al, 2010). 
% Inputs: DataIn ==> the input data matrix.
%            Enum ==> total number of electrodes
%            Eind ==> the indices of the electrodes used for GFP
%            calculation
%            TLim ==> vector indicating time limits in ms [time_low
%            time_high]
%            Time=> time vector
% Output : GFP => calculated GFP for each input time interval
%**************************************************************************************************************
disp('THIS ONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
varargin=varargin';
if nargin ==1
    DataIn=varargin{1,1};
    i=1:size(DataIn,2);
    Enum=size(DataIn,1);
    Eind=[1:Enum];
    
elseif nargin>1
    
    DataIn=varargin{1,1};
    Enum=varargin{2,1};
    Eind=varargin{3,1};
    TLim=varargin{4,1}
    Time=varargin{5,1}
    
    i=find(Time>=TLim(1) & Time<=TLim(2));   %find the time indices of the current time interval
    
end

X=zeros(Enum,1);
umean=zeros(length(i),1); 
GFP=zeros(length(i),1);    %initialise the GFP vector

for tcount = 1:length(i)   %for each time interval
    
    umean(tcount)=sum(DataIn(Eind,i(tcount)),1)/Enum;   %calculate the mean voltage of all electrodes of the current time point (map) -the average reference
    
    for ecount=1:Enum       %for each electrode
        
        X(ecount)=(DataIn(ecount,i(tcount)) - umean(tcount))^2;   %Ui - umean
        
    end
    
    GFP(tcount)=sqrt(sum(X)/Enum);   % calculation of the GFP for the current time interval 
    
end
    
    
   
    