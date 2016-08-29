function [SC,GMD,GFP1,GFP2]=GMD_calc(DataIn1, DataIn2,Enum,Eind,TLim1,TLim2,Time)
%*************************************************************************************************
% Date: 15-2-2015                 Programmer: D.BOLGER
% This function calculates the Global Map Dissimilarity (GMD) between two
% maps. The GMD is a measure of topographic differences of scalp potential
% maps. GMD calculation is an initial step in defining whether the cortical
% activity at the scalp for the two datasets being compared was generated
% by different sources. 
%The maps can based on voltage values from different datasets
% (different condition, population) or the same dataset but different time
% interval or two successive time points. 
% This function calls the GFP_calc function as the two maps being compared
% are first normalised by their respective GFP values; the potential values
% at each electrode in each map are divided by the map's GFP. This ensures
% that only topographical differences are taken into account when
% calculating the GMD. 
% Input: DataIn1 ==> input data matrix 1 (map 1)
%          DataIn2 ==> input data matrix 2 (map 2)
%          Enum ==> Number of electrodes
%          Eind ==> indices of the electrodes.
%          TLim1 ==> upper and lower limits of time interval of map 1 (ms)
%          TLim2==> upper and lower limits of time interval of map 2 (ms)
%          Time ==> time vector. 
% Output: GMD==> GMD of the two input maps. 
%             SC ==> Spatial Correlation (C) - the Pearson cross
%             correlation coefficient. 
%**************************************************************************************************

%% Verify that the number of channels and time frames of the two input maps

if size(DataIn1,2)~=size(DataIn2,2)
    error('Number of timeframes is different');
    return; 
elseif size(DataIn1,1)~=size(DataIn2,1)
    error('Number of channels is different');
    return;
end

%% CHECK THAT THE TIME DIMENSION OF INPUT DATA MATRIX IS SAME LENGTH AS TIME VECTOR

ind_std=find(Time>=TLim1(1) & Time<=TLim1(2));
ind_dev=find(Time>=TLim2(1) & Time<=TLim2(2));

%% Configure the input datasets if the two input datasets are the same (GMD
%calculation over time

     if isequal(DataIn1,DataIn2)==1
         
        disp('The two datasets are the same');
     
       DataIn2=DataIn2(:,2:end);
       DataIn1=DataIn1(:,1:end-1);
       T=Time(1:length(Time)-1);
     
     else
         
         disp('The two datasets are different.');
         T=Time;
     end
    
%% Normalise the two maps by dividing the values at each electrode by map GFP.


GFP1=GFP_calc(DataIn1, Enum, Eind,TLim1, T);   %calculate the GFP for each time point of map1
GFP2=GFP_calc(DataIn2, Enum, Eind, TLim2, T);  %calculate the GFP for each time point of map2
size(GFP1)
size(GFP2)

norm1=repmat(GFP1', Enum,1);        %  Not sure if it is necessary to do that here..
norm2=repmat(GFP2',Enum,1);
size(norm1)
size(norm2)

% DataIn1_norm=DataIn1(Eind,ind_std)./norm1;  %normalise map1
% DataIn2_norm=DataIn2(Eind,ind_dev)./norm2;  %normaise map2

%% Calculate the GMD of the input maps.

i1=find(T >=TLim1(1) & T<=TLim1(2));  % Find the indices of the time interval of map 1
i2=find(T >=TLim2(1) & T<=TLim2(2));

U=zeros(Enum,1);
V=zeros(Enum,1);
DiffUV=zeros(Enum,1);
GMD=zeros(length(i1),1); 
SC=zeros(length(i1),1);

% We assume that the two intervals are of the same length (same number of
% indices)
for Tcount=1:length(i1)
    
    M1=DataIn1(:,Tcount);                  %DataIn1(Eind,i1(Tcount)); 
    M2=DataIn2(:,Tcount);                  %DataIn2(Eind,i2(Tcount));
    
    umean=sum(M1,1)/Enum;
    vmean=sum(M2,1)/Enum;
    
    for Ecount=1:Enum
        
      U(Ecount)=(M1(Eind(Ecount))-umean)/GFP1(Tcount);
      V(Ecount)=(M2(Eind(Ecount))-vmean)/GFP2(Tcount)';
      DiffUV(Ecount)=(U(Ecount)-V(Ecount))^2;
      
    end
    
    GMD(Tcount)=sqrt(sum(DiffUV)/Enum);
    SC(Tcount)=1-((GMD(Tcount)^2)/2); 
end 
      
    
    GMD=cat(1,0,GMD);
    SC=cat(1,0,SC);
end  









