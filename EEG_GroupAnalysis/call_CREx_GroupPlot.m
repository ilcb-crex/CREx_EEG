% Date: December 2015                Programmed by : D. Bolger
%
% Setting the parameters of the Group Plot configuration structure:
% GPcfg.Condnom should be a cell array of 2 or 3 strings corresponding to
% the condition names.
% GPcfg.TimeWind can be a two-element vector [min max] or 'all' to include
% all time points. 
% This script calls the function CREx_GroupPlot()
% The function CREx_GroupPlot() plots ERPs for two conditions for the defined channels (eoi),
% calculates the permutation test over the time window defined (TimeWind)
% and plot this agains the global field power (GFP) and the global map
% dissimilarity (GMD) calculated between the two conditions. 
% NOTE: 
% - If more than two conditions are defined, the permutation test will
% not be carried out. 
% - To ensure that this function runs, conditions need to be in a folder
% whose name is the condition name defined in GPcfg.Condnom. Any
% sub-groups should be contained in folders with the same title as that
% defined in GPcfg.Condnom1. The function CREx_GroupPlot() will load all
% *.set files found in the defined search path, GPcfg.dirs. 
%****************************************************************************************************************************

clear all
close all

AllChans=load('C:\Users\bolger\Documents\MATLAB_tools\DeeScripts\CREx_GroupAnalysis\ChanInfo.mat'); % Load the channel info file.

GPcfg=[];
GPcfg.currdir='F:\BLRI\EEG\Projets_en_cours\Projet_MotInter\ExpEEG_Phase1\Data_Biosemi\P1AUD_Results\'    % The current directory in which the condition folders can be found...needs to be changed
GPcfg.Condnom1= {'MotsCC','MotsCC'};                                                                                                          % 'MAct' 'Num' 'Trans' Declare the conditions to be compared; need to have same names as condition folders.
GPcfg.Condnom2={'NOL' 'DYS'};                                                                                                                 % If there are subgroups...
GPcfg.Subnum=24;                                                                                                                                      % Number of subjects.
GPcfg.Channum=64;                                                                                                                                    % The number of channels on which to base the analysis.
GPcfg.TimeWind='all';                                                                                                                           % The time window of interest (in ms) or 'all' to carry out analysis on entire epoch. 
GPcfg.eoi=AllChans.Chaninfo;                                                                                                         %  Define the electrodes of interest.
GPcfg.ylims=[-10 10];

if isempty(GPcfg.Condnom2)==0
    GPcfg.Condnom={strcat(GPcfg.Condnom1{1,1},'-',GPcfg.Condnom2{1,2}),strcat(GPcfg.Condnom1{1,2},'-',GPcfg.Condnom2{1,2})};
    GPcfg.dirs={strcat(GPcfg.currdir,GPcfg.Condnom1{1,1},'\',GPcfg.Condnom2{1,1},'\') strcat(GPcfg.currdir,GPcfg.Condnom1{1,2},'\',GPcfg.Condnom2{1,2},'\')}
elseif isempty(GPcfg.Condnom2)==1
    GPcfg.Condnom=GPcfg.Condnom1;
    GPcfg.dirs={strcat(GPcfg.currdir,GPcfg.Condnom1{1,1},'\cleaned\') strcat(GPcfg.currdir,GPcfg.Condnom1{1,2},'\cleaned\')};
end

[AllCond, GFPcond, perm_pval,condData]=CREx_GroupPlot_v2(GPcfg);                                                                % Call of function CREx_GroupPlot to plot ERPs with results of permutation test against GFP and GMD



