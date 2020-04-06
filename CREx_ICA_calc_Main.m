%% CALL OF FUNCTION TO CARRY OUT ICA ON CONTINUOUS DATA
% It uses functions from the ADJUST toolbox to detect thos ICs that
% correspond to artifacts.
% It can also apply a PCA before carrying out ICA to reduce the number of components and speed up ICA computation for continuous data.
% The number of PCA components is calculated automatically based on
% explained variance (99% explained variance).
%% 

paramfile_nom = 'parameters-L2EEGVR.txt';  % The title of the parameters file.
% The path to the parameters path. This should be the only path that needs
% to be changed.
paramfile_path = fullfile(filesep,'Users','bolger','Documents','work','Projects','Projet-L2-VREEG','Data',paramfile_nom); % Put the path to your parameters file here.

fid2 = fopen(paramfile_path);  % Define a file identifier.
mydata = textscan(fid2,'%s %s');  % Scan textfile...

for i = 1:length(mydata{1,1})                     % generate a parameters structure from the parameters text file
    Params.(genvarname(mydata{1,1}{i})) = mydata{1,2}(i);
end
currsuj = 'S4MMP';

%% Define Options 

dopca = 0; % or 0 if not doing PCA before ICA.
doICA = 0; % or 0 if not doing ICA calculation (and only rejection of components).
doRej = 1; % or 1 if only rejecting components and carrying out back-projection of retained components
doVis = 1; % or 0 if you do not want to visualise and data before and after IC rejection.

CREx_ICA_calc(dopca, doICA, doRej,doVis, Params, currsuj);

