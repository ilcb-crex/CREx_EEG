%%Date: 09-04-2019 .     Programmed by: D. Bolger
% Function to carry out ICA on a selected EEG dataset.
% Updated March/April 2020.
% To run ICA, it requires that the ADJUST toolbox is installed.
% It also runs the ADJUST toolbox to assist the user in identifying ICA
% components that correspond to ocular and muscular artefacts.
%
% If dopca is 1, this implies carrying out a PCA before ICA; the number of
% PCA components will be determined automatically by calculating the
% explained variance and determining the number of PCs that explain 99%+ of
% the variance. A plot of the explained variance versus number of ICs is
% generated. 
% The script also allows the user to compare the EEG signals before and
% after IC rejection to confirm or other wise if those ICs should be
% retained or rejected.
%
% Inputs: dopca: 1 to do PCA before ICA, 0 not to carry out PCA.
%         isICA: 1 to carry out ICA, 0 not to carry out ICA.
%         isrej: 1 to back-project retained ICs, 0 not to carry out back
%         projection.
%         isVis: 1 to visualise the EEG signals before and after IC
% rejection.
%         Params: the parameters structure containing the paths to where
%         the data should be saved (Params.Datadir) and the channel
%         locations file (Params.Chanlocdir).
%         sujcurr: the title of the current subject (minus the *.set). 
%**************************************************************************

function CREx_ICA_calc(dopca, isICA, isrej, isVis, Params, sujcurr)


[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;                % Open eeglab session
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_loadset();           %Allows you to select and load the dataset manually.
icalen=length(EEG.chanlocs)-8; % Do not take into account the external electrodes.

pathcurr = fullfile(Params.Datadir{1,1},sujcurr,filesep);

allfiles = dir(pathcurr);
indx = ismember({allfiles(:).name},strcat(sujcurr,'-info.txt'));
txtnom= allfiles(indx).name;
fid = fopen(txtnom,'a+');

if isICA==1
    
    if size(EEG.data,3)>1          %If input data is epoched
        
        
        EEG = eeg_epoch2continuous(EEG);
        
    end
    
    %% Estimate the number of PCA components to define by calculating the explained variance.
    
    if dopca==1
        
        [coeffs, score, latent, tsquared, explained, ~] = pca(EEG.data(1:icalen,:));
        cumsum_exp = cumsum(explained); %find the cumulative sum of the explained variance for all the PCs.
        iex = dsearchn(cumsum_exp, 99);
        cumsum_dummy = nan(size(cumsum_exp,1),1);
        cumsum_dummy(iex)= cumsum_exp(iex);
        
        figure;
        plot(1:size(explained,1), cumsum_exp,'LineWidth',2)
        hold on
        plot(1:size(explained,1), cumsum_dummy,'r.','MarkerSize',36);
        ylabel('Explained Variance (%)')
        xlabel('Number of Components');
        set(gca,'XGrid','on','YGrid','on');
        
        %% CARRY OUT ICA USING THE RUNICA FUNCTION
        
        fprintf(fid,'\nCarrying out PCA with %4d components before ICA\n', iex);
        
        display('--------------Carrying out ica: patience!--------------------');
        [weights, sphere,compvars,bias,signs,lrates,activations]=runica(EEG.data(1:icalen,:), 'extended',1, 'pca',iex);
        
        
        EEG.icaweights = weights;
        EEG.icasphere = sphere;
        EEG.icachansind = 1:icalen;
        EEG.icaact = activations;
        EEG.icawinv = pinv(weights*sphere); % Pseudoinverse as matrix is not square if pca was applied first.
        
    else % In case of no prior PCA
        
        fprintf(fid,'\nCarrying out ICA with %4d components\n', icalen);
        
        disp('-------------------ICA: No prior dimensionality reduction using pca---------------------------\n');
        
        [weights, sphere,~,~,~,~,activations]=runica(EEG.data(1:icalen,:), 'extended',1);
        
        EEG.icaweights = weights;
        EEG.icasphere = sphere;
        EEG.icachansind = 1:icalen;
        EEG.icaact = activations;
        icaprojdata = icaproj(EEG.data(1:icalen,:),weights,1:icalen);
        EEG.icawinv = inv(weights*sphere);
        
    end
    
    %% SAVE THE DATASET WITH ICA COMPONENTS CALCULATED TO A NEW DATASET.
    
    ica_nom = strcat(EEG.setname,'-ica');
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(ica_nom),'gui','off'); % current set = xx;
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(ica_nom),'filepath',pathcurr);
    eeglab redraw
    
    %% RUN THE ADJUST TOOLBOX TO AUTOMATICALLY DETECT IC COMPONENTS CORRESPONDING TO ARTIFACTS
    
    report_nom = strcat(EEG.setname,'_icareport');
    report_out = fullfile(pathcurr,report_nom); % Adjust algorithm seems to work better when it computes ICAs itself.
    %EEG.icaact = activations;
    
    % Calculate component statistics.
    % This function creates a text file summarizing the results of the
    % statistical analyses and identifying components corresponding to
    % artifacts.
    
    [art, horiz, vert, blink, disc,soglia_DV, diff_var, soglia_K, med2_K, meanK,...
        soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD,soglia_GDSF, med2_GDSF,...
        GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin,max_var,activs] = ADJUST (EEG,report_out);
    
    %% Use the function from ADJUST toolbox to indicate those components that correspond to artifacts.
    % All components are displayed in a figure as topographies and those corresponding to
    % artifacts according to the statistics calculated in the ADJUST toolbox.
    % Click just below each topography to open a window showing the details of
    % the statistics calculated for that component.
    
    chanfile = Params.Chanlocdir{1,1};
    typecomp = 1;
    windhandle = nan;
    
    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    nrows = 5;
    ncols = ceil(size(EEG.icawinv,2)/nrows);
    
    for counter = 1:size(EEG.icawinv,2)
        
        ax1(counter) = subplot(nrows, ncols, counter);
        
        if sum(ismember(art,counter))==0
            currtitre = horzcat('Component ',num2str(counter));
            fontcol = 'k';
            
        elseif sum(ismember(art,counter))==1
            currtitre = horzcat('Component ',num2str(counter),':Artifact');
            fontcol = 'r';
        end
        
        axes('pos',ax1(counter).Position)
        topoplot(EEG.icawinv(:,counter),EEG.chanlocs)
        
        set(f1,'CurrentAxes',ax1(counter));
        set(gcf,'Color',[1 1 1]);
        
        ax1(counter).Title.String = currtitre;
        ax1(counter).Title.Color = fontcol;
        set(ax1(counter),'XColor','w','YColor','w')
        set(ax1(counter),'HitTest','on','SelectionHighlight','on','UserData',{EEG, typecomp, counter, windhandle, horiz, vert, blink, disc,...
            soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD,...
            soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, max_var, soglia_D, maxdin, activations},'NextPlot','add');
        set(ax1(counter),'ButtonDownFcn',@showICinfo)
        
    end
    waitfor(f1,'Close')
    %% NOTE THE COMPONENTS THAT YOU WANT TO REJECT
    % Create an edit box in which to define the components to be rejected.
    % The indices of those ICs to reject is stored in the
    % EEG.reject.gcompreject structure.
    
    prompt={'Enter indices of components to reject: '};
    dlg_title='Reject Components';
    num_lignes=1;
    deflts={'1'};
    rejica_ui = inputdlg(prompt,dlg_title,num_lignes,deflts);
    ica2rej = str2double(rejica_ui{1,1});
    EEG.reject.gcompreject(ica2rej) = 1;
    
    fprintf(fid,'The following ICs were marked for rejection: %2d\n', ica2rej);
    fclose(fid)
    
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(ica_nom),'filepath',pathcurr);
    eeglab redraw
    
end

if isVis==1 % Compare data before and after IC rejection
    
   if sum(EEG.reject.gcompreject)>0
    sphw = cell(1,2);   % Initialise cell array of sphere and weight data.
    sphw{1,1} = EEG.icasphere;
    sphw{1,2} = EEG.icaweights;
    
    winv = EEG.icawinv; % inverse or pseudoinverse
    compon_ind = 1:size(EEG.icaweights,1);
    ic2rej = find(EEG.reject.gcompreject);
    comps2proj = compon_ind(~ismember(compon_ind,ic2rej)); 
    
    disp('-----Back projecting ICs retained-----------------------');
    [projs, pvaf] = compvar(EEG.data(1:icalen,:), sphw, winv, comps2proj);
    
    eegplot(EEG.data(1:icalen,:), 'srate',EEG.srate,'eloc_file',EEG.chanlocs(1:icalen),'dispchans',icalen,'title',...
        'Data with and without ICs removed','xgrid','off','ygrid','off','data2',projs); 
   elseif sum(EEG.reject.gcompreject)==0
       msg = 'No Component indices defined in EEG structure!';
       title = 'Oh no!!';
       f = warndlg(msg,title);
   end
    
end
    
if isrej ==1   %just carry out component rejection
    
    
    sphw = cell(1,2);   % Initialise cell array of sphere and weight data.
    sphw{1,1} = EEG.icasphere;
    sphw{1,2} = EEG.icaweights;
    
    winv = EEG.icawinv; % inverse or pseudoinverse
    compon_ind = 1:size(EEG.icaweights,1);
    
    % Create an edit box in which to define the components to be rejected.
    if isVis==1
        waitfor(gcf,'Close')
        prompt={'Enter indices of components to reject: '};
        dlg_title='Reject Components';
        num_lignes=1;
        deflts={num2str(find(EEG.reject.gcompreject))};   % automatically displays those components already marked for rejection.
        rejica_ui = inputdlg(prompt,dlg_title,num_lignes,deflts);
        ica2rej = str2double(rejica_ui{1,1});
    else
        prompt={'Enter indices of components to reject: '};
        dlg_title='Reject Components';
        num_lignes=1;
        deflts={num2str(find(EEG.reject.gcompreject))};   % automatically displays those components already marked for rejection.
        rejica_ui = inputdlg(prompt,dlg_title,num_lignes,deflts);
        ica2rej = str2double(rejica_ui{1,1});
    end
    
    fprintf(fid,'The following ICs were rejected: %2d\n', ica2rej);
    fclose(fid)
    
    
    comps2proj = compon_ind(~ismember(compon_ind,ica2rej));  % Find the indices of the components to retain and project.
    disp('-----Back projecting ICs retained-----------------------');
    [projs, pvaf] = compvar(EEG.data(1:icalen,:), sphw, winv, comps2proj);  % Call of function to back project components.
    
    disp(horzcat('Those components back-projected account for ',num2str(pvaf),'% of the variance'))
    
    %% PLOT THE DATA WITH AND WITHOUT SELECTED ICA COMPONENTS REMOVED
    
    eegplot(EEG.data(1:icalen,:), 'srate',EEG.srate,'eloc_file',EEG.chanlocs(1:icalen),'dispchans',icalen,'title',...
        'Data with and without ICs removed','plottitle','blue=before IC rej red=after IC rej','xgrid','off','ygrid','off','data2',projs); 
    
    EEG.data = projs;
    icarej_nom = strcat(EEG.setname,'-icarej');
    
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(icarej_nom),'gui','off'); % current set = xx;
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(icarej_nom),'filepath',pathcurr);
    eeglab redraw
    
end % end if isRej if loop

end
