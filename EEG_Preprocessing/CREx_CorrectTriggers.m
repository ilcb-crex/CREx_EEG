function [ALLEEG, EEG,CURRENT]= CREx_CorrectTriggers(ALLEEG,EEG,CURRENTSET)

% Programmed by: Deirdre Bolger
% Function to correct trigger codes that are in an incorrect format when
% imported into EEGLAB.
% It saves a new EEGLAB dataset with the extension '*_trigcorrect'.
% Use as follows in Matlab command window:
% [ALLEEG, EEG,CURRENT]= CREx_CorrectTriggers(ALLEEG,EEG,CURRENTSET)
%************************************************************************


%% CHECK IF THE TRIGGERS CODES ARE IN STRING FORM

if ischar(EEG.event(2).type)==1
    display('**************************************Triggers are in string form*******************************************')
    
    for ctrig1=1:length(EEG.event)
        EEG.event(ctrig1).type=str2double(EEG.event(ctrig1).type);
    end
end

%% CHECK THAT THE TRIGGERS ARE IN DOUBLE FORMAT AND CORRECT IF NECESSARY

Y=[EEG.event.type];
Y=Y';
if numel(num2str(Y(3)))>3      %implies not in correct form need to correct the triggers
    Trig_corr=TrigFix(Y);
    
    for ctrig2=1:length(Trig_corr)
        EEG.event(ctrig2).type=Trig_corr(ctrig2);
    end
    
    [ALLEEG,EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    
else   %end of numel if
    display('Triggers corrected');
end

newtitle=strcat(EEG.setname,'_trigcorrect');
EEG = pop_saveset( EEG, 'filename',newtitle,'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

%**************************Call of subfunction to correct tiggers codes

    function Events_Correct=TrigFix(Event_Data)
        
        Events_Correct=zeros(length(Event_Data),1);
        
        for counter =1:length(Event_Data)
            if isnan(Event_Data(counter))==1
                continue;
            else
                x=num2str(dec2bin(Event_Data(counter)));
                Events_Correct(counter,1)=bin2dec(x(9:16));
            end
        end
    end

end
