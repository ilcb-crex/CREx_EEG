# EpochChannel_QualityCheck
Suite of functions allowing the user to check the quality of EEG data in an interactive manner; both continuous and epoched, in an interactive manner. It aids the user in deciding whether they need to exclude electrodes or take out individual epochs.It allows the user to make use of information from the frequency spectra (P-welch) of each electrode, a robust z-score calculation (Bigdeley-Shamlo et al, 2015) and the total amount of above-threshold time per electrode (for continuous data) and per trial (for segmented data) to help the user decide which electrodes/trials should be rejected.

To run, EEGLAB is required.
To run type "EpochChan_dlg(EEG) : the input variable, EEG, is the current dataset (structure) loaded into EEGLAB.
Note: If the program does not find any bad channels or data, according to the user-set criteria, it will stop suddenly. This will be rectified in the future.

This tool is still being constantly updated. 
