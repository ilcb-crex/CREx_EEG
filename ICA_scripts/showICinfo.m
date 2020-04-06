function showICinfo(hdl,~)

D=get(hdl,'UserDat');
EEG = D{1,1};
typecomp = D{1,2};
numpcomp = D{1,3};
windhandle = D{1,4};
soglia_DV = D{1,9};
diff_var = D{1,10};
soglia_K = D{1,11};
med2_K = D{1,12};
meanK = D{1,13};
soglia_SED = D{1,14};
med2_SED = D{1,15};
SED = D{1,16};
soglia_SAD = D{1,17};
med2_SAD = D{1,18};
SAD = D{1,19};
soglia_GDSF = D{1,20};
med2_GDSF = D{1,21};
GDSF = D{1,22};
soglia_V = D{1,23};
med2_V = D{1,24};
max_var = D{1,25};
soglia_D =D{1,26};
maxdin = D{1,27};
activations = D{1,28};


is_horiz = sum(ismember(D{1,5},numpcomp));
is_vert = sum(ismember(D{1,6},numpcomp));
is_blink = sum(ismember(D{1,7},numpcomp));
is_disc = sum(ismember(D{1,8},numpcomp));

com = pop_prop_ADJ_crex(EEG, typecomp,numpcomp , windhandle, is_horiz, is_vert, is_blink, is_disc,...
     soglia_DV, diff_var(1,numpcomp), soglia_K, med2_K,  meanK(1,numpcomp), soglia_SED, med2_SED, SED(1,numpcomp), soglia_SAD, med2_SAD,...
     SAD(1,numpcomp), soglia_GDSF, med2_GDSF, GDSF(1,numpcomp), soglia_V, med2_V, max_var(1,numpcomp), soglia_D, maxdin(1,numpcomp), activations);

end
