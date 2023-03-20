function T= drow_table(delta,CK_RMSE,PF_RMSE,rr)
if strcmp(rr, 'delta1')
    g='delta';
elseif strcmp(rr, 'alpha1')
    g='alpha';
elseif strcmp(rr, 'alpha1T')
    g='Time'; 
else
    g='M';  
end
for i=1:length(delta)
ftable{1,i} =[g,' = (',num2str(delta(i)),')'];
com_rat(i)   =CK_RMSE(i);
rms(i)       =PF_RMSE(i);
end
%rr={'averge comunication ';'RMSE'}
T = table(ftable',com_rat',rms');
T.Properties.VariableNames={g 'DECKF' 'DEPF'};
end