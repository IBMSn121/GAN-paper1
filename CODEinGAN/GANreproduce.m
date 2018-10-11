function GANreproduce
warning off; idxStart = 0; idxAnaTp = 0; inxEmp = 1;
test_str = '\fontsize{16}{16}\selectfont Starting of GANreproduce (version 1.4)?';
Select1 = ' Yes start now ';
Select2 = ' No thank you ';
options.Default = Select1;
options.Interpreter = 'latex';
Answer1 = questdlg(test_str, 'Optional List', Select1, Select2, options);
switch Answer1
    case Select1
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp(['We are starting now.'])
        idxStart = 1; %pause(1);
    case Select2
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
        idxStart = 0; idxAnaTp = 0; inxEmp = 1;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Starting of Analysis%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%%%1. Select of Analysis
%%%
while  idxStart == 1;
test_str1 = '\fontsize{16}{16}\selectfont Which analyses do you want to conduct?';
Select1 = ' Simulation ';
Select2 = ' Real Data ';
Select3 = ' Linear Trend ';
Select4 = ' Surv-BACH2 ';
Select5 = ' Surv-Trio ';
Select6 = ' No thank you and close the dialog ';
Select7a = ' Next... ';
Select7b = ' Next.... ';
Select7c = ' Select again ';

options.Default = Select7a; idxStart = 0; idxAnaTp = 0;
Answer2a = questdlg(test_str1, 'Optional List', Select1, Select2, Select7a, options);
switch Answer2a
    case Select1
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp(['We are now conducting Simulation Analysis.'])
        idxAnaTp = 1; idxStart = 0;
    case Select2
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp(['We are now conducting Real Data Analysis.'])
        idxAnaTp = 2; idxStart = 0;
    case Select7a
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp(['Next selections'])
        idxAnaTp = 991;
end

if idxAnaTp == 991
options.Default = Select7b; idxAnaTp = 0;
Answer2b = questdlg(test_str1, 'Optional List', Select3, Select4, Select7b, options);
switch Answer2b
    case Select3
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp(['We are now conducting Linear Trend Analysis.'])
        idxAnaTp = 3; idxStart = 0;
    case Select4
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp(['We are now conducting Surv-BACH2 Analysis.'])
        idxAnaTp = 4; idxStart = 0;
    case Select7b
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp(['Next selections'])
        idxAnaTp = 992;
end
end%if 991

if idxAnaTp == 992
options.Default = Select7c; idxAnaTp = 0;
Answer2c = questdlg(test_str1, 'Optional List', Select5, Select6, Select7c, options);
switch Answer2c
    case Select5
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp(['We are now conducting Surv-Trio Analysis.'])
        idxAnaTp = 5; idxStart = 0;
		
    case Select6
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
        idxStart = 0; idxAnaTp = 0; inxEmp = 1;
		
    case Select7c
        disp(['                    '])
        disp(['                    '])
	    disp(['                    '])
        disp(['Next selections'])
        idxStart = 1;
end
end%if 992
end; %while for idxStart



%%%
%%%
%%%Double Check Conditions
if idxAnaTp > 0
test_str1 = '\fontsize{16}{16}\selectfont DATAinGAN was downloaded and put in C:/ ? ';
Select1 = ' Yes ensured ';
Select2 = ' Not sure ';

options.Default = Select1; inxEmp = 1;
AnswerF1 = questdlg(test_str1, 'Check necessary files 1', Select1, Select2,options);
switch AnswerF1
    case Select1
	    disp(['                    '])
        disp(['Analysis Preparation 1: all data files in DATAinGAN have been checked for loading.'])
        disp(['                    '])
        inxEmp = 0;
    case Select2
        disp(['                    '])
        disp(['Please ensure and try again'])
        disp(['                    '])
        idxAnaTp = 0; inxEmp = 1;
end
end

if idxAnaTp*(inxEmp == 0 ) > 0
test_str1 = '\fontsize{16}{16}\selectfont All GAN m-files in Current Folder of Matlab? ';
Select1 = ' Yes ensured ';
Select2 = ' Not sure ';

options.Default = Select1; inxEmp = 1;
AnswerF1 = questdlg(test_str1, 'Check necessary files 2', Select1, Select2,options);
switch AnswerF1
    case Select1
	    disp(['                    '])
        disp(['Analysis Preparation 2: all GAN m-files in CODEinGAN have been checked for loading.'])
        disp(['                    '])
        inxEmp = 0;
    case Select2
        disp(['                    '])
        disp(['Please ensure and try again'])
        disp(['                    '])
        idxAnaTp = 0; inxEmp = 1;
end
end
%%%for Double Check
%%%
%%%




%%%
%%%Analysis of Simu; simuAna.m
%%%idxAnaTp=1;
idxAna1=0;
if idxAnaTp*(inxEmp == 0 )==1;
Select11 = ' Methods comparisons ';
Select12 = ' Robust AWTE in large samples ';
test_str1 = '\fontsize{16}{16}\selectfont Selecting simulation types:';
options.Default = Select11;
options.Interpreter = 'latex';
Answer1 = questdlg(test_str1, 'Optional List', Select11, Select12,options);
switch Answer1
    case Select11
	    disp(['                    '])
        disp(['''Methods comparisons'' selected.'])
		disp(['                    '])
        idxAna1=1;
    case Select12
	    disp(['                    '])
        disp(['''Robust AWTE in large samples'' selected.'])
		disp(['                    '])
        idxAna1=2;	
end
end%if for pre-check Simu types 

%Simu type 1
if idxAna1 == 1;


chkAna=(idxAnaTp == 1)*1; idxNoiseL=0; 
while  chkAna == 1;
opts.Resize = 'on'; 
opts.Interpreter = 'latex';
prompt = '\fontsize{16}{16}\selectfont Please input Noise Level: 1 (lower), 2 (median), 3 (higher), 4 (worst) or 0 (quit)';
title = 'Select Noise Level';
dims = [1 120];
definput = {'2'};
AnswerNL = inputdlg(prompt,title,dims,definput,opts); idxNoiseL=str2double(AnswerNL);
if (idxNoiseL ~= 1 & idxNoiseL ~= 2 & idxNoiseL ~= 3 & idxNoiseL ~= 4 & idxNoiseL ~= 0); chkAna=1; disp(['Not available input, please select again']); end;
if (idxNoiseL == 1 | idxNoiseL == 2 | idxNoiseL == 3 | idxNoiseL == 4); chkAna=0; disp(['Your selected Noise Level = ' num2str(idxNoiseL) ', i.e., ' ' Sigma2 = ' num2str(idxNoiseL/2)]); end;  
if (idxNoiseL == 0); chkAna=0; disp('Have a nice day. Goodbye!'); end; inxEmp=double(isempty(AnswerNL));
if (inxEmp == 1); chkAna=0; disp('Have a nice day. Goodbye!'); end; 
end; %while for chkAna

chkAna=(idxAnaTp == 1)*(1-inxEmp)*(idxNoiseL > 0); idxSimuRn=0;
while  chkAna == 1;
opts.Resize = 'on'; 
opts.Interpreter = 'latex';
prompt = '\fontsize{16}{16}\selectfont Please input Simulation Runs: 1 (Median, taking 30 mins), 2 (Full, taking 300 mins) or 0 (quit)';
title = 'Select Simulation Runs';
dims = [1 135];
definput = {'2'};
AnswerSR = inputdlg(prompt,title,dims,definput,opts); idxSimuRn=str2double(AnswerSR); SimuRn=100*(idxSimuRn==1)+1000*(idxSimuRn==2);
if (idxSimuRn ~= 1 & idxSimuRn ~= 2 & idxSimuRn ~= 0); chkAna=1; disp(['Not available input, please select again']); end;
if (idxSimuRn == 2); chkAna=0; disp(['You selected Full Runs, i.e., Simulation Runs = ' num2str(SimuRn)] ); end; 
if (idxSimuRn == 1); chkAna=0; disp(['You selected Median Runs, i.e., Simulation Runs = ' num2str(SimuRn)] ); end;  
if (idxSimuRn == 0); chkAna=0; disp('Have a nice day. Goodbye!'); end; inxEmp=double(isempty(AnswerSR));
if (inxEmp == 1); chkAna=0; disp('Have a nice day. Goodbye!'); end; 
end; %while for chkAna

if idxAnaTp*(1-inxEmp)*(idxSimuRn > 0)*(idxNoiseL > 0)==1; simuAna(idxNoiseL/2,SimuRn); end;


end; %idxAna1 = 1 for Simu type 1


%Simu type 2
if idxAna1 == 2;


chkAna=(idxAnaTp == 1)*(1-inxEmp);
while  chkAna == 1;
opts.Resize = 'on'; 
opts.Interpreter = 'latex';
prompt = '\fontsize{16}{16}\selectfont Please input Simulation Runs: 1 (Median, taking 3 mins), 2 (Full, taking 30 mins) or 0 (quit)';
title = 'Select Simulation Runs';
dims = [1 135];
definput = {'2'};
AnswerSR = inputdlg(prompt,title,dims,definput,opts); idxSimuRn=str2double(AnswerSR); SimuRn=100*(idxSimuRn==1)+1000*(idxSimuRn==2);
if (idxSimuRn ~= 1 & idxSimuRn ~= 2 & idxSimuRn ~= 0); chkAna=1; disp(['Not available input, please select again']); end;
if (idxSimuRn == 2); chkAna=0; disp(['You selected Full Runs, i.e., Simulation Runs = ' num2str(SimuRn)] ); end; 
if (idxSimuRn == 1); chkAna=0; disp(['You selected Median Runs, i.e., Simulation Runs = ' num2str(SimuRn)] ); end;  
if (idxSimuRn == 0); chkAna=0; disp('Have a nice day. Goodbye!'); end; inxEmp=double(isempty(AnswerSR));
if (inxEmp == 1); chkAna=0; disp('Have a nice day. Goodbye!'); end; 
end; %while for chkAna

if idxAnaTp*(1-inxEmp)*(idxSimuRn > 0)==1; simuAwAna(SimuRn); end;


end; %idxAna1 = 2 for Simu type 2




%%%
%%%Analysis of Real; realAna.m
%%%idxAnaTp=2; inxEmp=0;

idxRD1=1; chkAna=(idxAnaTp == 2)*(1-inxEmp)*(idxRD1 > 0); idxRD2=0;
while  chkAna == 1;
opts.Resize = 'on'; 
opts.Interpreter = 'latex';
prompt = '\fontsize{16}{16}\selectfont Two datasets realDS10gn and realUS5gn have been in C:/DATAinGAN? 1 (Yes) or 0 (quit)';
title = 'Check Dataset Existence';
dims = [1 135];
definput = {'1'};
AnswerRD2 = inputdlg(prompt,title,dims,definput,opts); idxRD2=str2double(AnswerRD2); inxEmp=double(isempty(AnswerRD2));
if (idxRD2 ~= 1 & idxRD2 ~= 0); chkAna=1; disp(['Not available input, please check again']); end; 
if (idxRD2 == 1); chkAna=0; disp(['You have checked realDS10gn and realUS5gn are in C:\DATAinGAN'] ); end; 
if (idxRD2 == 0); chkAna=0; disp('Have a nice day. Goodbye!'); end;
if (inxEmp == 1); chkAna=0; disp('Have a nice day. Goodbye!'); end;   
end; %while for chkAna


if idxRD1*idxRD2*idxAnaTp==2; dXX=load('C:\DATAinGAN\realUS5gn.txt'); dYY=load('C:\DATAinGAN\realDS10gn.txt'); realAna(dXX,dYY); end;
%???realAnav2 for summary of Fig S1 (i.e., coeff&Pv) OK
%???
%???realAnav2.附表修正.NA
%???
%???realAnav2.本文修正.NA
%???
%???Quick guideline.pdf for GANreproduce




%%%
%%%Analysis of LT3X2; LT3X2Ana.m
%%%idxAnaTp=3; inxEmp=0;
idxLT1=1; chkAna=(idxAnaTp == 3)*(1-inxEmp)*(idxLT1 > 0); idxLT2=0;
while  chkAna == 1;
opts.Resize = 'on';
opts.Interpreter = 'latex';
prompt = '\fontsize{16}{16}\selectfont Ten datasets \{ ltaDS10gnD''s, ltaUS5gnD''s \} have been in C:/DATAinGAN? 1 (Yes) or 0 (quit)';
title = 'Check Dataset Existence';
dims = [1 150];
definput = {'1'};
AnswerLT2 = inputdlg(prompt,title,dims,definput,opts); idxLT2=str2double(AnswerLT2); inxEmp=double(isempty(AnswerLT2));
if (idxLT2 ~= 1 & idxLT2 ~= 0); chkAna=1; disp(['Not available input, please check again']); end; 
if (idxLT2 == 1); chkAna=0; disp(['You have checked ltaDS10gnDs and ltaUS5gnDs are in C:\DATAinGAN'] ); end; 
if (idxLT2 == 0); chkAna=0; disp('Have a nice day. Goodbye!'); end;
if (inxEmp == 1); chkAna=0; disp('Have a nice day. Goodbye!'); end;   
end; %while for chkAna


dXX1=load('C:\DATAinGAN\ltaUS5gnD1.txt'); dYY1=load('C:\DATAinGAN\ltaDS10gnD1.txt');
dXX2=load('C:\DATAinGAN\ltaUS5gnD2.txt'); dYY2=load('C:\DATAinGAN\ltaDS10gnD2.txt');
dXX3=load('C:\DATAinGAN\ltaUS5gnD3.txt'); dYY3=load('C:\DATAinGAN\ltaDS10gnD3.txt');
dXX4=load('C:\DATAinGAN\ltaUS5gnD4.txt'); dYY4=load('C:\DATAinGAN\ltaDS10gnD4.txt');
dXX60=load('C:\DATAinGAN\ltaUS5gnD60.txt'); dYY60=load('C:\DATAinGAN\ltaDS10gnD60.txt');
if idxLT1*idxLT2*idxAnaTp==3; LT3X2Ana(dXX1,dYY1,dXX2,dYY2,dXX3,dYY3,dXX4,dYY4,dXX60,dYY60); end;




%%%
%%%Analysis of SurvB; survBAna1.m survBAna2.m survBAna3.m
%%%idxAnaTp=4; inxEmp=0;
chkAna=(idxAnaTp == 4)*(1-inxEmp); idxSB1=0;
while  chkAna == 1;
opts.Resize = 'on';
opts.Interpreter = 'latex';
prompt = '\fontsize{16}{16}\selectfont Four DLBCL datasets \{sDLBCL''s\} have been in C:/DATAinGAN? 1 (Yes) or 0 (quit)';
title = 'Check Dataset Existence';
dims = [1 150];
definput = {'1'};
AnswerSB1 = inputdlg(prompt,title,dims,definput,opts); idxSB1=str2double(AnswerSB1); inxEmp=double(isempty(AnswerSB1));
if (idxSB1 ~= 1 & idxSB1 ~= 0); chkAna=1; disp(['Not available input, please check again']); end; 
if (idxSB1 == 1); chkAna=0; disp(['You have checked 4 DLBCL datasets are in C:\DATAinGAN'] ); 
test_strA = 'Please note that Log-rank test was reckoned by package of Fan Lin: www.mathworks.com/matlabcentral/fileexchange/20388';
Select1 = ' OK ';
AnswerA = questdlg(test_strA, 'Info Note', Select1, Select1);
disp(['Note that Log-rank test was reckoned by package of Fan Lin'] );
end; 
if (idxSB1 == 0); chkAna=0; disp('Have a nice day. Goodbye!'); end;
if (inxEmp == 1); chkAna=0; disp('Have a nice day. Goodbye!'); end;   
end; %while for chkAna


%%%Conduct 3-types of analyses for Clinical controversy
chkAna=idxSB1; idxSB2=0; inxEmpa=0;
while  chkAna == 1;
opts.Resize = 'on';
opts.Interpreter = 'latex';
prompt = '\fontsize{16}{16}\selectfont Now perform: 1 (Controversy), 2 (BACH2\(^{+}\) better), 3 (BACH2\(^{+}\) worse) or 0 (quit)?';
title = 'Check Analysis Types';
dims = [1 150];
definput = {'1'};
AnswerSB2 = inputdlg(prompt,title,dims,definput,opts); idxSB2=str2double(AnswerSB2); inxEmpa=double(isempty(AnswerSB2));
if (idxSB2 ~= 1 & idxSB2 ~= 2 & idxSB2 ~= 3 & idxSB2 ~= 0); chkAna=1; disp(['Not available input, please check again']); end; 

if (idxSB2 == 1);
disp(['Now we perform Clinical Controversy for BACH2''s role']);

dDLBCL1=load('C:\DATAinGAN\sDLBCL1.txt'); dDLBCL2=load('C:\DATAinGAN\sDLBCL2.txt');
dDLBCL3=load('C:\DATAinGAN\sDLBCL3.txt'); dDLBCL4=load('C:\DATAinGAN\sDLBCL4.txt');

survBAna1(dDLBCL1,dDLBCL2,dDLBCL3,dDLBCL4); pause(10);
test_strB = '\fontsize{16}{16}\selectfont Continue other Surv-BACH2 Analyses? ';
options.Default = Select1; chkAna=0;
options.Interpreter = 'latex';
Select1 = ' Yes ';
Select2 = ' No ';
AnswerB = questdlg(test_strB, 'Info Note', Select1, Select2, options);
switch AnswerB
    case Select1
	    disp(['                    '])
        disp(['Surv-BACH2 Analysis Continued.'])
		disp(['                    '])
		chkAna=1;
    case Select2
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
		disp(['                    '])
        chkAna=0;
end
end; %if idxSB2 == 1

if (idxSB2 == 2);
disp(['Now we perform BACH2 higher and survival better in patients with lower SPIB expression']);

dDLBCL1=load('C:\DATAinGAN\sDLBCL1.txt'); dDLBCL2=load('C:\DATAinGAN\sDLBCL2.txt');
dDLBCL3=load('C:\DATAinGAN\sDLBCL3.txt'); dDLBCL4=load('C:\DATAinGAN\sDLBCL4.txt');

survBAna2(dDLBCL1,dDLBCL2,dDLBCL3,dDLBCL4); pause(8);
test_strB = '\fontsize{16}{16}\selectfont Continue other Surv-BACH2 Analyses? ';
options.Default = Select1; chkAna=0;
options.Interpreter = 'latex';
Select1 = ' Yes ';
Select2 = ' No ';
AnswerB = questdlg(test_strB, 'Info Note', Select1, Select2, options);
switch AnswerB
    case Select1
	    disp(['                    '])
        disp(['Surv-BACH2 Analysis Continued.'])
		disp(['                    '])
		chkAna=1;
    case Select2
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
		disp(['                    '])
        chkAna=0;
end
end; %if idxSB2 == 2

if (idxSB2 == 3);
disp(['Now we perform BACH2 higher and survival worse in patients with higher SPIB expression']);

dDLBCL1=load('C:\DATAinGAN\sDLBCL1.txt'); dDLBCL2=load('C:\DATAinGAN\sDLBCL2.txt');
dDLBCL3=load('C:\DATAinGAN\sDLBCL3.txt'); dDLBCL4=load('C:\DATAinGAN\sDLBCL4.txt');

survBAna3(dDLBCL1,dDLBCL2,dDLBCL3,dDLBCL4); pause(8);
test_strB = '\fontsize{16}{16}\selectfont Continue other Surv-BACH2 Analyses? ';
options.Default = Select1; chkAna=0;
options.Interpreter = 'latex';
Select1 = ' Yes ';
Select2 = ' No ';
AnswerB = questdlg(test_strB, 'Info Note', Select1, Select2, options);
switch AnswerB
    case Select1
	    disp(['                    '])
        disp(['Surv-BACH2 Analysis Continued.'])
		disp(['                    '])
		chkAna=1;
    case Select2
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
		disp(['                    '])
        chkAna=0;
end
end; %if idxSB2 == 3


if (idxSB2 == 0); chkAna=0; disp('Have a nice day. Goodbye!'); end;
if (inxEmpa == 1); chkAna=0; disp('Have a nice day. Goodbye!'); end;   
end; %while for chkAna





%%%
%%%Analysis of SurvT; 
%%%idxAnaTp=5; inxEmp=0;
chkAna=(idxAnaTp == 5)*(1-inxEmp); idxLT1=0;
while  chkAna == 1;
opts.Resize = 'on';
opts.Interpreter = 'latex';
prompt = '\fontsize{16}{16}\selectfont All 13 blood cancer datasets \{ sDLBCL''s, sOSBC''s \} have been in C:/DATAinGAN? 1 (Yes) or 0 (quit)';
title = 'Check Dataset Existence';
dims = [1 150];
definput = {'1'};
AnswerST1 = inputdlg(prompt,title,dims,definput,opts); idxLT1=str2double(AnswerST1); inxEmp=double(isempty(AnswerST1));
if (idxLT1 ~= 1 & idxLT1 ~= 0); chkAna=1; disp(['Not available input, please check again']); end; 
if (idxLT1 == 1); chkAna=0; disp(['You have checked all 13 blood cancer datasets are in C:\DATAinGAN'] ); end; 
if (idxLT1 == 0); chkAna=0; disp('Have a nice day. Goodbye!'); end;
if (inxEmp == 1); chkAna=0; disp('Have a nice day. Goodbye!'); end;   
end; %while for chkAna

%sDLBCLs         (survTAna1.m)  DLBCL
%sOSBC2a sOSBC2b (survTAna2.m)  FL
%sOSBC3a sOSBC3b (survTAna3.m)  MM
%sOSBC4          (survTAna4.m)  CLL
%sOSBC5          (survTAna5.m)  MCL
%sOSBC6a sOSBC6b (survTAna6.m)  AML
%sOSBC7          (survTAna7.m)  ALL
%%%Conduct multiple blood cancer prognosis
chkAna=idxLT1; idxLT2=0; inxEmpa=0;
while  chkAna == 1;
opts.Resize = 'on';
opts.Interpreter = 'latex';
prompt = '\fontsize{16}{16}\selectfont Now perform: 1 (DLBCL), 2 (FL), 3 (MM), 4 (CLL), 5 (MCL), 6 (AML), 7 (ALL) or 0 (quit)?';
title = 'Check Analysis Types';
dims = [1 150];
definput = {'1'};
AnswerSB2 = inputdlg(prompt,title,dims,definput,opts); idxLT2=str2double(AnswerSB2); inxEmpa=double(isempty(AnswerSB2));
if (idxLT2 ~= 1 & idxLT2 ~= 2 & idxLT2 ~= 3 & idxLT2 ~= 4 & idxLT2 ~= 5 & idxLT2 ~= 6 & idxLT2 ~= 7 & idxLT2 ~= 0); chkAna=1; disp(['Not available input, please check again']); end; 



if (idxLT2 == 1);
disp(['Now we perform Trio''s prognosis for DLBCL']);

dDLBCL1=load('C:\DATAinGAN\sDLBCL1.txt'); dDLBCL2=load('C:\DATAinGAN\sDLBCL2.txt');
dDLBCL3=load('C:\DATAinGAN\sDLBCL3.txt'); dDLBCL4=load('C:\DATAinGAN\sDLBCL4.txt');
survTAna1(dDLBCL1,dDLBCL2,dDLBCL3,dDLBCL4); pause(8);

test_strB = '\fontsize{16}{16}\selectfont Continue other Surv-Trio Analyses? ';
options.Default = Select1; chkAna=0;
options.Interpreter = 'latex'; 
Select1 = ' Yes ';
Select2 = ' No ';
AnswerB = questdlg(test_strB, 'Info Note', Select1, Select2, options);
switch AnswerB
    case Select1
	    disp(['                    '])
        disp(['Surv-Trio Analysis Continued.'])
		disp(['                    '])
		chkAna=1;
    case Select2
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
		disp(['                    '])
        chkAna=0;
end
end; %if idxLT2 == 1



if (idxLT2 == 2);
disp(['Now we perform Trio''s prognosis for FL']);

dFL1=load('C:\DATAinGAN\sOSBC2a.txt'); dFL2=load('C:\DATAinGAN\sOSBC2b.txt');
survTAna2(dFL1,dFL2); pause(8);

test_strB = '\fontsize{16}{16}\selectfont Continue other Surv-Trio Analyses? ';
options.Default = Select1; chkAna=0;
options.Interpreter = 'latex'; 
Select1 = ' Yes ';
Select2 = ' No ';
AnswerB = questdlg(test_strB, 'Info Note', Select1, Select2, options);
switch AnswerB
    case Select1
	    disp(['                    '])
        disp(['Surv-Trio Analysis Continued.'])
		disp(['                    '])
		chkAna=1;
    case Select2
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
		disp(['                    '])
        chkAna=0;
end
end; %if idxLT2 == 2



if (idxLT2 == 3);
disp(['Now we perform Trio''s prognosis for MM']);

dMM1=load('C:\DATAinGAN\sOSBC3a.txt'); dMM2=load('C:\DATAinGAN\sOSBC3b.txt');
survTAna3(dMM1,dMM2); pause(8);

test_strB = '\fontsize{16}{16}\selectfont Continue other Surv-Trio Analyses? ';
options.Default = Select1; chkAna=0;
options.Interpreter = 'latex'; 
Select1 = ' Yes ';
Select2 = ' No ';
AnswerB = questdlg(test_strB, 'Info Note', Select1, Select2, options);
switch AnswerB
    case Select1
	    disp(['                    '])
        disp(['Surv-Trio Analysis Continued.'])
		disp(['                    '])
		chkAna=1;
    case Select2
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
		disp(['                    '])
        chkAna=0;
end
end; %if idxLT2 == 3



if (idxLT2 == 4);
disp(['Now we perform Trio''s prognosis for CLL']);

dCLL=load('C:\DATAinGAN\sOSBC4.txt'); 
survTAna4(dCLL); pause(8);

test_strB = '\fontsize{16}{16}\selectfont Continue other Surv-Trio Analyses? ';
options.Default = Select1; chkAna=0;
options.Interpreter = 'latex'; 
Select1 = ' Yes ';
Select2 = ' No ';
AnswerB = questdlg(test_strB, 'Info Note', Select1, Select2, options);
switch AnswerB
    case Select1
	    disp(['                    '])
        disp(['Surv-Trio Analysis Continued.'])
		disp(['                    '])
		chkAna=1;
    case Select2
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
		disp(['                    '])
        chkAna=0;
end
end; %if idxLT2 == 4



if (idxLT2 == 5);
disp(['Now we perform Trio''s prognosis for MCL']);

dMCL=load('C:\DATAinGAN\sOSBC5.txt'); 
survTAna5(dMCL); pause(8);

test_strB = '\fontsize{16}{16}\selectfont Continue other Surv-Trio Analyses? ';
options.Default = Select1; chkAna=0;
options.Interpreter = 'latex'; 
Select1 = ' Yes ';
Select2 = ' No ';
AnswerB = questdlg(test_strB, 'Info Note', Select1, Select2, options);
switch AnswerB
    case Select1
	    disp(['                    '])
        disp(['Surv-Trio Analysis Continued.'])
		disp(['                    '])
		chkAna=1;
    case Select2
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
		disp(['                    '])
        chkAna=0;
end
end; %if idxLT2 == 5



if (idxLT2 == 6);
disp(['Now we perform Trio''s prognosis for AML']);

dAML1=load('C:\DATAinGAN\sOSBC6a.txt'); dAML2=load('C:\DATAinGAN\sOSBC6b.txt');
survTAna6(dAML1,dAML2); pause(8);

test_strB = '\fontsize{16}{16}\selectfont Continue other Surv-Trio Analyses? ';
options.Default = Select1; chkAna=0;
options.Interpreter = 'latex'; 
Select1 = ' Yes ';
Select2 = ' No ';
AnswerB = questdlg(test_strB, 'Info Note', Select1, Select2, options);
switch AnswerB
    case Select1
	    disp(['                    '])
        disp(['Surv-Trio Analysis Continued.'])
		disp(['                    '])
		chkAna=1;
    case Select2
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
		disp(['                    '])
        chkAna=0;
end
end; %if idxLT2 == 6



if (idxLT2 == 7);
disp(['Now we perform Trio''s prognosis for ALL']);

dALL=load('C:\DATAinGAN\sOSBC7.txt'); 
survTAna7(dALL); pause(8);

test_strB = '\fontsize{16}{16}\selectfont Continue other Surv-Trio Analyses? ';
options.Default = Select1; chkAna=0;
options.Interpreter = 'latex'; 
Select1 = ' Yes ';
Select2 = ' No ';
AnswerB = questdlg(test_strB, 'Info Note', Select1, Select2, options);
switch AnswerB
    case Select1
	    disp(['                    '])
        disp(['Surv-Trio Analysis Continued.'])
		disp(['                    '])
		chkAna=1;
    case Select2
	    disp(['                    '])
        disp('Have a nice day. Goodbye!')
		disp(['                    '])
        chkAna=0;
end
end; %if idxLT2 == 7


if (idxLT2 == 0); chkAna=0; disp('Have a nice day. Goodbye!'); end;
if (inxEmpa == 1); chkAna=0; disp('Have a nice day. Goodbye!'); end;   
end; %while for chkAna





