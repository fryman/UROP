%cmd Run dos
% lfp_read  ////load event file
% lfp_add// add files to analuse
% lfp_SpikeNames
% lfp_Spikes
% help lfp_spikeAnalysis
%  lfp_ras([],1, [-1 1])
% lfp_showEvents(10)
% lfp_his([],1, [-0.2 1]) 
% lfp_read2(varargin)
% lfp_disp(1:50, [], [-3 5],'avg','err2' )
% lfp_selectByRule('HasEvent(Stim_VTA)')
% lfp_AlignmentRef=Stim_VTA
% lfp_disp([], [], [-.5 2],'avg','err2' )
% myevtidx = ismember(lfp_Events(:,2), [8 16 32 64]);
%diff(lfp_Events(myevtidx,1))
%diff(lfp_Events(ismember(lfp_Events(:,2), [8 16 32]),1))
%mystarts = ismember(lfp_Events(:,2), [8 16 32]);
%lfp_createEvents(@(a) a, 100, lfp_Events(mystarts,1) + 1.9)
% find(lfp_SelectedTrials)
plot (x,y,'.', 'MarkerSize',2)
[ts x y a targ pts hdr] = dg_Nlx2MatVT('C:\Users\Alexander\Desktop\Rat4\2011-04-29_22-24-04_rat4_tr_17_45\vt1.nvt');
%C:\Users\Alexander\Desktop\Tradeoff project\Matlab Analuses>bulkprocess fiflefirst.m BulcPreprocesing logfile
bulkprocess list.m leif_stimProcessing log_leif1 
bulkprocess My_Temp_file_list.m AF_BulcPreprocesing logfile
bulkprocess ListPDF.m AF_triger_plotPDFAligment2 log
bulkprocess fiflefirst.m  logfile
 bulkprocess  listRat1_2_4_all.m AF_funCal logfileAugust165
 
 bulkprocess list2.m AF_trialTimes LogFileSept
 bulkprocess list2.m AF_funCal LogFileAug26 AF_trialTimes
 AF_trialTimesChandMixSeporate
 
 bulkprocess list2.m AF_trialTimesChandMixSeporate LogFileSept16
TrakerCliner (sessiondir)
%Dan sugest use 
matlab -nojvm
path(pathdef)
 TS = Nlx2MatTT(filename,1,0,0,0,0,1);
lfp_EventColors{556}='r'


hold on
plot (get(gca,'XLim'),258 * [1 1])


 x2 = dg_superdeglitch(x, 5, 1);
    y2 = dg_superdeglitch(y, 5, 1);
    x3 = dg_superdeglitch(x2, 5, 5);
    y3 = dg_superdeglitch(y2, 5, 5);
    x5 = dg_superdeglitch(x3, 15, 5);
    y5 = dg_superdeglitch(y3, 15, 5);



% \\chunky.mit.edu\smbshare\software\lgibb\Striosome_Project_Spring_2011
dbstop if error
help debug
dbstack
dbquit
» whos

HELP MEMORY
dg_concatsessions(dir1, dir2, dir3)
lfp_saveValue(varname, varargin)
bulkprocess list.m TrakerCliner2 listlog
bulkprocess list2.m AF_funCal listlog
Lfp spec
LFp disp 
lftp spike to wave 
lfp wave smosss 
other analuses 
lfp_saveValue(varname, varargin)
Plan for today look to Effect on stimulation on 2 daqys on 3 rats 
look what happens befor start lic and after look to start of trial else 