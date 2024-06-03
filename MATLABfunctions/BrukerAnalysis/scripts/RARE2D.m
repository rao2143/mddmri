clear allwd = cd;% DataDir = '/Users/daniel/Dropbox/NMRdata/Darmstadt';%DataDir = '/opt/topspin2/data/DT/nmr';%DataDir = '/Users/daniel/NMRdata/AVII500/DT';%DataDir = '/Users/daniel/Dropbox';%DataDir = '/Users/daniel/Documents/Spaces/Presentations';%DataDir = '/opt/data/DT/nmr/';%DataDir = '/Users/daniel/Dropbox/NMRdata/SNC/HD3#600';%DataDir = '/opt/nmrdata/daniel';DataDir = '/Users/daniel/Dropbox/NMRdata/AV4_500/daniel';ExpNam = {'20230419_mousebrain'};% ExpNam = {'MIC5_10mm_setup'}; expno = [172 182 222 232 242 252 262 272 282 292 294:296 298 302:303 312 314 322 323 332 333 342 343 352 353 362 363 372 373]; expno = 373;% ExpNam = {'mousebrain_20220530'}; expno = [3 5 7 13 22]; expno = 22;%ExpNam = {'mousebrain_20220523'}; expno = [3]; expno = 3;%ExpNam = {'mrfood_trteaxde_dor_20220218'}; expno = [3 13 23 32 42:43 52:53]; expno = 53;%ExpNam = {'composite10mm_trteaxde_dor_20220216'}; expno = [3 9 13 17 23]; %expno = 3;%ExpNam = {'popc_trteaxde_dor_20220214'}; expno = [2:4 7:8 11:13 22 37]; expno = 37;%ExpNam = {'mousetumor_trteaxde_dor_20220106'}; %expno = [6]; expno = 163;%ExpNam = {'MouseTumors_20170619'; 'LiquidCrystals_20160920'; 'LiquidCrystals_20160422'; 'GradientTest'; 'DTrare2d_setup'};% ExpNam = {'axderare2dms_test'}; expno = 8;%ExpNam = {'qdor_test'};%ExpNam = {'composite10mm_temp_20210928'}; expno = [2:10:102]; expno = 102;%ExpNam = {'trteaxde_composite'};% ExpNam = {'DTrare2d_setup'}; expno = [6];%ExpNam = {'ToothPickPhantom'}; expno = [2];% ExpNam = {'DTI_HHnerve'};%ExpNam = {'MIC5_10mm_setup'};%ExpNam = {'uFA_RARE'}; expno = 100;% ExpNam = {'isoaniso_test'};%ExpNam = {'DiffVACSY'}; expno = 92;%ExpNam = {'qVAS_Starch'}; expno = 55;%ExpNam = {'AOToct'}; expno = 117;%ExpNam = {'AOToct_temp4'}; expno = [3 9 13 19 23 29 33 39];%ExpNam = {'AOToct_temp5'};    %ExpNam = {'AOToct_temp6'}; expno = [13:10:263 17:10:273];%ExpNam = {'Yeast_tripgse'}; expno = [78];%ExpNam = {'C14E5'}; expno = [66 71 78 83 84:159 165 166:189 195]; expno = 200:220; expno = 344;%ExpNam = {'C14E5_2'}; expno = 1;%ExpNam = {'C14E5_3'}; expno = [2:3 16 18 20 22]; expno = 22;%ExpNam = {'MIC5_2H_setup'}; expno = 21; expno = 4; expno =5;%ExpNam = {'SjolundOpt'}; expno = 6;%ExpNam = {'MIC5_2H_setup'}; expno = [14 15]; expno = 5;%ExpNam = {'MIC5_31P_setup'}; expno = 4;%ExpNam = {'C14E5_5'};%ExpNam = {'DDeltaMap'};%ExpNam = {'AOToct8'};%ExpNam = {'AOToct9'};%ExpNam = {'C14E5_6'};%ExpNam = {'C14E5_7'};%ExpNam = {'Saupe_SliceTest'};%ExpNam = {'AOToct_Eq1'};%ExpNam = {'AOToct_Eq2'};%ExpNam = {'AOToct_Eq4'};%ExpNam = {'AOToct_Eq6'};%ExpNam = {'AOToct_Eq7'};%ExpNam = {'MIC5test'}; expno = 53;%ExpNam = {'DTrare2d_setup'}; expno = 2;%ExpNam = {'LiquidCrystals_20160422'}; expno = 23;%ExpNam = {'LiquidCrystals_20160920'};% lb = 150e-6; si = 2*64; si1 = 2*32; %for z-y image% lb = 75e-6; si = 128; si1 = si; %for x-y image%lb = 150e-6; si = 64; si1 = si; %for x-y image%lb = 300e-6; si = 16; si1 = si; %for x-y image%DataDir = '/Users/daniel/NMRdata/AVII500/';%ExpNam = {'Norah'};fontsize = 15;cd(DataDir)%GetExpnamsfor ndir = 1:length(ExpNam)    ExpNam{ndir}    cd(ExpNam{ndir})    %clear expno    if exist('expno') == 0        GetExpnos    end    for nexp = 1:length(expno)        ConvertAcqus = 'Y';    ConvertProcs = 'N';    MakeTextfile = 'N';        if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0            ConvertAcqus = 'Y';        end        res = fExpnoInfo2(ConvertAcqus,MakeTextfile,ConvertProcs,expno(nexp));        eval(['load ' num2str(expno(nexp)) '/NMRacqus'])        if all([any(strcmp(NMRacqus.pulprog,{'DT_rare2d','DT_rare2d_AVIII'})); ...                exist([num2str(expno(nexp)) '/fid'])])%             try                sRARE2D                Imax = abs(I);                Imax(nudim.i/2+(0:2),nudim.j/2+(0:2)) = 0;                Imax = max(reshape(Imax,numel(Imax),1));                figure(1), clf                                        axes('position',[.15 .15 .85 .85],'FontSize',fontsize)                Iplot = abs(I)/Imax; clim = [0 1]; colormap(gray(256));                %Iplot = real(I)/Imax; clim = [0 1]; colormap(hot(256));                %Iplot = angle(I); clim = pi*[-1 1]; colormap(jet(256));                imagesc(r.i*1e3,r.j*1e3,Iplot',clim)                set(gca,'YDir','normal')                axis equal, axis tight                xlabel('read / mm'), ylabel('phase / mm')                delete([num2str(expno(nexp)) '/ReportFig*'])                delete([num2str(expno(nexp)) '/acqu*sconv'])                delete([num2str(expno(nexp)) '/*Dat.mat'])                delete([num2str(expno(nexp)) '/*.jpg'])                set(gcf, 'PaperPosition', [0 0 11 8.5],'PaperSize', [11 8.5]);                 eval(['print ' num2str(expno(nexp)) '/ReportFig -loose -dpdf'])                %NMRacqus.te%             catch%                 warning(['Error for ' ExpNam{ndir} ' expno ' num2str(expno(nexp)) ])%             end        end    end    cd ..endcd(wd)