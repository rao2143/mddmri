clear all

wd = cd;

% cd(['/opt/nmrdata/daniel/'])
% cd(['/Users/daniel/Dropbox/NMRdata/AV4_500/daniel/'])
% cd('/Users/daniel/Dropbox/NMRdata/Faellanden')
DataDir = '/Users/daniel/Dropbox/NMRdata/AV4_500/daniel';

% ExpNam = {'20230428_brainlipid_temp';'20230501_smalec_temp'};
% expno = [12:10:192]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; %Imax = 5e10;
% % % expno = [22 102 182]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; Imax = .2e10;
% % expno = 4+[11:10:191]; PT_type = 'DP'; ppmmin = -1; ppmmax = 11; Imax = .05e10;

% ExpNam = {'20230428_brainlipid_temp'};
% expno = 0+[11:10:191]; PT_type = 'DP'; ppmmin = -11; ppmmax = 21; Imax = 3e10;
% 
% ExpNam = {'20230427_schweineschmalz_temp'};
% expno = [12:10:192]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; %Imax = 5e10;
% % expno = [22 102 182]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; Imax = .2e10;
% % expno = 0+[11:10:191]; PT_type = 'DP'; ppmmin = -1; ppmmax = 11; Imax = 5e10;

ExpNam = {'20230426_griebenschmalz_temp'};
expno = [12:10:192]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; %Imax = 5e10;
% expno = [22 102 182]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; Imax = .2e10;
% expno = [11:10:191]; PT_type = 'DP'; ppmmin = -100; ppmmax = 111; Imax = .05e10;
% expno = 4+[11:10:191]; PT_type = 'DP'; ppmmin = -1; ppmmax = 11; Imax = 5e10;

% ExpNam = {'20230420_fakemeat_temp';'20230421_fakemeat_temp'};
% % ExpNam = {'20230420_fakemeat_temp'};
% expno = [12:10:192]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% expno = [11:10:191]; PT_type = 'DP'; ppmmin = -1; ppmmax = 11; Imax = .2e10;

% ExpNam = {'20230420_fakemeat'};
% expno = [2]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [1:10:11]; PT_type = 'DP';

% ExpNam = {'20230329_salamifryingfat_temp'};
% %expno = [12:10:172]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% expno = [11:10:171]; PT_type = 'DP';

% ExpNam = {'20230306_brainlipid_temp'};
% expno = [12:10:132]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [11:10:131]; PT_type = 'DP';

% ExpNam = {'20230302_finnishfat_temp'};
% expno = [12:10:172]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [11:10:171]; PT_type = 'DP';
% 
% ExpNam = {'20230301_smarowidlo_temp'};
% expno = [12:10:192]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [11:10:191]; PT_type = 'DP';

% ExpNam = {'20230228_smalec_temp'};
% expno = [12:10:182]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [11:10:181]; PT_type = 'DP';

% ExpNam = {'KEMM57_VT23'};
% % expno = [12:10:32 45 55 62 72 82 95 105 112]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% expno = [105 55 32 22  95 12 45 82 112 72 62]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [11:10:111]; PT_type = 'DP';

% ExpNam = {'20230224_fakemeat'};
% expno = [2]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [12 2]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [1:10:11]; PT_type = 'DP';

% ExpNam = {'20230223_coffeepowder_temp'};
% expno = [12:10:132]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [11:10:131]; PT_type = 'DP';
% % 
% ExpNam = {'20230222_kidneybeanwet_temp'};
% expno = [12:10:132]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [11:10:141]; PT_type = 'DP';

% ExpNam = {'20230208_dijonmustard_temp'};
% expno = [12:10:162]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [11:10:161]; PT_type = 'DP';

% ExpNam = {'20230207_cremefraiche_temp'};
% expno = [12:10:212]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [11:10:211]; PT_type = 'DP';

% ExpNam = {'20230206_butter_temp'};
% expno = [22 52 152]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [12:10:202]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % % expno = [11:10:201]; PT_type = 'DP';

% ExpNam = {'20230203_sesame_wet_temp'}; 
% expno = [12:10:92]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205;
% % expno = [11:10:91]; PT_type = 'DP';

% ExpNam = {'20221125_SCdry_temp'};
% expno = [14:10:124 144]; PT_type = 'DPCPINEPT';
% expno = [12:10:122 142]; PT_type = 'DPCP';
% expno = [11:10:121 141]; PT_type = 'DP';

% ExpNam = {'ZHAA03_dTopgaard_210322_SC_MariaG_20221117'};
% expno = [100]; PT_type = 'DPCPINEPT';

% ExpNam = {'20220110_POPCGM3'};
% expno = [7 23 175 209 293 325]; PT_type = 'DPCPINEPT';
% expno = [65 115 119 154 204 236 272 304]; PT_type = 'DPCP';

% ExpNam = {'20220222_POPCGM3AS'};
% expno = [7 19 23 26]; PT_type = 'DPCPINEPT';

% ExpNam = {'20220223_POPCGM3AS'};
% expno = [12 61 98 137]; PT_type = 'DPCPINEPT';
% % expno = [9 58 95 134]; PT_type = 'DPCP';
% % expno = [11 60 97 136]; PT_type = 'DP';

% ExpNam = {'20220509_POPCGM3'};
% expno = [106 124 136 148]; PT_type = 'DPCPINEPT';
% % expno = -1+[106 124 136 148]; PT_type = 'DP';
% expno = -2+[106 124 136 148]; PT_type = 'DP';

% ExpNam = {'20220706_CTADNA'};
% expno = [11]; PT_type = 'DPCPINEPT';

% ExpNam = {'20220725_CTADNA'};
% expno = [9]; PT_type = 'DPCPINEPT';
% % expno = [12]; PT_type = 'DPCP';

% ExpNam = {'20220802_CTADNA_temp'};
% expno = 62; PT_type = 'DPCPINEPT';
% expno = 65; PT_type = 'DPCP';

% ExpNam = {'20220803_CTADNA_temp'};
% expno = 14:10:114; PT_type = 'DPCPINEPT';
% % expno = 12:10:112; PT_type = 'DPCP';
% % expno = 11:10:101; PT_type = 'DP';

% ExpNam = {'20220803_CTADNA_temp2'};
% expno = 14:10:104; PT_type = 'DPCPINEPT';
% % expno = 12:10:102; PT_type = 'DPCP';
% % expno = 11:10:101; PT_type = 'DP';

% ExpNam = {'20220805_CTADNA'};
% expno = 14:10:44; PT_type = 'DPCPINEPT';
% % expno = 12:10:42; PT_type = 'DPCP';
% % expno = 11:10:41; PT_type = 'DP';

% ExpNam = {'20220811_CTADNAdry_temp'};
% expno = 14:10:124; PT_type = 'DPCPINEPT';
% expno = 12:10:122; PT_type = 'DPCP';
% expno = 11:10:121; PT_type = 'DP';

% ExpNam = {'20220813_CTADNAdry_temp'};
% expno = 14:10:124; PT_type = 'DPCPINEPT';
% % expno = 12:10:122; PT_type = 'DPCP';
% % expno = 11:10:121; PT_type = 'DP';

% ExpNam = {'20220815_CTADNAdry_temp'};
% expno = 14:10:94; PT_type = 'DPCPINEPT';
% expno = 12:10:92; PT_type = 'DPCP';
% expno = 11:10:91; PT_type = 'DP';

% ExpNam = {'20220901_DPPCdry_temp'};
% expno = 14:10:104; PT_type = 'DPCPINEPT';
% % expno = 12:10:102; PT_type = 'DPCP';
% % expno = 11:10:101; PT_type = 'DP';

% ExpNam = {'20220903_DPPCdry_temp'};
% expno = 14:10:204; PT_type = 'DPCPINEPT';
% expno = 12:10:202; PT_type = 'DPCP';
% expno = 11:10:201; PT_type = 'DP';

% ExpNam = {'20220912_dispulp_temp'};
% expno = 12:10:212; PT_type = 'DPCPINEPT';
% % expno = 11:10:211; PT_type = 'DP';

% ExpNam = {'20220920_dispulp_temp'};
% expno = 12:10:62; PT_type = 'DPCPINEPT';
% % expno = 11:10:61; PT_type = 'DP';

% ExpNam = {'20220803_CTADNA_temp2'};
% expno = 14:10:104; PT_type = 'DPCPINEPT';
% % expno = 12:10:102; PT_type = 'DPCP';
% % expno = 11:10:101; PT_type = 'DP';

% ExpNam = {'20220222_POPCGM3AS'};
% expno = 19;
% ppmmin = 0; ppmmax = 220;

% ExpNam = {'C10E4relax'};
% expno = 22;
% expno = 84;
% 
% 
% cd('/Users/daniel/NMRdata/AVII500/Nils_C8G2/')
% ExpNam = {'C8G2_RH0_3','C8G2_RH11_3','C8G2_RH23_3','C8G2_RH33_3','C8G2_RH43_3','C8G2_RH68_3','C8G2_RH75_3','C8G2_RH96_3'};
% ExpNam = {'C8G2_RH96_3'};
% expno = 15:10:85;
% expno = 25;
% %ppmmin = -5; ppmmax = 700;
% 
% DataDir = '/Users/daniel/NMRdata/AVII500/DT';
% cd(DataDir)
% ExpNam = {'C8G2RH97'};
% expno = 54;
% 
% 
% DataDir = '/Users/daniel/NMRdata/AVII500/DT';
% cd(DataDir)
% ExpNam = {'Pstarch_temp2'};
% expno = [25 65 95];
% %expno = [35:10:85];
% ppmmin = 55; ppmmax = 109;
% 
% DataDir = '/Users/daniel/NMRdata/AVII500/DT';
% cd(DataDir)
% ExpNam = {'CelluloseDissolution'};
% expno = [15];
% ppmmin = 51; ppmmax = 119;
% 
% DataDir = '/Users/daniel/NMRdata/AVII500/DT';
% cd(DataDir)
% ExpNam = {'Pstarch_temp1'};
% expno = [35:10:125];
% ppmmin = 55; ppmmax = 109;

DataDir = '/Users/daniel/Dropbox/NMRdata/AV4_500/jenny';
ExpNam = {'yeastlipids'};
expno = [5:10:35]; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; %Imax = 5e10;
expno = [4:10:34]; PT_type = 'DP'; ppmmin = -10; ppmmax = 20; %Imax = 5e10;

cd(DataDir);
% GetExpnams

fs = 2*5;
lw = 2*.5;

nexp = length(expno);

bottom_cm = 2*0.75;
dheight_cm = 2*1;
papersize_cm = [2*8.3 bottom_cm+dheight_cm*nexp];
if nexp == 1
    papersize_cm = [2*8.3 bottom_cm+dheight_cm*3];
% elseif nexp == 2
%     papersize_cm = 2*[8.3 bottom_cm+dheight_cm*4];
end
bottom = bottom_cm/papersize_cm(2);

minI = -.1;
maxI = 1.1;
dI = maxI - minI;

for ndir = 1:length(ExpNam)
    ExpNam{ndir}
    cd(ExpNam{ndir})

    figure(1), clf
    axes('position',[0 bottom 1 1-bottom],'LineWidth',lw,'FontSize',fs)

    load([num2str(expno(1)) '/NMRacqus.mat']);
    for nexp = length(expno):-1:1
        load([num2str(expno(nexp)) '/NMRacqus.mat']);
        NMRacqus.te
        
        try
        if strcmp(PT_type,'DPCPINEPT')
            DP = load([num2str(expno(nexp)) '/SpecDat.mat']);
            CP = load([num2str(expno(nexp)+1) '/SpecDat.mat']);
            INEPT = load([num2str(expno(nexp)+2) '/SpecDat.mat']);

            if exist('Imax') == 1
                Imax_temp = Imax;
            else
                Imax_temp = 1*max([DP.Spec.I; CP.Spec.I; INEPT.Spec.I]);
            end

            figure(1)
            h = plot(DP.Spec.ppm,DP.Spec.I/Imax_temp/dI + (nexp-1),'k-');
            set(h,'LineWidth',1.5*lw,'Color',.8*[1 1 1])
            hold on

            h = plot(INEPT.Spec.ppm,INEPT.Spec.I/Imax_temp/dI + (nexp-1),'k-');
            set(h,'LineWidth',1*lw,'Color',[1 0 0])

            h = plot(CP.Spec.ppm,CP.Spec.I/Imax_temp/dI + (nexp-1),'k-');
            set(h,'LineWidth',.5*lw,'Color',[0 0 1])
            
            delete([num2str(expno(nexp)+0) '/ReportFig.*'])
            delete([num2str(expno(nexp)+1) '/ReportFig.*'])
            delete([num2str(expno(nexp)+2) '/ReportFig.*'])
        elseif strcmp(PT_type,'DPCP')
            DP = load([num2str(expno(nexp)) '/SpecDat.mat']);
            CP = load([num2str(expno(nexp)+1) '/SpecDat.mat']);

            if exist('Imax') == 1
                Imax_temp = Imax;
            else
                Imax_temp = 1*max([DP.Spec.I; CP.Spec.I]);
            end

            figure(1)
            h = plot(DP.Spec.ppm,DP.Spec.I/Imax_temp/dI + (nexp-1),'k-');
            set(h,'LineWidth',1.5*lw,'Color',.8*[1 1 1])
            hold on

            h = plot(CP.Spec.ppm,CP.Spec.I/Imax_temp/dI + (nexp-1),'k-');
            set(h,'LineWidth',.5*lw,'Color',[0 0 1])
            
            delete([num2str(expno(nexp)+0) '/ReportFig.*'])
            delete([num2str(expno(nexp)+1) '/ReportFig.*'])
        elseif strcmp(PT_type,'DP')
            DP = load([num2str(expno(nexp)) '/SpecDat.mat']);

            if exist('Imax') == 1
                Imax_temp = Imax;
            else
                Imax_temp = 1*max([DP.Spec.I]);
            end

            figure(1)
            h = plot(DP.Spec.ppm,DP.Spec.I/Imax_temp/dI + (nexp-1),'k-');
            set(h,'LineWidth',1*lw,'Color',0*[1 1 1])
            hold on
            
            delete([num2str(expno(nexp)+0) '/ReportFig.*'])
        end
        if exist('ppmmin') == 0
            ppmmin = min(DP.Spec.ppm);
            ppmmax = max(DP.Spec.ppm);
        end
        
        TextString = ['(' num2str(expno(nexp)) ') ' num2str(NMRacqus.te,3) ' K'];
        text(ppmmax-1,nexp-1+.15,TextString,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fs)

        catch
            warning(['no data found in expno ' num2str(expno(nexp))])
        end

        
    end
    TextString = ExpNam{ndir};
    text(ppmmax-1,length(expno)+minI,TextString,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',fs)
    
    figure(1)
    hold off
    position = get(gca,'Position');
    height_cm = position(4)*papersize_cm(2);
    width_cm = position(3)*papersize_cm(1);
    ticklength_cm = .2;
    ticklength = min([ticklength_cm/height_cm ticklength_cm/width_cm]);
    set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,...
        'Box','off','TickDir','out','TickLength',[ticklength 0.01],'Ycolor','none')
    axis([ppmmin ppmmax 0+minI length(expno)+minI])
    xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'])
    
    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize_cm],'PaperSize', papersize_cm);
    set(gcf,'Render','painters')
%     set(gcf,'color','none','inverthardcopy','off');

    eval(['print -dpdf -loose ' num2str(expno(1)) '/ReportFig'])

    cd ..
end

cd(wd)
