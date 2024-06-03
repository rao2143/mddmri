clear all

wd = cd;

% DataDir = '/opt/nmrdata/daniel';
%DataDir = '/Users/daniel/Dropbox/NMRdata/AVII500';
%DataDir = '/Users/daniel/Dropbox';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
% DataDir = '/Users/daniel/Dropbox/NMRdata/AV4_500/mariag'; ExpNam = {'210322_SC_batch'};
DataDir = '/Users/daniel/Dropbox/NMRdata/AV4_500/daniel';
cd(DataDir)

ExpNam = {'20230428_brainlipid_temp';'20230501_smalec_temp'};
% ExpNam = {'20230427_schweineschmalz_temp'};
% ExpNam = {'20230426_griebenschmalz_temp'};
% ExpNam = {'20230420_fakemeat_temp';'20230421_fakemeat_temp'};
% ExpNam = {'20230329_salamifryingfat_temp'};
% ExpNam = {'20230306_brainlipid_temp'};
%ExpNam = {'20230224_fakemeat';'20230223_coffeepowder_temp';'20230222_kidneybeanwet_temp';'20230208_dijonmustard_temp';'20230207_cremefraiche_temp'; '20230206_butter_temp'; '20230203_sesame_wet_temp'};
% ExpNam = {'KEMM57_VT23'};
% ExpNam = {'20230302_finnishfat_temp';'20230301_smarowidlo_temp';'20230228_smalec_temp';'20230227_bacon';'20230224_fakemeat';'20230224_dispulp_temp_cycle2';'20230224_dispulp_temp';'20230223_coffeepowder_temp';'20230222_kidneybeanwet_temp';'20230213_dispulp_temp';'20230210_dispulp_temp'};
% ExpNam = {'20230208_dijonmustard_temp'};
% ExpNam = {'20230207_cremefraiche_temp'; '20230206_butter_temp'; '20230205_CPMAS_Tjumptest'; '20230203_sesame_wet_temp'};
% ExpNam = {'20221221_SCwet'};
% ExpNam = {'20221125_SCdry_temp'};
% ExpNam = {'20220922_CTADNAwet_edumbosdhetcor'};
% ExpNam = {'20220803_CTADNA_temp2'};
% ExpNam = {'20220222_POPCGM3AS'};
% ExpNam = {'20220919_CTADNAwet_edumbosdhetcor';'20220907_CTADNAdry_edumbosdhetcor';'20220905_DPPC_edumbosdhetcor';'20220903_DPPCdry_temp';'20220901_DPPCdry_temp';'20220816_CTADNAdry_sdhetcor';'20220815_CTADNAdry_temp';'20220813_CTADNAdry_temp';'20220811_CTADNAdry_temp';'20220805_CTADNA';'20220803_CTADNA_temp2';'20220803_CTADNA_temp';'20220802_CTADNA_temp';'20220725_CTADNA';'20220706_CTADNA';'20220509_POPCGM3'};
% ExpNam = {'20220920_dispulp_temp'};
% ExpNam = {'20220912_dispulp_temp'};
% ExpNam = {'composite10mm_dor_temp_20220404'};
% ExpNam = {'composite10mm_trteaxde_dor_20220216'};
% ExpNam = {'20220223_POPCGM3AS'}; 
% ExpNam = {'20220110_POPCGM3'};
%ExpNam = {'MOACRARE2D_test'};
%ExpNam = {'MOACRARE2D_160213'};
%ExpNam = {'Sofia_FEXSY'};
%ExpNam = {'DTD2'};
%ExpNam = {'DTD'};
% ExpNam = {'AOToct_Eq7'};
% ExpNam = {'AOToct_Eq6'};
% ExpNam = {'AOToct_Eq5'};
% ExpNam = {'AOToct_Eq4'};
% ExpNam = {'AOToct_Eq3'};
%ExpNam = {'AOToct_temp11'};
%ExpNam = {'AOToct_Eq1'};
%ExpNam = {'AOToct_temp10'};
%ExpNam = {'Pstarch_temp1'};
%ExpNam = {'Pstarch_temp2'};
%ExpNam = {'AOToct_temp9'};
%ExpNam = {'AOToct_temp8'};
%ExpNam = {'AOToct_temp7'};
%ExpNam = {'AOToct_temp'};
%ExpNam = {'AOToct_temp3'};
%ExpNam = {'AOToct_temp'};
%ExpNam = {'AOTisooctane'};
%ExpNam = {'Lalpha_3Dimag'};
%ExpNam = {'MLV_3Dimag'};

for ndir = 1:length(ExpNam)
    ExpNam{ndir}
    cd(ExpNam{ndir})

    GetExpnos

    maxexpno = max(expno);
    maxexpno = 1+floor(maxexpno/10)*10;

    % expno = 11:10:maxexpno;
    te = zeros(size(expno));
    year = zeros(size(expno));
    month = zeros(size(expno));
    day = zeros(size(expno));
    hour = zeros(size(expno));
    minute = zeros(size(expno));
    second = zeros(size(expno));

    for nexp = 1:length(expno)
        %if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0
    %     if expno(nexp)>2461
    %         res = fExpnoInfo2('Y','N','N',expno(nexp));
    %     end
        eval(['load ' num2str(expno(nexp)) '/NMRacqus'])
        te(nexp) = NMRacqus.te;
        year(nexp) = NMRacqus.year;
        month(nexp) = NMRacqus.month;
        day(nexp) = NMRacqus.day;
        hour(nexp) = NMRacqus.hour;
        minute(nexp) = NMRacqus.minute;
        second(nexp) = NMRacqus.second;
    end
    
    time = second+60*(minute+60*(hour+24*(day+31*(month+12*year))));
    time = time-time(1);

    fs = 20; lw = 2;
    
    figure(1), clf
    axh1 = axes('position',[.2 .2 .75 .35],'FontSize',fs);
    ph1 = plot(expno,te,'-o');
    xlabel('expno'), ylabel('te')

    axh2 = axes('position',[.2 .6 .75 .35],'FontSize',fs);
    ph2 = plot(expno,time/60/60/24,'-o');
    ylabel('time / days')

    set([ph1 ph2],'LineWidth',lw)
    set([axh1 axh2],'LineWidth',lw,'TickDir','out','Box','off')
    print -dpdf TempvsExpno
    
    cd ..
end

cd(wd)