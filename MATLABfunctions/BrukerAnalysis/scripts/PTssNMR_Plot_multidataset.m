clear all

wd = cd;

% cd(['/opt/nmrdata/daniel/'])
% cd(['/Users/daniel/Dropbox/NMRdata/AV4_500/daniel/'])
% cd('/Users/daniel/Dropbox/NMRdata/Faellanden')
DataDir = '/Users/daniel/Dropbox/NMRdata/AV4_500/daniel';

cd(DataDir);
% GetExpnams

% ExpNam = {'20230224_fakemeat';'20230420_fakemeat'};
% expno_c = {[2];[2]}; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; Imax = 8e9;
% % expno_c = {[1];[1]};  PT_type = 'DP'; ppmmin = -1; ppmmax = 11; Imax = 2.5e10;

ExpNam = {'20230427_schweineschmalz_temp';'20230426_griebenschmalz_temp';'20230228_smalec_temp';'20230501_smalec_temp';'20230301_smarowidlo_temp';'20230329_salamifryingfat_temp'};
% expno_c = {[22];[22];[22];[22];[22];[22]}; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; %Imax = 8e9;
% expno_c = {[52];[52];[52];[52];[52];[52]}; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; %Imax = 8e9;
% expno_c = {[162];[162];[162];[162];[162];[162]}; PT_type = 'DPCPINEPT'; ppmmin = -5; ppmmax = 205; %Imax = 8e9;
expno_c = {[152];[152];[152];[152];[152];[152]}; PT_type = 'DP'; ppmmin = -5; ppmmax = 205; %Imax = 8e9;
% expno_c = {[1];[1]};  PT_type = 'DP'; ppmmin = -1; ppmmax = 11; Imax = 2.5e10;

fs = 2*5;
lw = 2*.5;

nexp = 0;
for ndir = 1:length(ExpNam)
    expno = expno_c{ndir};
    nexp = nexp + length(expno);
end

% bottom_cm = 2*0.75;
% dheight_cm = 2*1;
% papersize_cm = [2*8.3 bottom_cm+dheight_cm*nexp];
% if nexp == 1
%     papersize_cm = [2*8.3 bottom_cm+dheight_cm*3];
% elseif nexp == 2
%     papersize_cm = 2*[8.3 bottom_cm+dheight_cm*4];
% end
% bottom = bottom_cm/papersize_cm(2);

bottom_cm = 1.5;
papersize_cm = [20 10];
bottom = bottom_cm/papersize_cm(2);

minI = -.1;
maxI = 1.1;
dI = maxI - minI;

figure(1), clf
axes('position',[0 bottom 1 1-bottom],'LineWidth',lw,'FontSize',fs)
nexp_cum = 0;
for ndir = 1:length(ExpNam)
    ExpNam{ndir}
    expno = expno_c{ndir};
    cd(ExpNam{ndir})


    load([num2str(expno(1)) '/NMRacqus.mat']);
    for nexp = length(expno):-1:1
        load([num2str(expno(nexp)) '/NMRacqus.mat']);
        NMRacqus.te
        nexp_cum = nexp_cum+1;
        
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
            h = plot(DP.Spec.ppm,DP.Spec.I/Imax_temp/dI + (nexp_cum-1),'k-');
            set(h,'LineWidth',1.5*lw,'Color',.8*[1 1 1])
            hold on

            h = plot(INEPT.Spec.ppm,INEPT.Spec.I/Imax_temp/dI + (nexp_cum-1),'k-');
            set(h,'LineWidth',1*lw,'Color',[1 0 0])

            h = plot(CP.Spec.ppm,CP.Spec.I/Imax_temp/dI + (nexp_cum-1),'k-');
            set(h,'LineWidth',.5*lw,'Color',[0 0 1])
            
        elseif strcmp(PT_type,'DPCP')
            DP = load([num2str(expno(nexp)) '/SpecDat.mat']);
            CP = load([num2str(expno(nexp)+1) '/SpecDat.mat']);

            if exist('Imax') == 1
                Imax_temp = Imax;
            else
                Imax_temp = 1*max([DP.Spec.I; CP.Spec.I]);
            end

            figure(1)
            h = plot(DP.Spec.ppm,DP.Spec.I/Imax_temp/dI + (nexp_cum-1),'k-');
            set(h,'LineWidth',1.5*lw,'Color',.8*[1 1 1])
            hold on

            h = plot(CP.Spec.ppm,CP.Spec.I/Imax_temp/dI + (nexp_cum-1),'k-');
            set(h,'LineWidth',.5*lw,'Color',[0 0 1])
            
        elseif strcmp(PT_type,'DP')
            DP = load([num2str(expno(nexp)) '/SpecDat.mat']);

            if exist('Imax') == 1
                Imax_temp = Imax;
            else
                Imax_temp = 1*max([DP.Spec.I]);
            end

            figure(1)
            h = plot(DP.Spec.ppm,DP.Spec.I/Imax_temp/dI + (nexp_cum-1),'k-');
            set(h,'LineWidth',1*lw,'Color',0*[1 1 1])
            hold on
            
        end
        if exist('ppmmin') == 0
            ppmmin = min(DP.Spec.ppm);
            ppmmax = max(DP.Spec.ppm);
        end
        
        TextString = [ExpNam{ndir} ' (' num2str(expno(nexp)) ') ' num2str(NMRacqus.te,3) ' K'];
        text(ppmmax-1,nexp_cum-1+.15,TextString,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fs)

        catch
            warning(['no data found in expno ' num2str(expno(nexp))])
        end

        
    end
%     TextString = ExpNam{ndir};
%     text(ppmmax-1,length(expno)+minI,TextString,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',fs)
    

    cd ..
end

figure(1)
hold off
position = get(gca,'Position');
height_cm = position(4)*papersize_cm(2);
width_cm = position(3)*papersize_cm(1);
ticklength_cm = .2;
ticklength = min([ticklength_cm/height_cm ticklength_cm/width_cm]);
set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,...
    'Box','off','TickDir','out','TickLength',[ticklength 0.01],'Ycolor','none')
axis([ppmmin ppmmax 0+minI nexp_cum+minI])
xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'])

set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize_cm],'PaperSize', papersize_cm);
set(gcf,'Render','painters')
%     set(gcf,'color','none','inverthardcopy','off');

print(fullfile(DataDir,ExpNam{1},num2str(expno_c{1}(1)),'ReportFig'),'-dpdf','-loose')

figure(1)
hold off
position = get(gca,'Position');
height_cm = position(4)*papersize_cm(2);
width_cm = position(3)*papersize_cm(1);
ticklength_cm = .1;
ticklength = min([ticklength_cm/height_cm ticklength_cm/width_cm]);
set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,...
    'Box','off','TickDir','out','TickLength',[ticklength 0.01],'Ycolor','none')
axis([ppmmin ppmmax 0+minI nexp_cum+minI])
%     xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'])

set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize_cm],'PaperSize', papersize_cm);
set(gcf,'Render','painters')

print(fullfile(DataDir,ExpNam{1},num2str(expno_c{1}(1)),'ReportFig_ppt'),'-dpng','-loose','-r300')

cd(wd)
