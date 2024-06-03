clear all

wd = cd;

% cd(['/opt/nmrdata/daniel/'])
cd(['/Users/daniel/Dropbox/NMRdata/AV4_500/daniel/'])

ExpNam = {'20221220_SCdry'};
expno = [14]; PT_type = 'DPCPINEPT';
ydisp = .5; yscale = 3;

% cd('/Users/daniel/Dropbox/NMRdata/Faellanden')
% ExpNam = {'ZHAA03_dTopgaard_210322_SC_MariaG_20221117'};
% expno = [100]; PT_type = 'DPCPINEPT';
% ydisp = .5; yscale = 30;

% cd('/Users/daniel/Dropbox/NMRdata/AV4_500/mariag'), ExpNam = {'210322_SC_batch'};
% expno = [15]; PT_type = 'DPCPINEPT';
% ydisp = .5; yscale = 3;

fs = 2*5;
lw = 2*.5;

nexp = length(expno);



bottom_cm = 1.5;
dheight_cm = 1;
papersize_cm = 2*[8.3 bottom_cm+dheight_cm*nexp];
if nexp == 1
    papersize_cm = 2*[8.3 bottom_cm+dheight_cm*5];
end
bottom = bottom_cm/papersize_cm(2);

for ndir = 1:length(ExpNam)
    ExpNam{ndir}
    cd(ExpNam{ndir})

    figure(1), clf
    axes('position',[0 bottom 1 1-bottom],'LineWidth',lw,'FontSize',fs)

    load([num2str(expno(1)) '/NMRacqus.mat']);
    for nexp = length(expno):-1:1
        load([num2str(expno(nexp)) '/NMRacqus.mat']);
        NMRacqus.te
        
        if strcmp(PT_type,'DPCPINEPT')
            DP = load([num2str(expno(nexp)) '/SpecDat.mat']);
            CP = load([num2str(expno(nexp)+1) '/SpecDat.mat']);
            INEPT = load([num2str(expno(nexp)+2) '/SpecDat.mat']);

            Imax = 2*1.1*max([DP.Spec.I; CP.Spec.I; INEPT.Spec.I]);

            figure(1)
            h = plot(DP.Spec.ppm,DP.Spec.I/Imax + (nexp-1),'k-');
            set(h,'LineWidth',1.5*lw,'Color',.8*[1 1 1])
            hold on

            h = plot(INEPT.Spec.ppm,INEPT.Spec.I/Imax + (nexp-1),'k-');
            set(h,'LineWidth',1*lw,'Color',[1 0 0])

            h = plot(CP.Spec.ppm,CP.Spec.I/Imax + (nexp-1),'k-');
            set(h,'LineWidth',.5*lw,'Color',[0 0 1])
            
            h = plot(DP.Spec.ppm,yscale*DP.Spec.I/Imax + ydisp + (nexp-1),'k-');
            set(h,'LineWidth',1.5*lw,'Color',.8*[1 1 1])
            h = plot(INEPT.Spec.ppm,yscale*INEPT.Spec.I/Imax + ydisp + (nexp-1),'k-');
            set(h,'LineWidth',1*lw,'Color',[1 0 0])
            h = plot(CP.Spec.ppm,yscale*CP.Spec.I/Imax + ydisp + (nexp-1),'k-');
            set(h,'LineWidth',.5*lw,'Color',[0 0 1])

            delete([num2str(expno(nexp)+0) '/ReportFig.*'])
            delete([num2str(expno(nexp)+1) '/ReportFig.*'])
            delete([num2str(expno(nexp)+2) '/ReportFig.*'])
        elseif strcmp(PT_type,'DPCP')
            DP = load([num2str(expno(nexp)) '/SpecDat.mat']);
            CP = load([num2str(expno(nexp)+1) '/SpecDat.mat']);

            Imax = 1.1*max([DP.Spec.I; CP.Spec.I]);

            figure(1)
            h = plot(DP.Spec.ppm,DP.Spec.I/Imax + (nexp-1),'k-');
            set(h,'LineWidth',1.5*lw,'Color',.8*[1 1 1])
            hold on

            h = plot(CP.Spec.ppm,CP.Spec.I/Imax + (nexp-1),'k-');
            set(h,'LineWidth',.5*lw,'Color',[0 0 1])
            
            delete([num2str(expno(nexp)+0) '/ReportFig.*'])
            delete([num2str(expno(nexp)+1) '/ReportFig.*'])
        elseif strcmp(PT_type,'DP')
            DP = load([num2str(expno(nexp)) '/SpecDat.mat']);

            Imax = 1.1*max([DP.Spec.I]);

            figure(1)
            h = plot(DP.Spec.ppm,DP.Spec.I/Imax + (nexp-1),'k-');
            set(h,'LineWidth',1*lw,'Color',0*[1 1 1])
            hold on
            
            delete([num2str(expno(nexp)+0) '/ReportFig.*'])
        end

        if exist('ppmmin') == 0
            ppmmin = min(DP.Spec.ppm);
            ppmmax = max(DP.Spec.ppm);
        end
%         TextString = [num2str(NMRacqus.te,3) ' K'];
%         text(ppmmax-1,nexp-1+.15,TextString,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fs)
        TextString = ['\times' num2str(yscale)];
        text(ppmmax-1,.5+nexp-1+.15,TextString,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fs)
        
    end
    figure(1)
    hold off
    set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,...
        'Box','off','TickDir','out','Ycolor','none')
    Imin = -.05;
    Imax = length(expno)+ .01;
    axis([ppmmin ppmmax Imin Imax])
    xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'])
    
    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize_cm],'PaperSize', papersize_cm);
    set(gcf,'Render','painters')
%     set(gcf,'color','none','inverthardcopy','off');

    eval(['print -dpdf -loose ' num2str(expno(1)) '/ReportFig'])

    cd ..
end

cd(wd)
