clear all%cd('/Users/daniel/NMRdata/AVII500/Celine/')cd('/Users/daniel/Dropbox/nmrdata/AVII500/Celine')ExpNam = {'DMPC_Dat'}; expno = [5 7 8 9 10 15]; expno = 15;ppmlimits = [-60 60];Intensity = 1.5e6;T2 = 1;deltaiso = 2;deltaaniso = 40; %ppmeta = .1;baseline = -.01*Intensity;masr = 1250;% ExpNam = {'Coaggregate_DMPS_151007'}; expno = [11 15 19 27 31 35 43]; %expno = 31;% ppmlimits = [-60 -2 7 60];% Intensity = 10e5;% T2 = 0.5;% deltaiso = 2.5;% deltaaniso = 45; %ppm% eta = .01;% masr = 1250;% ExpNam = {'Pure_DMPS'}; expno = [4 8 12 31 35]; %expno = 4;% ppmlimits = [-60 -4 1 60];% Intensity = 1e7;% T2 = 1;% deltaiso = -1;% deltaaniso = 40; %ppm% eta = .001;% % expno = [23 27];% % ppmlimits = [-60 1 5 60];% % Intensity = 1e7;% % T2 = 2;% % deltaiso = 2.8;% % deltaaniso = 40; %ppm% % eta = .001;fs = 10;lw = 1;MkReportFig = 1;for ndir = 1:length(ExpNam)    ExpNam{ndir}    cd(ExpNam{ndir})    for nexp = 1:length(expno)        load([num2str(expno(nexp)) '/NMRacqus.mat']); NMRacqus.masr = masr;        load([num2str(expno(nexp)) '/SpecDat.mat'])        ppm = Spec.ppm;        I = real(Spec.I);        fitpoints = [];        seriesno = [];        for nppm = 1:2:length(ppmlimits)-1            ind = find([ppm>ppmlimits(nppm) & ppm<ppmlimits(nppm+1)]);            fitpoints = [fitpoints; ind];            seriesno = [seriesno; (nppm+1)/2*ones(size(ind))];        end        %figure(1), clf, plot(ppm,I,ppm(fitpoints),I(fitpoints),'r'), axis([min(ppm) max(ppm) max(I(fitpoints))*[-.1 1.1]]),     set(gca,'XDir','reverse','YTick',[],'LineWidth',lw), return        Xin = Spec.ppm(fitpoints);        Yin = I(fitpoints);        Para = [NMRacqus.bf1*1e6 NMRacqus.masr];        Pin = [Intensity T2 deltaiso deltaaniso eta baseline]; Funam = 'fCSAMAS_BL';        LB = [.1*Intensity .1*T2 deltaiso-.5 deltaaniso-20 0 -.5*Intensity];        %UB = [10*Intensity 10*T2 deltaiso+.5 deltaaniso+20 .01 .5*Intensity];        UB = [10*Intensity 10*T2 deltaiso+.5 deltaaniso+20 1 .5*Intensity];        Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = abs(Pin);         Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,LB./Pnorm,UB./Pnorm,[],Pnorm,Xnorm,Ynorm,Para);        %Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,[],[],[],Pnorm,Xnorm,Ynorm,Para);        Ycalc = feval(Funam,Pout,Xin,ones(size(Pin)),1,1,Para); error = Yin - Ycalc;        FitDat.Ycalc = Ycalc;        FitDat.Xcalc = Xin;        FitDat.residual = error;        chisq = sum(error.^2)/length(error)        %figure(1), clf, plot(Xin,Yin,Xin,Ycalc,Xin,error), cd .., return        FitDat.Intensity = Pout(1);        FitDat.T2 = Pout(2);        FitDat.deltaiso = Pout(3);        FitDat.deltaaniso = Pout(4);        FitDat.eta = Pout(5);        FitDat.ppm = ppm;        FitDat.I = I;        figure(1), clf        axes('position',[-.01 .18 1.02 .72])        plot(ppm,I,'k','LineWidth',2)        hold on        xlabel('\it\delta\rm(^3^1P) / ppm','FontSize',fs)        set(gca,'XDir','reverse','YTick',[],'LineWidth',lw)        axis([[min(ppm) max(ppm)]+.02*NMRacqus.sw_h/NMRacqus.bf1*[-1 1] max(max(max(FitDat.Ycalc)))*[-.2 1.5]]), grid        for nrange = 1:max(seriesno)            ind = nrange == seriesno;            plot(Xin(ind),FitDat.Ycalc(ind),'r','LineWidth',1)            plot(Xin(ind),error(ind)/max(abs(error))*max(FitDat.Ycalc)*.1 - max(FitDat.Ycalc)*.1,'k'), grid        end                hold off        set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,'Box','off','TickDir','out','Ycolor',[1 1 1])        title([num2str(expno(nexp)) ' \it\delta\rm_{iso} = ' num2str(FitDat.deltaiso,3) ...            ' ppm \it\delta\rm_{aniso} = ' num2str(FitDat.deltaaniso,3) ' ppm \it\eta\rm = ' num2str(FitDat.eta,3) ...            '  \itT\rm_2 = ' num2str(FitDat.T2,3) ' s'])        set(gca, 'FontSize',fs,'LineWidth',lw)        eval(['save ' num2str(expno(nexp)) '/DeconvDat FitDat'])                fname = [num2str(expno(nexp)) '/SpecExp.txt'];        DatMat = [Spec.ppm I];        save(fname, 'DatMat', '-ascii', '-tabs')        fname = [num2str(expno(nexp)) '/SpecFit.txt'];        DatMat = [Xin Yin Ycalc error seriesno];        save(fname, 'DatMat', '-ascii', '-tabs')                if MkReportFig            set(gcf, 'PaperPosition', [0 0 11 8.5],'PaperSize', [11 8.5]);             eval(['print ' num2str(expno(nexp)) '/ReportFig -loose -dpdf'])        end    endend