clear allwd = cd;%DataDir = '/Users/daniel/Dropbox';DataDir = '/Users/daniel/NMRdata/AVII500/DT/';%DataDir = '/Users/daniel/NMRdata/AVII200/';%DataDir = '/opt/topspin2/data/DT/nmr';%DataDir = '/opt/topspin/data/DT/nmr';cd(DataDir)ExpNam = {'Yeast_tripgse'};expno = [52:58 67 69 73 75 77]; Nb = 256;expno = 77; Nb = 256;AutoPhase = 1; si = 2*1024; lb = 20;CheckBasline = 0;FindPeaks = 1;CheckPeaks = 0;PlotInterm = 0;thresh = .15;td1start = 2;Imin = 2e-2;signal = 'area'; %signal = 'intensity';cd(ExpNam{1})if exist('expno') == 0    GetExpnosendfor nexp = 1:length(expno)    ConvertAcqus = 'Y';    ConvertProcs = 'N';    MakeTextfile = 'N';    if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0        ConvertAcqus = 'Y';    end    res = fExpnoInfo2(ConvertAcqus,MakeTextfile,ConvertProcs,expno(nexp));    eval(['load ' num2str(expno(nexp)) '/NMRacqus'])    if any(strcmp(NMRacqus.pulprog,{'DT_tbppgste','DT_tbppgsteT2'})) == 1        Spec2Dtd1%%        if any(strcmp(NMRacqus.pulprog,{'DT_tbppgste'})) == 1            Nvd = 1;        elseif any(strcmp(NMRacqus.pulprog,{'DT_tbppgsteT2'})) == 1            Nvd = NMRacqus.l1;        end        Ndir = td1/Nb/Nvd;        %         I_PA = reshape(Itd1,[si Nb Ndir Nvd]);%         I_PA = squeeze(sum(I_PA,3));%         %figure(1), clf, plot(1:si,squeeze(I_PA(:,4,:))), return%%        PeakPick        %figure(1), clf, plot(1:numel(PP.Apeak(:,1)),PP.Apeak(:,1)), return%%        gamma = 26.75e7;        if strcmp(NMRacqus.nuc1,'2H') == 1            gamma = 4.1065e7;        elseif strcmp(NMRacqus.nuc1,'23Na') == 1            gamma = 7.0761e7;        end        Gmax = 3;        if any(strcmp(NMRacqus.probhd,{'5 mm BBO BB-1H/D XYZ-GRD Z107255/0001',...                '5 mm TXI 1H/D-13C/15N XYZ-GRD Z8588/0006'})) == 1            Gmax = 0.5;        end                %gradient ramps in indirect dimension        fid = fopen([num2str(expno(nexp)) '/rax.txt']);        ramp.ax = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/ray.txt']);        ramp.ay = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/raz.txt']);        ramp.az = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rbx.txt']);        ramp.bx = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rby.txt']);        ramp.by = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rbz.txt']);        ramp.bz = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rcx.txt']);        ramp.cx = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rcy.txt']);        ramp.cy = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rcz.txt']);        ramp.cz = fscanf(fid,'%f');        fclose(fid);        G.ax = ramp.ax*Gmax.*NMRacqus.cnst1/100;        G.bx = ramp.bx*Gmax.*NMRacqus.cnst1/100;        G.cx = ramp.cx*Gmax.*NMRacqus.cnst1/100;        G.ay = ramp.ay*Gmax.*NMRacqus.cnst2/100;        G.by = ramp.by*Gmax.*NMRacqus.cnst2/100;        G.cy = ramp.cy*Gmax.*NMRacqus.cnst2/100;        G.az = ramp.az*Gmax.*NMRacqus.cnst3/100;        G.bz = ramp.bz*Gmax.*NMRacqus.cnst3/100;        G.cz = ramp.cz*Gmax.*NMRacqus.cnst3/100;                %symmetry vector of b-matrix        symv.x = G.ax + G.bx + G.cx;        symv.y = G.ay + G.by + G.cy;        symv.z = G.az + G.bz + G.cz;        symv.norm = sqrt(symv.x.^2 + symv.y.^2 + symv.z.^2);        symv.x = symv.x./symv.norm;        symv.y = symv.y./symv.norm;        symv.z = symv.z./symv.norm;        symv.theta = acos(symv.z);        symv.phi = atan2(symv.y,symv.x);%         figure(1), clf%         subplot(1,2,1)%         plot(symv.theta,symv.phi,'o')%         subplot(1,2,2)% %         for n = 1:Ndir%                 nG = (n-1)*td1/Ndir + td1/Ndir/NDeltab;%                 plot3([0; G.ax(nG)],[0; G.ay(nG)],[0; G.az(nG)],'-o',[0; G.bx(nG)],[0; G.by(nG)],[0; G.bz(nG)],'-s',[0; G.cx(nG)],[0; G.cy(nG)],[0; G.cz(nG)],'-x')%                 hold on%                 plot3([0; symv.x(nG)],[0; symv.y(nG)],[0; symv.z(nG)],'k-o')%                 axis(1*[-1 1 -1 1 -1 1])%                 axis square  %         end%         pause(.1)%%        epsilon = NMRacqus.d2;        delta = 2*(NMRacqus.d2 + NMRacqus.d3);        Delta = 2*NMRacqus.d2 + 2*NMRacqus.d3 + 2*NMRacqus.d4 + 2*NMRacqus.d6 + NMRacqus.d5;                    if NMRacqus.l11            Delta = Delta + NMRacqus.d46 + 2*NMRacqus.d42 + NMRacqus.d43;        end        tdiff = Delta - delta/3;        %diffusion weighting matrix b        bmat.xx = gamma^2*delta^2*tdiff*(G.ax.*G.ax + G.bx.*G.bx + G.cx.*G.cx);        bmat.xy = gamma^2*delta^2*tdiff*(G.ax.*G.ay + G.bx.*G.by + G.cx.*G.cy);        bmat.xz = gamma^2*delta^2*tdiff*(G.ax.*G.az + G.bx.*G.bz + G.cx.*G.cz);        bmat.yy = gamma^2*delta^2*tdiff*(G.ay.*G.ay + G.by.*G.by + G.cy.*G.cy);        bmat.yz = gamma^2*delta^2*tdiff*(G.ay.*G.az + G.by.*G.bz + G.cy.*G.cz);        bmat.zz = gamma^2*delta^2*tdiff*(G.az.*G.az + G.bz.*G.bz + G.cz.*G.cz);                b = (bmat.xx + bmat.yy + bmat.zz)';        %figure(1), clf, plot((1:td1)',[b],'o'), return                %consistency check of the b-matrix eigenvalues        G.a = sqrt(G.ax.^2 + G.ay.^2 + G.az.^2);        G.b = sqrt(G.bx.^2 + G.by.^2 + G.bz.^2);        G.c = sqrt(G.cx.^2 + G.cy.^2 + G.cz.^2);        %figure(1), clf, plot((1:td1)',[G.c G.b G.a],'o'), return        %figure(1), clf, plot((1:td1)',[G.angle],'-o'), return%         G.r = abs(G.c + 1i*G.b);%         G.angle = angle(G.c + 1i*G.b);%         figure(1), clf, plot((1:td1)',[G.r],'-o'), return% %         F.a = cumsum(Gmod.a*dt);%         F.b = cumsum(Gmod.b*dt);%         F.c = cumsum(Gmod.c*dt);% %         b_check = 2*sum(F.c.^2*dt)*gamma^2*G.r.^2;        %figure(1), clf, plot(b,b_check), return        lambda1 = zeros(1,length(G.a));        lambda2 = zeros(1,length(G.a));        lambda3 = zeros(1,length(G.a));        for ntd1 = 1:length(G.a)            lambdas = eig([bmat.xx(ntd1) bmat.xy(ntd1) bmat.xz(ntd1)            bmat.xy(ntd1) bmat.yy(ntd1) bmat.yz(ntd1)            bmat.xz(ntd1) bmat.yz(ntd1) bmat.zz(ntd1)]);            lambda.mean(1,ntd1) = mean(lambdas);            Dlambdas = abs(lambdas-lambda.mean(1,ntd1));            [dummy,indx] = sort(Dlambdas,'descend');            lambda.c(1,ntd1) = lambdas(indx(1));            lambda.b(1,ntd1) = lambdas(indx(2));            lambda.a(1,ntd1) = lambdas(indx(3));                                    lambda.iso(1,ntd1) = 3*min(lambdas);                                end                lambda.Delta = (lambda.c - (lambda.b+lambda.a)/2)/3./lambda.mean;        lambda.eta = (lambda.b - lambda.a)./lambda.mean;        lambda.eA = lambda.c - lambda.iso/3;        lambda.eR = 2*((lambda.a+lambda.b)/2 - lambda.iso/3);        lambda.aniso = lambda.eA - lambda.eR;                G.angle = acos(sqrt((2*lambda.Delta + 1)/3));                bmat.b = b';        bmat.zeta = G.angle';        bmat.Delta = lambda.Delta';        bmat.eta = lambda.eta';        bmat.iso = lambda.iso';        bmat.aniso = lambda.aniso';        bmat.dir.x = symv.x;        bmat.dir.y = symv.y;        bmat.dir.z = symv.z;        bmat.dir.theta = symv.theta;        bmat.dir.phi = symv.phi;        %figure(1), clf, plot((1:Ndt)',G.x,'-'), return        %figure(1), clf, plot((1:td1)',[lambda.c' lambda.b' lambda.a' b'],'-'), return        %figure(1), clf, plot((1:td1)',[lambda.Delta' lambda.eta'],'-'), return%%        clear FitDat        FitDat.fitpoints = 1:Nb;        AcqDat = struct('bmat',bmat,'td1',td1,'Nb',Nb,'Ndir',Ndir,'Nvd',Nvd);               for npeak = 1:PP.Npeaks            Apeak_PA = reshape(PP.Apeak(:,npeak),[Nb Ndir Nvd]);            Apeak_PA = squeeze(sum(Apeak_PA,2));            Apeak_PA = reshape(Apeak_PA,Nb*Nvd,1);            PP.Apeak_PA(:,npeak) = Apeak_PA;            Ipeak_PA = reshape(PP.Ipeak(:,npeak),[Nb Ndir Nvd]);            Ipeak_PA = squeeze(sum(Ipeak_PA,2));            Ipeak_PA = reshape(Ipeak_PA,Nb*Nvd,1);            PP.Ipeak_PA(:,npeak) = Ipeak_PA;        end        %figure(1), clf, plot(FitDat.fitpoints, G.angle(FitDat.fitpoints),'-o'), return%         anglemin = min(G.angle);%         tol = 1e-3;%         indx_anglemin = find([G.angle(FitDat.fitpoints)>(anglemin-tol)...%             & G.angle(FitDat.fitpoints)<(anglemin+tol)]);%         FitDat.fitpoints = FitDat.fitpoints(indx_anglemin);        %PlotInterm = 1;        for npeak = 1:PP.Npeaks            %npeak = 3;            Xin = b(FitDat.fitpoints)'; Xin2 = G.angle(FitDat.fitpoints)';            Yin = PP.Ipeak_PA(FitDat.fitpoints,npeak);             if strcmp(signal,'area')==1                Yin = PP.Apeak_PA(FitDat.fitpoints,npeak);            end            %figure(1), clf, semilogy(Xin,Yin,'o'), return            Pin = [Yin(1)*1.05   5/mean(Xin(:,1))*[.01 1]]; Funam = 'fDiffVACSY';             Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = Pin;             Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,zeros(size(Pin)),[],[],Xin2,Pnorm,Xnorm,Ynorm);            Yout = feval(Funam,Pout,Xin,Xin2); error = Yin - Yout;            Pout_oblate = Pout;            Yout_oblate = Yout;            chisq_oblate = mean(error.^2);                        Dpar_oblate = Pout(2);            Dperp_oblate = Pout(3);            MD_oblate = (Dpar_oblate + 2*Dperp_oblate)/3;            DD_oblate = Dpar_oblate - Dperp_oblate;                        MD_prolate = MD_oblate;            DD_prolate = -DD_oblate;            Dpar_prolate = MD_oblate+2*DD_prolate/3;            Dperp_prolate = MD_oblate-DD_prolate/3;            Pin = [Pout(1) Dpar_prolate Dperp_prolate];            Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = Pin;             Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,zeros(size(Pin)),[],[],Xin2,Pnorm,Xnorm,Ynorm);            Yout = feval(Funam,Pout,Xin,Xin2); error = Yin - Yout;            chisq_prolate = mean(error.^2);                        if chisq_oblate < chisq_prolate                Pout = Pout_oblate;                Yout = Yout_oblate;            end                        if PlotInterm%%                figure(11), clf                axes('position',[.1 .27 .8 .65])                semilogy(Xin,Yin,'o')                hold on                semilogy(Xin,Yout,'-')                title(['expno=' num2str(expno(nexp)) '  peak=' num2str(npeak) '  '  Funam  '   Pout= ' num2str(Pout,3) ])                axis([min(min(Xin)) max(max(Xin)) max(max(Yout))*[1e-2 1.2] ])                xlabel('time'), ylabel('intensity')                axes('position',[.1 .05 .8 .1])                plot(Xin,error,'o'), grid                ylabel('residual')                pause(1)                %return            end            FitDat.Xin{npeak} = Xin;            FitDat.Xin2{npeak} = Xin2;            FitDat.Yin{npeak} = Yin;            FitDat.Yout{npeak} = Yout;            FitDat.Y0(:,npeak) = Pout(1);            FitDat.Dpar(:,npeak) = Pout(2);            FitDat.Dperp(:,npeak) = Pout(3);        end        %%        fs = 20;        lw = 1;        figure(1), clf        axes('position',[-.02 .18 .8 .3],'FontSize',fs*.8)        plot(PP.ppm,PP.Ispec,'k','LineWidth',lw)        set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,'Box','off','TickDir','out','Ycolor',[1 1 1])        axis('tight')        ylim = get(gca,'YLim');        ylim = .05*diff(ylim)*[-1 1] + ylim;        set(gca,'YLim',ylim)        xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'],'FontSize',fs)        hold on        plot(PP.peakppm(PP.td1start,:),PP.Ipeak(PP.td1start,:),'bo','LineWidth',lw)        [X,Y] = fSchmidt(symv.x,symv.y,symv.z);        latitude.theta = pi/180*[30:30:150 179];        latitude.phi = linspace(0,2*pi,100);        [latitude.phi,latitude.theta] = ndgrid(latitude.phi,latitude.theta);        latitude.z = cos(latitude.theta);        latitude.x = sin(latitude.theta).*cos(latitude.phi);        latitude.y = sin(latitude.theta).*sin(latitude.phi);        [latitude.X,latitude.Y] = fSchmidt(latitude.x,latitude.y,latitude.z);        longitude.theta = pi/180*linspace(30,180,100);        longitude.phi = pi/180*[30:30:360];        [longitude.theta,longitude.phi] = ndgrid(longitude.theta,longitude.phi);        longitude.z = cos(longitude.theta);        longitude.x = sin(longitude.theta).*cos(longitude.phi);        longitude.y = sin(longitude.theta).*sin(longitude.phi);        [longitude.X,longitude.Y] = fSchmidt(longitude.x,longitude.y,longitude.z);        axes('position',[.79 .18 .2 .3],'FontSize',fs*.8)        plot(X,Y,'ko')        hold on        plot(latitude.X,latitude.Y,'b-')        plot(longitude.X,longitude.Y,'b-')        axis tight equal off        dleft = .93/Npeaks;        width = .9*dleft;        height = .3;        left = 1-width;        bottom = .55;        if Npeaks == 1, width = .8*1/2; end        for npeak = 1:Npeaks            axes('position',[left bottom width/2 height],'FontSize',fs*.5)            Xval = FitDat.Xin{npeak};            X2val = (3*cos(FitDat.Xin2{npeak}).^2 - 1)/2;            Yinval = FitDat.Yin{npeak}/FitDat.Y0(npeak);            Youtval = FitDat.Yout{npeak}/FitDat.Y0(npeak);             Xval = reshape(Xval,sqrt(Nb),sqrt(Nb));            X2val = reshape(X2val,sqrt(Nb),sqrt(Nb));            Yinval = reshape(Yinval,sqrt(Nb),sqrt(Nb));            Youtval = reshape(Youtval,sqrt(Nb),sqrt(Nb));            bzz = (2*X2val+1)/3.*Xval;            bxx = (Xval - bzz)/2;            bmin_array = zeros(sqrt(Nb),sqrt(Nb),2);            bmin_array(:,:,1) = bzz;            bmin_array(:,:,2) = bxx;            biso = 3*min(bmin_array,[],3);            beA = bzz - biso/3;            beR = 2*(bxx - biso/3);            baniso = beA - beR;                        Xval = biso;            X2val = baniso;                        semilogy(Xval,Yinval,'bo','LineWidth',lw)            hold on            semilogy(Xval,Youtval,'b-','LineWidth',lw)            axis('tight')            ylim = get(gca,'YLim');            %ylim = .05*abs(diff(ylim))*[-1 1] + abs(ylim); ylim(1) = 0;            %ylim = ylim(2)*[1e-2 1.2];            ylim = 1*[Imin 1.2];            %xlim = get(gca,'XLim');            %xlim = .05*diff(xlim)*[-1 1] + xlim;            xlim = max(reshape(Xval,numel(Xval),1))*[-.1 1.1];            set(gca,'XLim',xlim,'YLim',ylim,'Box','off','YTick',[1e-3 1e-2 1e-1 1],...                'TickDir','out','TickLength',.03*[1 1],'Box','off','LineWidth',1.5*lw)            if npeak == Npeaks                xlabel('b_{iso} / m^-^2s','FontSize',fs*.6)                ylabel('I / I_0','FontSize',fs*.6)            end            if npeak < Npeaks                set(gca,'XTickLabel',[],'YTickLabel',[])            end            Dpar = FitDat.Dpar(npeak);            Dperp = FitDat.Dperp(npeak);            MD = (Dpar + 2*Dperp)/3;            title({[num2str(PP.peakppm(PP.td1start,npeak),3) ' ppm'];...                ['MD=' num2str(MD,2) ' m^2s^-^1']},'FontSize',fs*.5)            axes('position',[left+width/2 bottom width/2 height],'FontSize',fs*.5)            semilogy(X2val',Yinval','bo','LineWidth',lw)            hold on            semilogy(X2val',Youtval','b-','LineWidth',lw)            set(gca,'XLim',xlim,'YLim',ylim,'Box','off','YTick',[1e-3 1e-2 1e-1 1],...                'TickDir','out','TickLength',.03*[1 1],'Box','off','LineWidth',1.5*lw)            set(gca,'YTickLabel',[])            title({['AD=' num2str(Dpar,2) ' m^2s^-^1'];...                [' RD=' num2str(Dperp,2) ' m^2s^-^1']},'FontSize',fs*.5)            if npeak == Npeaks                xlabel('b_{aniso} / m^-^2s','FontSize',fs*.6)            end            if npeak < Npeaks                set(gca,'XTickLabel',[])            end            left = left-dleft;        end%%                                    eval(['save ' num2str(expno(nexp)) '/FitDat PP FitDat AcqDat'])        delete([num2str(expno(nexp)) '/ReportFig.*'])        eval(['print -depsc -loose ' num2str(expno(nexp)) '/ReportFig'])        %%    end    end