clear allwd = cd;%DataDir = '/Users/daniel/Dropbox';DataDir = '/Users/daniel/NMRdata/AVII500/DT/';%DataDir = '/opt/topspin2/data/DT/nmr';cd(DataDir)%ExpNam = {'DiffVACSY'}; expno = 124;%basl = [.02 .1 .65 .98];%lb = 20;%si = 1*1024;%ll = 1959; ul = 2307;% ExpNam = {'AOT_1H2H'};%expno = 110; GridSamp = 1; Ndir = 7; NDeltab = 4;%expno = 111; GridSamp = 1; Ndir = 30; NDeltab = 4;%expno = 112; GridSamp = 1; Ndir = 30; NDeltab = 7;%expno = 113; GridSamp = 1; Ndir = 30; NDeltab = 10;%expno = 114; GridSamp = 1; Ndir = 30; NDeltab = 13;%expno = 115; GridSamp = 1; Ndir = 30; NDeltab = 4;%expno = 116; GridSamp = 1; Ndir = 30; NDeltab = 4;%expno = 117; GridSamp = 1; Ndir = 30; NDeltab = 4;%expno = 118; GridSamp = 1; Ndir = 30; NDeltab = 4;%expno = 119; GridSamp = 1; Ndir = 39; NDeltab = 13;%expno = 120; GridSamp = 1; Ndir = 39; NDeltab = 13;%basl = [.02 .2 .4 .6 .8 .98];%basl = [.5 .63 .77 .95];%lb = 10; si = 2*1024;%phc0 = 200; phc1 = -70;%phc0 = 160; phc1 = 0;%expno = 123; GridSamp = 1; Ndir = 15; NDeltab = 4;%expno = 124; GridSamp = 1; Ndir = 15; NDeltab = 4;%expno = 125; GridSamp = 1; Ndir = 7; NDeltab = 4;%phc0 = 150; phc1 = -70;%expno = 126; GridSamp = 1; Ndir = 7; NDeltab = 4;%expno = 128; GridSamp = 1; Ndir = 39; NDeltab = 4;%basl = [.5 .63 .77 .95];%lb = 10; si = 2*1024;%phc0 = 130; phc1 = -70;%expno = 130; GridSamp = 1; Ndir = 5; NDeltab = 4;%expno = 131; GridSamp = 1; Ndir = 5; NDeltab = 4;%expno = 132; GridSamp = 1; Ndir = 5; NDeltab = 4;%expno = 133; GridSamp = 1; Ndir = 5; NDeltab = 4;%expno = 134; GridSamp = 1; Ndir = 39; NDeltab = 4;%expno = 135; GridSamp = 1; Ndir = 39; NDeltab = 13;%expno = 204; GridSamp = 1; Ndir = 5; NDeltab = 1;%expno = 205; GridSamp = 1; Ndir = 5; NDeltab = 1;%expno = 206; GridSamp = 1; Ndir = 5; NDeltab = 1;%phc0 = 180; phc1 = -70;%expno = 207; GridSamp = 1; Ndir = 39; NDeltab = 13;%expno = 208; GridSamp = 1; Ndir = 39; NDeltab = 13;% expno = 209; GridSamp = 1; Ndir = 39; NDeltab = 13;%expno = 210; GridSamp = 1; Ndir = 39; NDeltab = 13;% expno = 211; GridSamp = 1; Ndir = 39; NDeltab = 13;% basl = [.5 .61 .77 .95];% lb = 10; si = 2*1024;% phc0 = 165; phc1 = -70;%expno = 160; GridSamp = 1; Ndir = 5; NDeltab = 4;%expno = 161; GridSamp = 1; Ndir = 5; NDeltab = 4; FIDeff = .7; %expno = 162; GridSamp = 1; Ndir = 39; NDeltab = 4; FIDeff = .5; %expno = 163; GridSamp = 1; Ndir = 39; NDeltab = 4;% expno = 165; GridSamp = 1; Ndir = 39; NDeltab = 13; % expno = 171; GridSamp = 1; Ndir = 39; NDeltab = 13; % basl = [.4 .55 .8 .95];% lb = 50; si = 1*1024;% phc0 = 180; phc1 = 0;%expno = 185; GridSamp = 1; Ndir = 5; NDeltab = 1; %expno = 186; GridSamp = 1; Ndir = 5; NDeltab = 1; %expno = 187; GridSamp = 1; Ndir = 5; NDeltab = 1; %expno = 188; GridSamp = 1; Ndir = 5; NDeltab = 1; %expno = 189; GridSamp = 1; Ndir = 5; NDeltab = 1; %expno = 190; GridSamp = 1; Ndir = 5; NDeltab = 1; %expno = 191; GridSamp = 1; Ndir = 5; NDeltab = 4; %expno = 192; GridSamp = 1; Ndir = 5; NDeltab = 4; %expno = 193; GridSamp = 1; Ndir = 5; NDeltab = 4; %expno = 194; GridSamp = 1; Ndir = 5; NDeltab = 4; %expno = 195; GridSamp = 1; Ndir = 39; NDeltab = 4; % expno = 196; GridSamp = 1; Ndir = 39; NDeltab = 13; % basl = [.61 .64 .72 .75]; skip = [0 .61 .67 1];% lb = 1; si = 8*1024;% phc0 = 250; phc1 = -90; pivot = .66;%expno = 224; GridSamp = 1; Ndir = 39; NDeltab = 13; % expno = 225; GridSamp = 1; Ndir = 39; NDeltab = 13; % basl = [.64 .65 .73 .75]; %skip = [0 .61 .67 1];% lb = 1; si = 16*1024;% phc0 = 277; phc1 = -90; pivot = .6752;% AutoPhase = 0; %AutoPhaseAll = 0;% ExpNam = {'qVAS_Starch'};% expno = 15; GridSamp = 1; Ndir = 1; NDeltab = 4; % expno = 16; GridSamp = 1; Ndir = 1; NDeltab = 4; % expno = 17; GridSamp = 1; Ndir = 1; NDeltab = 4; % expno = 18; GridSamp = 1; Ndir = 1; NDeltab = 4; % expno = 19; GridSamp = 1; Ndir = 1; NDeltab = 4; % expno = 20; GridSamp = 1; Ndir = 5; NDeltab = 7; % expno = 26; GridSamp = 1; Ndir = 1; NDeltab = 4; % expno = 27; GridSamp = 1; Ndir = 5; NDeltab = 4; % expno = 28; GridSamp = 1; Ndir = 5; NDeltab = 4; % expno = 29; GridSamp = 1; Ndir = 15; NDeltab = 7; % expno = 30; GridSamp = 1; Ndir = 15; NDeltab = 7; % expno = 31; GridSamp = 1; Ndir = 15; NDeltab = 7; % expno = 32; GridSamp = 1; Ndir = 15; NDeltab = 7; % expno = 34; GridSamp = 1; Ndir = 15; NDeltab = 7; % expno = 35; GridSamp = 1; Ndir = 7; NDeltab = 4; % expno = 36; GridSamp = 1; Ndir = 7; NDeltab = 4; % expno = 37; GridSamp = 1; Ndir = 7; NDeltab = 4; % expno = 59; GridSamp = 1; Ndir = 7; NDeltab = 4; % expno = 60; GridSamp = 1; Ndir = 7; NDeltab = 4; % expno = 61; GridSamp = 1; Ndir = 7; NDeltab = 4; % expno = 62; GridSamp = 1; Ndir = 7; NDeltab = 4; % expno = 63; GridSamp = 1; Ndir = 15; NDeltab = 1; % expno = 64; GridSamp = 1; Ndir = 15; NDeltab = 4; % expno = 65; GridSamp = 1; Ndir = 55; NDeltab = 13; % basl = [.02 .4 .6 .98]; %skip = [0 .61 .67 1];% phc0 = 30; phc1 = 0; pivot = .5;% lb = 20; si = 512;%AutoPhase = 1; %AutoPhaseAll = 0;% ExpNam = {'AOTisooctane'};% expno = 31; GridSamp = 1; Ndir = 7; NDeltab = 4; % expno = 32; GridSamp = 1; Ndir = 17; NDeltab = 7; % basl = [.05 .15 .4 .6 .85 .95]; %skip = [0 .61 .67 1];% AutoPhase = 1; %phc0 = 30; phc1 = 0; pivot = .5;% lb = 5; si = 1024;% expno = 205; GridSamp = 1; Ndir = 33; NDeltab = 7; % lb = 5; si = 4*1024;% ExpNam = {'AOToct_temp'};% expno = 13:10:103;% expno= 13;% GridSamp = 1; Ndir = 15; NDeltab = 4; % basl = [.02 .2 .4 .6 .8 .98]; %skip = [0 .61 .67 1];% AutoPhase = 1; %phc0 = 30; phc1 = 0; pivot = .5;% lb = 5; si = 2*1024;% % ExpNam = {'AOToct_temp2'};% expno = 13:10:213;% %expno = 13;% GridSamp = 1; Ndir = 15; NDeltab = 4; % basl = [.02 .2 .4 .6 .8 .98]; %skip = [0 .61 .67 1];% AutoPhase = 1; %phc0 = 30; phc1 = 0; pivot = .5;% lb = 5; si = 2*1024;% ExpNam = {'DatDMPC'};% expno = 3;% GridSamp = 1; Ndir = 7; NDeltab = 4;% lb = 100; basl = [.02 .3 .7 .98];% AutoPhase = 1;% % ExpNam = {'AOToct_temp3'};% expno = 13:10:83;% expno = 93:10:113;% expno = 143:10:183;% GridSamp = 1; Ndir = 15; NDeltab = 4; % basl = [.02 .2 .4 .6 .8 .98]; %skip = [0 .61 .67 1];% AutoPhase = 1; %phc0 = 30; phc1 = 0; pivot = .5;% lb = 5; si = 2*1024;% ExpNam = {'AOToct'};% expno = 103; expno = 113;% ExpNam = {'AOToct_temp4'}; expno = [6 16 26 36];% GridSamp = 1; Ndir = 15; NDeltab = 4; % basl = [.02 .15 .4 .6 .85 .98]; %skip = [0 .61 .67 1];% AutoPhase = 1; phc0 = 75; phc1 = 0; pivot = .5;% lb = 30; si = 2*1024;% ExpNam = {'AOToct2'}; expno = 26; expno = 27;% GridSamp = 1; Ndir = 15; NDeltab = 4; % basl = [.03 .2 .8 .97]; lb = 5; si = 4*1024;% ll = [2698 2783 2825]; ul = [2708 2793 2835];% AutoPhase = 1; AutoPhaseAll = 0;% ExpNam = {'AOT_1H2H_comp'}; expno = 6; expno = 7;% GridSamp = 1; Ndir = 5; NDeltab = 3; % basl = [.03 .2 .8 .97]; lb = 20; si = 4*1024;% %ll = [2698 2783 2825]; ul = [2708 2793 2835];% phc0 = 270; phc1 = -80; AutoPhase = 0; %AutoPhaseAll = 0;% % ExpNam = {'AOT_1H2H_comp'}; % expno = 10; GridSamp = 1; Ndir = 17; NDeltab = 3; % expno = 11; GridSamp = 1; Ndir = 17; NDeltab = 3; % expno = 12; GridSamp = 1; Ndir = 17; NDeltab = 3; % expno = 13; GridSamp = 1; Ndir = 39; NDeltab = 3; % expno = 15; GridSamp = 1; Ndir = 17; NDeltab = 3; % expno = 17; GridSamp = 1; Ndir = 17; NDeltab = 4; % expno = 18; GridSamp = 1; Ndir = 39; NDeltab = 4; % basl = [.03 .2 .8 .97]; lb = 20; si = 4*1024;% %ll = [2698 2783 2825]; ul = [2708 2793 2835];% phc0 = 230; phc1 = -80; AutoPhase = 0; %AutoPhaseAll = 0;% % ExpNam = {'AOToct7'}; expno = 8;% GridSamp = 1; Ndir = 17; NDeltab = 4; % basl = [.03 .2 .8 .97]; lb = 20; si = 4*1024;% %ll = [2698 2783 2825]; ul = [2708 2793 2835];% phc0 = 230; phc1 = -80; AutoPhase = 1; %AutoPhaseAll = 0;% ExpNam = {'AOToct_temp9'}; expno = 16:10:136; expno = [146];% GridSamp = 1; Ndir = 17; NDeltab = 4; % basl = [.03 .2 .8 .97]; lb = 20; si = 4*1024;% %ll = [2698 2783 2825]; ul = [2708 2793 2835];% phc0 = 230; phc1 = -80; AutoPhase = 1; %AutoPhaseAll = 0;% ExpNam = {'AOTdecanol'};% expno = 5; GridSamp = 1; Ndir = 17; NDeltab = 4; % expno = 6; GridSamp = 1; Ndir = 7; NDeltab = 4; % expno = 7; GridSamp = 1; Ndir = 17; NDeltab = 4; % expno = 8; GridSamp = 1; Ndir = 39; NDeltab = 4; % expno = 10:12; GridSamp = 1; Ndir = 17; NDeltab = 4; % basl = [.03 .2 .8 .97]; lb = 20; si = 4*1024;% %ll = [2698 2783 2825]; ul = [2708 2793 2835];% phc0 = 210; phc1 = -160; AutoPhase = 1; %AutoPhaseAll = 0;ExpNam = {'Yeast_tripgse'};%expno = 11; GridSamp = 1; Ndir = 7; NDeltab = 4; %expno = 12; GridSamp = 1; Ndir = 7; NDeltab = 4; %expno = 13:17; GridSamp = 1; Ndir = 17; NDeltab = 4; expno = 18; GridSamp = 1; Ndir = 17; NDeltab = 4; %expno = 19; GridSamp = 1; Ndir = 17; NDeltab = 4; phc0 = 195; phc1 = 10; %expno = 20; GridSamp = 1; Ndir = 17; NDeltab = 4; phc0 = 65; phc1 = 10; %expno = 22; GridSamp = 1; Ndir = 17; NDeltab = 4; phc0 = 70; phc1 = -150; %expno = 24:27; GridSamp = 1; Ndir = 7; NDeltab = 4; %expno = 28; GridSamp = 1; Ndir = 17; NDeltab = 4; %expno = 35; GridSamp = 1; Ndir = 17; NDeltab = 4; %expno = 28:40; GridSamp = 1; Ndir = 17; NDeltab = 4; %expno = 43; GridSamp = 1; Ndir = 7; NDeltab = 4; %expno = 46:51; GridSamp = 1; Ndir = 17; NDeltab = 4; %expno = 61:64; GridSamp = 1; Ndir = 7; NDeltab = 2; basl = [.03 .15 .8 .97]; lb = 20; si = 2*1024;%ll = [2698 2783 2825]; ul = [2708 2793 2835];AutoPhase = 1; %AutoPhaseAll = 0;CheckBasline = 0;FindPeaks = 1;CheckPeaks = 0;PlotInterm = 0;thresh = .1;td1start = 2;Imin = .8e-2;signal = 'area'; %signal = 'intensity';cd(ExpNam{1})if exist('expno') == 0    GetExpnosendfor nexp = 1:length(expno)    ConvertAcqus = 'N';    ConvertProcs = 'N';    MakeTextfile = 'N';    if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0        ConvertAcqus = 'Y';    end    res = fExpnoInfo2(ConvertAcqus,MakeTextfile,ConvertProcs,expno(nexp));    eval(['load ' num2str(expno(nexp)) '/NMRacqus'])    if any(strcmp(NMRacqus.pulprog,{'DT_qVAS','DT_tbppgste'})) == 1        Spec2Dtd1        PeakPick%%        gamma = 26.75e7;        if strcmp(NMRacqus.nuc1,'2H') == 1            gamma = 4.1065e7;        elseif strcmp(NMRacqus.nuc1,'23Na') == 1            gamma = 7.0761e7;        end        Gmax = 3;        if any(strcmp(NMRacqus.probhd,{'5 mm BBO BB-1H/D XYZ-GRD Z107255/0001',...                '5 mm TXI 1H/D-13C/15N XYZ-GRD Z8588/0006'})) == 1            Gmax = 0.5;        end                %gradient ramps in indirect dimension        fid = fopen([num2str(expno(nexp)) '/rax.txt']);        ramp.ax = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/ray.txt']);        ramp.ay = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/raz.txt']);        ramp.az = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rbx.txt']);        ramp.bx = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rby.txt']);        ramp.by = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rbz.txt']);        ramp.bz = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rcx.txt']);        ramp.cx = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rcy.txt']);        ramp.cy = fscanf(fid,'%f');        fclose(fid);        fid = fopen([num2str(expno(nexp)) '/rcz.txt']);        ramp.cz = fscanf(fid,'%f');        fclose(fid);        G.ax = ramp.ax*Gmax.*NMRacqus.cnst1/100;        G.bx = ramp.bx*Gmax.*NMRacqus.cnst1/100;        G.cx = ramp.cx*Gmax.*NMRacqus.cnst1/100;        G.ay = ramp.ay*Gmax.*NMRacqus.cnst2/100;        G.by = ramp.by*Gmax.*NMRacqus.cnst2/100;        G.cy = ramp.cy*Gmax.*NMRacqus.cnst2/100;        G.az = ramp.az*Gmax.*NMRacqus.cnst3/100;        G.bz = ramp.bz*Gmax.*NMRacqus.cnst3/100;        G.cz = ramp.cz*Gmax.*NMRacqus.cnst3/100;                %symmetry vector of b-matrix        symv.x = G.ax + G.bx + G.cx;        symv.y = G.ay + G.by + G.cy;        symv.z = G.az + G.bz + G.cz;        symv.norm = sqrt(symv.x.^2 + symv.y.^2 + symv.z.^2);        symv.x = symv.x./symv.norm;        symv.y = symv.y./symv.norm;        symv.z = symv.z./symv.norm;        symv.theta = acos(symv.z);        symv.phi = atan2(symv.y,symv.x);%%                        if any(strcmp(NMRacqus.pulprog,{'DT_qVAS'})) == 1            %qVAS gradient time-modulation            fid = fopen([num2str(expno(nexp)) '/qVASa.txt']);            Gmod.a = fscanf(fid,'%f');            fclose(fid);            fid = fopen([num2str(expno(nexp)) '/qVASb.txt']);            Gmod.b = fscanf(fid,'%f');            fclose(fid);            fid = fopen([num2str(expno(nexp)) '/qVASc.txt']);            Gmod.c = fscanf(fid,'%f');            fclose(fid);            Ndt = length(Gmod.c);            tau = NMRacqus.d3;            dt = tau/Ndt;                    %figure(1), clf, plot((1:Ndt)',[Gmod.c Gmod.b Gmod.a],'-'), return            G.x = repmat(Gmod.a,[1, td1]).*repmat(G.ax',[Ndt, 1]) + ...                repmat(Gmod.b,[1, td1]).*repmat(G.bx',[Ndt, 1]) + ...                repmat(Gmod.c,[1, td1]).*repmat(G.cx',[Ndt, 1]);            G.y = repmat(Gmod.a,[1, td1]).*repmat(G.ay',[Ndt, 1]) + ...                repmat(Gmod.b,[1, td1]).*repmat(G.by',[Ndt, 1]) + ...                repmat(Gmod.c,[1, td1]).*repmat(G.cy',[Ndt, 1]);            G.z = repmat(Gmod.a,[1, td1]).*repmat(G.az',[Ndt, 1]) + ...                repmat(Gmod.b,[1, td1]).*repmat(G.bz',[Ndt, 1]) + ...                repmat(Gmod.c,[1, td1]).*repmat(G.cz',[Ndt, 1]);            %dephasing vector F            F.x = cumsum(G.x*dt,1);            F.y = cumsum(G.y*dt,1);            F.z = cumsum(G.z*dt,1);            F.r = sqrt(F.x.^2 + F.y.^2 + F.z.^2);            %diffusion weighting matrix b            %factor 2 from the double DW blocks            bmat.xx = 2*gamma^2*sum(F.x.*F.x*dt,1);            bmat.xy = 2*gamma^2*sum(F.x.*F.y*dt,1);            bmat.xz = 2*gamma^2*sum(F.x.*F.z*dt,1);            bmat.yy = 2*gamma^2*sum(F.y.*F.y*dt,1);            bmat.yz = 2*gamma^2*sum(F.y.*F.z*dt,1);            bmat.zz = 2*gamma^2*sum(F.z.*F.z*dt,1);        else            epsilon = NMRacqus.d2;            tau = 2*NMRacqus.d4 + 1e-6*NMRacqus.p2;            delta = 2*(NMRacqus.d2 + NMRacqus.d3);            Delta = 2*NMRacqus.d2 + 2*NMRacqus.d3 + 2*NMRacqus.d4 + 2*NMRacqus.d6 + NMRacqus.d5;                        if NMRacqus.l11                Delta = Delta + NMRacqus.d46 + 2*NMRacqus.d42 + NMRacqus.d43;            end            tdiff = Delta - delta/3 - tau/2 - epsilon/2 - epsilon^2/6/delta + epsilon^3/15/delta^2;                        %diffusion weighting matrix b            bmat.xx = gamma^2*delta^2*tdiff*(G.ax.*G.ax + G.bx.*G.bx + G.cx.*G.cx);            bmat.xy = gamma^2*delta^2*tdiff*(G.ax.*G.ay + G.bx.*G.by + G.cx.*G.cy);            bmat.xz = gamma^2*delta^2*tdiff*(G.ax.*G.az + G.bx.*G.bz + G.cx.*G.cz);            bmat.yy = gamma^2*delta^2*tdiff*(G.ay.*G.ay + G.by.*G.by + G.cy.*G.cy);            bmat.yz = gamma^2*delta^2*tdiff*(G.ay.*G.az + G.by.*G.bz + G.cy.*G.cz);            bmat.zz = gamma^2*delta^2*tdiff*(G.az.*G.az + G.bz.*G.bz + G.cz.*G.cz);        end                b = (bmat.xx + bmat.yy + bmat.zz)';        %figure(1), clf, plot((1:td1)',[b],'o'), return                %consistency check of the b-matrix eigenvalues        G.a = sqrt(G.ax.^2 + G.ay.^2 + G.az.^2);        G.b = sqrt(G.bx.^2 + G.by.^2 + G.bz.^2);        G.c = sqrt(G.cx.^2 + G.cy.^2 + G.cz.^2);        %figure(1), clf, plot((1:td1)',[G.c G.b G.a],'o'), return        lambda1 = zeros(1,td1);        lambda2 = zeros(1,td1);        lambda3 = zeros(1,td1);        for ntd1 = 1:td1            lambdas = eig([bmat.xx(ntd1) bmat.xy(ntd1) bmat.xz(ntd1)            bmat.xy(ntd1) bmat.yy(ntd1) bmat.yz(ntd1)            bmat.xz(ntd1) bmat.yz(ntd1) bmat.zz(ntd1)]);            lambda.mean(1,ntd1) = mean(lambdas);            Dlambdas = abs(lambdas-lambda.mean(1,ntd1));            [dummy,indx] = sort(Dlambdas,'descend');            lambda.c(1,ntd1) = lambdas(indx(1));            lambda.b(1,ntd1) = lambdas(indx(2));            lambda.a(1,ntd1) = lambdas(indx(3));                                    lambda.iso(1,ntd1) = 3*min(lambdas);                                end                lambda.Delta = (lambda.c - (lambda.b+lambda.a)/2)/3./lambda.mean;        lambda.eta = (lambda.b - lambda.a)./lambda.mean;        lambda.eA = lambda.c - lambda.iso/3;        lambda.eR = 2*((lambda.a+lambda.b)/2 - lambda.iso/3);        lambda.aniso = lambda.eA - lambda.eR;                G.angle = acos(sqrt((2*lambda.Delta + 1)/3));                bmat.b = b';        bmat.zeta = G.angle';        bmat.Delta = lambda.Delta';        bmat.eta = lambda.eta';        bmat.iso = lambda.iso';        bmat.aniso = lambda.aniso';        bmat.dir.x = symv.x;        bmat.dir.y = symv.y;        bmat.dir.z = symv.z;        bmat.dir.theta = symv.theta;        bmat.dir.phi = symv.phi;        %figure(1), clf, plot(bmat.b,3*lambda.mean), axis square, return        %figure(1), clf, plot((1:Ndt)',G.x,'-'), return        %figure(1), clf, plot((1:td1)',[lambda.c' lambda.b' lambda.a' b'],'-'), return        %figure(1), clf, plot((1:td1)',[lambda.Delta' lambda.eta'],'-'), return%%        NG = td1/Ndir;                            clear FitDat        if GridSamp            NG = td1/Ndir;            Nb = td1/Ndir/NDeltab;            FitDat.fitpoints = 1:NG;        else            bmin = min(b);            tol = 5e-1;            indx_bmin = find([b>bmin*(1-tol) & b<bmin*(1+tol)]);            Ndir = length(indx_bmin);            NG = td1/Ndir;            FitDat.fitpoints = td1start:NG;        end                AcqDat = struct('bmat',bmat,'td1',td1,'Nb',Nb,'Ndir',Ndir,'NDeltab',NDeltab);               for npeak = 1:PP.Npeaks            PP.Apeak_PA(:,npeak) = sum(reshape(PP.Apeak(:,npeak),NG,Ndir),2);            PP.Ipeak_PA(:,npeak) = sum(reshape(PP.Ipeak(:,npeak),NG,Ndir),2);        end        %PlotInterm = 1;        for npeak = 1:PP.Npeaks            %npeak = 2;            Xin = bmat.b(FitDat.fitpoints);            Xin2 = bmat.Delta(FitDat.fitpoints);            Yin = PP.Ipeak_PA(FitDat.fitpoints,npeak);             if strcmp(signal,'area')==1                Yin = PP.Apeak_PA(FitDat.fitpoints,npeak);            end            %figure(1), clf, semilogy(Xin,Yin,'o'), return            options = optimset('Display','off');            Pin = [Yin(1)*1.05   5/mean(Xin(:,1))*[.001 1]]; Funam = 'fDiffVACSY';             Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = abs(Pin);             Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,zeros(size(Pin)),[],options,Xin2,Pnorm,Xnorm,Ynorm);            Yout = feval(Funam,Pout,Xin,Xin2); error = Yin - Yout;            Pout_oblate = Pout;            Yout_oblate = Yout;            chisq_oblate = mean(error.^2);            %figure(1), clf, semilogy(Xin,Yin,'o',Xin,Yout,'x'), set(gca,'YLim',Pout(1)*[Imin, 1.1]), return                        Dpar_oblate = Pout(2);            Dperp_oblate = Pout(3);            Diso_oblate = (Dpar_oblate + 2*Dperp_oblate)/3;            DeltaD_oblate = (Dpar_oblate - Dperp_oblate)/3/Diso_oblate;                        Diso_prolate = Diso_oblate;            DeltaD_prolate = -DeltaD_oblate;            Dpar_prolate = Diso_prolate*(1 + 2*DeltaD_prolate);            Dperp_prolate = Diso_prolate*(1 - DeltaD_prolate);            Pin = [Pout(1) Dpar_prolate Dperp_prolate];            Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = abs(Pin);             Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,zeros(size(Pin)),[],options,Xin2,Pnorm,Xnorm,Ynorm);            Yout = feval(Funam,Pout,Xin,Xin2); error = Yin - Yout;            %figure(1), clf, semilogy(Xin,Yin,'o',Xin,Yout,'-'), set(gca,'YLim',Pout(1)*[Imin, 1.1]), return            chisq_prolate = mean(error.^2);                        if chisq_oblate < chisq_prolate                Pout = Pout_oblate;                Yout = Yout_oblate;            end                        if PlotInterm                figure(11), clf                axes('position',[.1 .27 .8 .65])                semilogy(Xin,Yin,'o')                hold on                semilogy(Xin,Yout,'-')                title(['expno=' num2str(expno(nexp)) '  peak=' num2str(npeak) '  '  Funam  '   Pout= ' num2str(Pout,3) ])                axis([min(min(Xin)) max(max(Xin)) max(max(Yout))*[1e-2 1.2] ])                xlabel('time'), ylabel('intensity')                axes('position',[.1 .05 .8 .1])                plot(Xin,error,'o'), grid                ylabel('residual')                pause(1)                %return            end            FitDat.Xin{npeak} = Xin;            FitDat.Xin2{npeak} = Xin2;            FitDat.Yin{npeak} = Yin;            FitDat.Yout{npeak} = Yout;            FitDat.Y0(:,npeak) = Pout(1);            FitDat.Dpar(:,npeak) = Pout(2);            FitDat.Dperp(:,npeak) = Pout(3);        end        %%        fs = 20;        lw = 1;        figure(1), clf        axes('position',[-.02 .18 .8 .3],'FontSize',fs*.8)        plot(PP.ppm,PP.Ispec,'k','LineWidth',lw)        set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,'Box','off','TickDir','out','Ycolor',[1 1 1])        axis('tight')        ylim = get(gca,'YLim');        ylim = .05*diff(ylim)*[-1 1] + ylim;        set(gca,'YLim',ylim)        xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'],'FontSize',fs)        hold on        plot(PP.peakppm(PP.td1start,:),PP.Ipeak(PP.td1start,:),'bo','LineWidth',lw)        [X,Y] = fSchmidt(symv.x,symv.y,symv.z);        latitude.theta = pi/180*[30:30:150 179];        latitude.phi = linspace(0,2*pi,100);        [latitude.phi,latitude.theta] = ndgrid(latitude.phi,latitude.theta);        latitude.z = cos(latitude.theta);        latitude.x = sin(latitude.theta).*cos(latitude.phi);        latitude.y = sin(latitude.theta).*sin(latitude.phi);        [latitude.X,latitude.Y] = fSchmidt(latitude.x,latitude.y,latitude.z);        longitude.theta = pi/180*linspace(30,180,100);        longitude.phi = pi/180*[30:30:360];        [longitude.theta,longitude.phi] = ndgrid(longitude.theta,longitude.phi);        longitude.z = cos(longitude.theta);        longitude.x = sin(longitude.theta).*cos(longitude.phi);        longitude.y = sin(longitude.theta).*sin(longitude.phi);        [longitude.X,longitude.Y] = fSchmidt(longitude.x,longitude.y,longitude.z);        axes('position',[.79 .13 .2 .3],'FontSize',fs*.8)        plot(X,Y,'ko')        hold on        plot(latitude.X,latitude.Y,'b-')        plot(longitude.X,longitude.Y,'b-')        axis tight equal off        title(['Ndir=' num2str(Ndir)],'FontSize',fs*.5)        dleft = .93/Npeaks;        width = .9*dleft;        height = .3;        left = 1-width;        bottom = .58;        if Npeaks == 1, width = .8*1/2; end        for npeak = 1:Npeaks            axes('position',[left bottom width/2 height],'FontSize',fs*.5)            Xval = FitDat.Xin{npeak};            X2val = FitDat.Xin2{npeak};            Yinval = FitDat.Yin{npeak}/FitDat.Y0(npeak);            Youtval = FitDat.Yout{npeak}/FitDat.Y0(npeak);            if GridSamp                                Xval = reshape(Xval,Nb,NDeltab);                X2val = reshape(X2val,Nb,NDeltab);                Yinval = reshape(Yinval,Nb,NDeltab);                Youtval = reshape(Youtval,Nb,NDeltab);            end            semilogy(Xval,Yinval,'bo','LineWidth',lw)            hold on            semilogy(Xval,Youtval,'b-','LineWidth',lw)            axis('tight')            ylim = get(gca,'YLim');            ylim = 1*[Imin 1.2];            xlim = max(b)*[-.1 1.1];            set(gca,'XLim',xlim,'YLim',ylim,'Box','off','YTick',[1e-3 1e-2 1e-1 1],...                'TickDir','out','TickLength',.03*[1 1],'Box','off','LineWidth',1.5*lw)            if npeak == Npeaks                xlabel('b / m^-^2s','FontSize',fs*.6)                ylabel('I / I_0','FontSize',fs*.6)            end            if npeak < Npeaks                set(gca,'XTickLabel',[],'YTickLabel',[])            end            Dpar = FitDat.Dpar(npeak);            Dperp = FitDat.Dperp(npeak);            MD = (Dpar + 2*Dperp)/3;            title({[num2str(PP.peakppm(PP.td1start,npeak),3) ' ppm'];...                ['MD=' num2str(MD,2) ' m^2s^-^1']},'FontSize',fs*.5)            axes('position',[left+width/2 bottom width/2 height],'FontSize',fs*.5)            semilogy(X2val',Yinval','bo','LineWidth',lw)            hold on            semilogy(X2val',Youtval','b-','LineWidth',lw)            xlim = [-.7 1.2];            set(gca,'XLim',xlim,'YLim',ylim,'Box','off','YTick',[1e-3 1e-2 1e-1 1],...                'TickDir','out','TickLength',.03*[1 1],'Box','off','LineWidth',1.5*lw)            set(gca,'YTickLabel',[])            title({['AD=' num2str(Dpar,2) ' m^2s^-^1'];...                [' RD=' num2str(Dperp,2) ' m^2s^-^1']},'FontSize',fs*.5)            if npeak == Npeaks                xlabel('\Delta_b','FontSize',fs*.6)            end            if npeak < Npeaks                set(gca,'XTickLabel',[])            end            left = left-dleft;        end%%                                    eval(['save ' num2str(expno(nexp)) '/FitDat PP FitDat AcqDat'])        delete([num2str(expno(nexp)) '/ReportFig.*'])        eval(['print -depsc -loose ' num2str(expno(nexp)) '/ReportFig'])        %%    end    end