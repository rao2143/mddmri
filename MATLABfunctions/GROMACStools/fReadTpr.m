function tpr = fReadTpr(tpr_fn)

tprdump_fn = fullfile(fileparts(tpr_fn),'tprdump.txt');

source_cmd = fSourceCommand()
status = unix([source_cmd '; gmx dump -s ' tpr_fn ' > ' tprdump_fn]);

fid = fopen(tprdump_fn, 'r');

templine = fgetl(fid);
while ~strcmp('topology:',templine)
    templine = fgetl(fid);
    templine = deblank(strjust(templine,'left'));
end

templine = fgetl(fid);
templine = fgetl(fid);
indxeq = find(templine == '=');
tempval = templine(indxeq+1:length(templine));
tpr.Natoms = str2num(tempval);
templine = fgetl(fid);
indxeq = find(templine == '=');
tempval = templine(indxeq+1:length(templine));
tpr.Nmoltypes = str2num(tempval);
tpr.molnams = cell(1,tpr.Nmoltypes);
tpr.Nmolspermoltype = zeros(1,tpr.Nmoltypes);
tpr.Natomspermoltype = zeros(1,tpr.Nmoltypes);

nmoltype = 0;
while nmoltype < tpr.Nmoltypes
    templine = strjust(templine,'left');
    while ~strcmp(templine(1:7),'moltype')
        templine = fgetl(fid);
        templine = strjust(templine,'left');
    end
    nmoltype = nmoltype + 1;
    indxcit = find(templine == '"');
    tpr.molnams{1,nmoltype} = templine((indxcit(1)+1):(indxcit(2)-1));
    templine = fgetl(fid);
    indxeq = find(templine == '=');
    tempval = templine(indxeq+1:length(templine));
    tpr.Nmolspermoltype(1,nmoltype) = str2num(tempval);
    templine = fgetl(fid);
end
nmoltype = 0;
while nmoltype < tpr.Nmoltypes
    templine = strjust(templine,'left');
    while ~strcmp(templine,strcat('moltype (',num2str(nmoltype),'):'))
        templine = fgetl(fid);
        templine = deblank(strjust(templine,'left'));
    end
    nmoltype = nmoltype + 1;
    templine = fgetl(fid); templine = fgetl(fid); templine = fgetl(fid);
    indxstart = find(templine == '(');
    indxend = find(templine == ')');
    tpr.Natomspermoltype(1,nmoltype) = str2num(templine((indxstart+1):(indxend-1)));
    templine = fgetl(fid);
end

molnr_v = [];
molnrmax = 0;
for nmoltype = 1:tpr.Nmoltypes
    molnr_temp = repmat(molnrmax+(1:tpr.Nmolspermoltype(nmoltype)),[tpr.Natomspermoltype(nmoltype) 1]);
    molnr_v = cat(1,molnr_v,molnr_temp(:));
    molnrmax = max(molnr_v);
end

tpr.molnr_v = molnr_v;

fclose(fid);
% delete(tprdump_fn)

if tpr.Natoms ~= sum(tpr.Nmolspermoltype.*tpr.Natomspermoltype)
    warning('Inconsistent number of atoms in fReadTpr.m')
end