function [t,g] = mkh_read_rotacfn(rotac_fn)

fid = fopen(rotac_fn,'r');

format = '%f32';
numperrow = 2;

datin = textscan(fid,format,'CommentStyle',{'@'},'Headerlines',0);
datin = datin{1};
Nrows = length(datin)/numperrow;
dat = reshape(datin,numperrow,Nrows)';

t = dat(:,1);
g = dat(:,2);

t = 1e-12*double(t); % From MD to SI units
g = double(g);
