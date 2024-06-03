clear all

% niftis_path = '/Users/daniel/Dropbox/NMRdata/dtd_example_py/rwi/pmaps_py';
% nifti_fn = fullfile(niftis_path,'dtd_s0.nii.gz');
% [I,nii_h] = mdm_nii_read(nifti_fn);
% hdr_py = nii_h;
% 
% niftis_path = '/Users/daniel/Dropbox/NMRdata/dtd_example_py/rwi/pmaps';
% nifti_fn = fullfile(niftis_path,'dtd_s0.nii.gz');
% [I,nii_h] = mdm_nii_read(nifti_fn);
% hdr_mat = nii_h;

pmap_name = 'dtd_mdiso_bin1.nii.gz';
pmap_name = 'dtd_mdiso_bin3.nii.gz';
pmap_name = 'dtd_fractions.nii.gz';
% pmap_name = 'dtd_s0.nii.gz';

niftis_path = '/Users/daniel/Dropbox/NMRdata/dtd_example_py/rwi/pmaps_py/96';
nifti_fn = fullfile(niftis_path,pmap_name);
[I_py,hdr_py] = mdm_nii_read(nifti_fn);

niftis_path = '/Users/daniel/Dropbox/NMRdata/dtd_example_py/rwi/pmaps';
nifti_fn = fullfile(niftis_path,pmap_name);
[I_mat,hdr_mat] = mdm_nii_read(nifti_fn);


I_py = double(I_py);
I_mat = double(I_mat);

sz = size(I_py);

fov.x = hdr_py.dim(2)*hdr_py.pixdim(2);
fov.y = hdr_py.dim(3)*hdr_py.pixdim(3);
fov.z = hdr_py.dim(4)*hdr_py.pixdim(4);

sz = sz(2:4);

xslice = ceil(sz(1)/2);
yslice = ceil(sz(2)/2);
zslice = ceil(sz(3)/2);

position_pmap_x_0 = [fov.x/(fov.x+fov.y) fov.y/(fov.y+fov.z) fov.y/(fov.x+fov.y) fov.z/(fov.y+fov.z)];
position_pmap_y_0 = [0 fov.y/(fov.y+fov.z) fov.x/(fov.x+fov.y) fov.z/(fov.y+fov.z)];
position_pmap_z_0 = [0 0 fov.x/(fov.x+fov.y) fov.y/(fov.y+fov.z)];

Nmaps = 3;


clim_temp = [0 1];
fh = figure(2); clf
set(fh,'Color','k','InvertHardCopy','off','Visible', 'on');
papersize = 8.3*[1*Nmaps (fov.y+fov.z)/(fov.y+fov.x)];
set(fh, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 

I_dif = I_py-I_mat;

varnams = {'I_py'; 'I_mat'; 'I_dif'};

for npmap = 1:length(varnams)
    eval(['I = ' varnams{npmap} ';'])
    position_pmap_x   = position_pmap_x_0.*[1/Nmaps 1 1/Nmaps 1] + [(npmap-1)/Nmaps 0 0 0];
    position_pmap_y   = position_pmap_y_0.*[1/Nmaps 1 1/Nmaps 1] + [(npmap-1)/Nmaps 0 0 0];
    position_pmap_z   = position_pmap_z_0.*[1/Nmaps 1 1/Nmaps 1] + [(npmap-1)/Nmaps 0 0 0];

    I_r = reshape(I(1,:,:,:),sz)/255;
    I_g = reshape(I(2,:,:,:),sz)/255;
    I_b = reshape(I(3,:,:,:),sz)/255; 

    im2d_x = cat(3,squeeze(I_r(xslice,:,:))',squeeze(I_g(xslice,:,:))',squeeze(I_b(xslice,:,:))');
    im2d_y = cat(3,squeeze(I_r(:,yslice,:))',squeeze(I_g(:,yslice,:))',squeeze(I_b(:,yslice,:))');
    im2d_z = cat(3,squeeze(I_r(:,:,zslice))',squeeze(I_g(:,:,zslice))',squeeze(I_b(:,:,zslice))');

    axh_pmap_x = axes('position',position_pmap_x);
    axh_pmap_y = axes('position',position_pmap_y);
    axh_pmap_z = axes('position',position_pmap_z);

    imagesc(axh_pmap_x,.99*im2d_x/clim_temp(2)) % Factor .99 to avoid black pixels in pdf
    imagesc(axh_pmap_y,.99*im2d_y/clim_temp(2))
    imagesc(axh_pmap_z,.99*im2d_z/clim_temp(2))
    if npmap == 3
        scale = 10;
        imagesc(axh_pmap_x,.99*im2d_x/clim_temp(2)*scale) % Factor .99 to avoid black pixels in pdf
        imagesc(axh_pmap_y,.99*im2d_y/clim_temp(2)*scale)
        imagesc(axh_pmap_z,.99*im2d_z/clim_temp(2)*scale)
    end

    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'YDir','normal')
    axis([axh_pmap_x; axh_pmap_y; axh_pmap_z],'off')
end

