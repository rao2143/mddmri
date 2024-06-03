clear all


nii_paths = cell(0);
nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/NIA_NIH/2mm_SoftGrads/test1234-3_DTproc/mdd/nopreproc/pmaps';

nii_ext = '.nii.gz';
nii_name = 'dtor1r2d_mdiso';

for nnii = 1:numel(nii_paths)
    nii_path = nii_paths{nnii};
    [~,data_name,~] = fileparts(fileparts(fileparts(nii_path)));

    nii_fn = fullfile(nii_path,[nii_name nii_ext]);

    opt = mdm_opt();
    opt = mplot_opt(opt);

    [I,nii_h] = mdm_nii_read(nii_fn);

    clim.mdiso = [0 3]*1e-9;
    clim_temp = clim.mdiso;

    figure(2), clf
    
    sz = size(I);

    clim.s2000 = .8*[0 1];
    fov.x = nii_h.dim(2)*nii_h.pixdim(2);
    fov.y = nii_h.dim(3)*nii_h.pixdim(3);
    fov.z = nii_h.dim(4)*nii_h.pixdim(4);
    
    roi_maxx = 63;
    roi_maxy = 57;
    roi_maxz = 30;
    
    im2d_x = squeeze(I(roi_maxx,:,:));
    im2d_y = squeeze(I(:,roi_maxy,:));
    im2d_z = squeeze(I(:,:,roi_maxz));

    axh_pmap_z = axes('position',[0 fov.z/(fov.y+fov.z) fov.x/(fov.x+fov.y) fov.y/(fov.y+fov.z)]);
    imagesc(axh_pmap_z,im2d_z')
    set(axh_pmap_z,'YDir','normal')

    axh_pmap_x = axes('position',[fov.x/(fov.x+fov.y) 0 fov.y/(fov.x+fov.y) fov.z/(fov.y+fov.z)]);
    imagesc(axh_pmap_x,im2d_x')
    set(axh_pmap_x,'YDir','normal')

    axh_pmap_y = axes('position',[0 0 fov.x/(fov.x+fov.y) fov.z/(fov.y+fov.z)]);
    imagesc(axh_pmap_y,im2d_y')
    set(axh_pmap_y,'YDir','normal')

    hold(axh_pmap_z,'on')
    phxz = plot(axh_pmap_z,roi_maxx*[1 1],[-1 sz(1)+1],'-');
    phyz = plot(axh_pmap_z,[-1 sz(2)+1],roi_maxy*[1 1],'-');
    hold(axh_pmap_y,'on')
    phzy = plot(axh_pmap_y,[-1 sz(2)+1],roi_maxz*[1 1],'-');
    phxy = plot(axh_pmap_y,roi_maxx*[1 1],[-1 sz(3)+1],'-');
    hold(axh_pmap_x,'on')
    phzx = plot(axh_pmap_x,[-1 sz(1)+1],roi_maxz*[1 1],'-');
    phyx = plot(axh_pmap_x,roi_maxy*[1 1],[-1 sz(3)+1],'-');
    set([phzy; phzx; phxy; phxz; phyz; phyx],'Color',[1 .2 .2])

    th_z = text(axh_pmap_z,1,1,num2str(roi_maxz));
    th_y = text(axh_pmap_y,1,1,num2str(roi_maxy));
    th_x = text(axh_pmap_x,1,1,num2str(roi_maxx));
    set([th_z; th_y; th_x],'Color',[1 1 .999],'VerticalAlignment','bottom')

    colormap(axh_pmap_x,'gray')
    colormap(axh_pmap_y,'gray')
    colormap(axh_pmap_z,'gray')
    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',clim_temp)
    axis([axh_pmap_x; axh_pmap_y; axh_pmap_z],'off')

    papersize = 17.5*[1 (fov.y+fov.z)/(fov.y+fov.x)];
    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    print(fullfile(fileparts(nii_path),['axialcoronalsagittal_' nii_name]),'-loose','-dpdf')
    
end


