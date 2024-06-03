function coord_out = fSubsampleCoord(coord_in,index)

coord_out.x = coord_in.x(:,index);
coord_out.y = coord_in.y(:,index);
coord_out.z = coord_in.z(:,index);

if isfield(coord_in,'t')
    coord_out.t = coord_in.t;
end
