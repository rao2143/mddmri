function [ rot_matrix ] = jm_rotMatrix( v )

% Computes a rotation matrix to rotate the unit vector in x-direction
% (1,0,0) onto an arbitrary vector (v1,v2,v3).

% Uses
%   jm_rotX
%   jm_rotY
%   jm_rotZ

v = v / norm(v);

if( v(2) == 0 && v(1) ~= 0 && v(3) ~= 0 )
    % in this case a single rotation around the y-axis is sufficient
    % careful: rotY turns -x onto z in a right-handed system
    eta = atan2(v(3), v(1));
    rot_matrix = jm_rotY( -eta );
elseif( v(3) == 0 && v(1) ~= 0 && v(2) ~= 0 )
    % in this case a single rotation around the z-axis is sufficient
    phi = atan2(v(2), v(1));
    rot_matrix = jm_rotZ( phi );
else
    % arbitrary vectors are first rotated into the xz-plane
    phi = atan2( norm(cross([v(1) v(2) 0],[1 0 0])), dot([v(1) v(2) 0],[1 0 0]) );
    if( v(2) > 0 )
        phi = -phi;
    end
    v_temp = jm_rotZ( phi ) * v;
    % and afterwards onto the x-axis
    eta = atan2(v_temp(3), v_temp(1));
    
    rot_matrix = jm_rotY( eta ) * jm_rotZ( phi );
    % transpose the rotation matrix for inverse rotation process
    rot_matrix = rot_matrix';
end

end

