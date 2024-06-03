function [ rot_matrix ] = jm_rotZ( angle_rad )

% 3D-matrix for a rotation of angle_rad around the z-axis

rot_matrix = [ [cos( angle_rad ) -sin( angle_rad ) 0 ];
               [sin( angle_rad )  cos( angle_rad ) 0 ];
               [0                 0                1 ] ];

end


