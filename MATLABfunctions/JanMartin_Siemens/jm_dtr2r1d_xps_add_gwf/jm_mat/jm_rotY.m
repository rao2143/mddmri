function [ rot_matrix ] = jm_rotY( angle_rad )

% 3D-matrix for a rotation of angle_rad around the y-axis

rot_matrix = [ [ cos( angle_rad ) 0 sin( angle_rad )];
               [ 0                1 0               ];
               [-sin( angle_rad ) 0 cos( angle_rad )] ];

end

