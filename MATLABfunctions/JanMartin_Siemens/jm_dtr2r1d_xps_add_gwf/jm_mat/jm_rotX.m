function [ rot_matrix ] = jm_rotX( angle_rad )

% 3D-matrix for a rotation of angle_rad around the x-axis

rot_matrix = [ [1 0                 0                ];
               [0 cos( angle_rad ) -sin( angle_rad ) ];
               [0 sin( angle_rad )  cos( angle_rad ) ] ];

end

