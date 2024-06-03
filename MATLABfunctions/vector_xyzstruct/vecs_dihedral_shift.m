function delta = vecs_dihedral_shift(r1,r1x,r1y,r2,r2x,r2y,r3,r3x,r3y,r4,r4x,r4y)% r11x = vecs_dif(r1,r1x);% r11y = vecs_dif(r1,r1y);r11x = vecs_dif(r1x,r1);r11y = vecs_dif(r1y,r1);r1x2x = vecs_dif(r1x,r2x);r1x2y = vecs_dif(r1x,r2y);r1y2x = vecs_dif(r1y,r2x);r1y2y = vecs_dif(r1y,r2y);r1x3x = vecs_dif(r1x,r3x);r1x3y = vecs_dif(r1x,r3y);r1y3x = vecs_dif(r1y,r3x);r1y3y = vecs_dif(r1y,r3y);r1x4x = vecs_dif(r1x,r4x);r1x4y = vecs_dif(r1x,r4y);r1y4x = vecs_dif(r1y,r4x);r1y4y = vecs_dif(r1y,r4y);theta1x4x = vecs_angle(r11x,r1x4x);theta1x4y = vecs_angle(r11x,r1x4y);theta1y4x = vecs_angle(r11y,r1y4x);theta1y4y = vecs_angle(r11y,r1y4y);theta1x3x = vecs_angle(r11x,r1x3x);theta1x3y = vecs_angle(r11x,r1x3y);theta1y3x = vecs_angle(r11y,r1y3x);theta1y3y = vecs_angle(r11y,r1y3y);theta1x2x = vecs_angle(r11x,r1x2x);theta1x2y = vecs_angle(r11x,r1x2y);theta1y2x = vecs_angle(r11y,r1y2x);theta1y2y = vecs_angle(r11y,r1y2y);delta_1x2x = geom2delta(theta1x4x,r1x2x);delta_1x2y = geom2delta(theta1x4y,r1x2y);delta_1y2x = geom2delta(theta1y4x,r1y2x);delta_1y2y = geom2delta(theta1y4y,r1y2y);delta_1x3x = geom2delta(theta1x3x,r1x3x);delta_1x3y = geom2delta(theta1x3y,r1x3y);delta_1y3x = geom2delta(theta1y3x,r1y3x);delta_1y3y = geom2delta(theta1y3y,r1y3y);delta_1x4x = geom2delta(theta1x4x,r1x4x);delta_1x4y = geom2delta(theta1x4y,r1x4y);delta_1y4x = geom2delta(theta1y4x,r1y4x);delta_1y4y = geom2delta(theta1y4y,r1y4y);% delta_1x4x = 1680*cos(theta1x4x).*exp(-2.671*vecs_length(r1x4x)*10);% delta_1x4y = 1680*cos(theta1x4y).*exp(-2.671*vecs_length(r1x4y)*10);% delta_1y4x = 1680*cos(theta1y4x).*exp(-2.671*vecs_length(r1y4x)*10);% delta_1y4y = 1680*cos(theta1y4y).*exp(-2.671*vecs_length(r1y4y)*10);% delta_1x3x = 1680*cos(theta1x4x).*exp(-2.671*vecs_length(r1x3x)*10);% delta_1x3y = 1680*cos(theta1x4y).*exp(-2.671*vecs_length(r1x3y)*10);% delta_1y3x = 1680*cos(theta1y4x).*exp(-2.671*vecs_length(r1y3x)*10);% delta_1y3y = 1680*cos(theta1y4y).*exp(-2.671*vecs_length(r1y3y)*10);% delta_1x2x = 1680*cos(theta1x4x).*exp(-2.671*vecs_length(r1x2x)*10);% delta_1x2y = 1680*cos(theta1x4y).*exp(-2.671*vecs_length(r1x2y)*10);% delta_1y2x = 1680*cos(theta1y4x).*exp(-2.671*vecs_length(r1y2x)*10);% delta_1y2y = 1680*cos(theta1y4y).*exp(-2.671*vecs_length(r1y2y)*10);delta = delta_1x4x + delta_1x4y + delta_1y4x + delta_1y4y + delta_1x3x + delta_1x3y + delta_1y3x + delta_1y3y + delta_1x2x + delta_1x2y + delta_1y2x + delta_1y2y;end    function delta = geom2delta(theta_chh,r_hh)        delta = 1680*cos(theta_chh).*exp(-2.671*vecs_length(r_hh)*10);    end