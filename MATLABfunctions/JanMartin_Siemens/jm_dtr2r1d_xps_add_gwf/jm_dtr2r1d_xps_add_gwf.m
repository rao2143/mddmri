function xps_gwf = jm_dtr2r1d_xps_add_gwf(fn_xps, fn_gwf, bSave)

% Settings
gradAmp = 80e-3;  % Maximum gradient amplitude available on the scanner
delta_t = 1e-5;  % Stepsize between time points
rps.maxomega = 100 * 2 * pi;  % Maximum frequency
rps.nomega = 100;  % Number of points in frequency space

% Load original xps
xps_gwf = mdm_xps_load(fn_xps);

% Load workspace with gradient waveforms
load(fn_gwf);

xps_gwf.gwf_t = repmat(gwf_t, xps_gwf.n, 1);
xps_gwf.gwf_x = zeros(xps_gwf.n, size(gwf_x_lin, 2));
xps_gwf.gwf_y = xps_gwf.gwf_x;
xps_gwf.gwf_z = xps_gwf.gwf_x;

tol = 0.01;  % Tolerance for b_delta
for ii = 1 : 1 : xps_gwf.n
    orientation = xps_gwf.u(ii, :)';
    orientation = orientation([2 1 3]);  % rotation matrix is calculated in PRS instead of RPS coordinates
    b_delta = xps_gwf.b_delta(ii);
    % Case 1: Linear
    if abs(b_delta - 1) < tol
        % Case 1.1: Short TE
        if xps_gwf.te(ii) < 0.06
            gwf_x = gwf_x_lin_short;
            gwf_y = gwf_y_lin_short;
            gwf_z = gwf_z_lin_short;
        % Case 1.2: Long TE
        else
            gwf_x = gwf_x_lin;
            gwf_y = gwf_y_lin;
            gwf_z = gwf_z_lin;
        end
    % Case 3: Spherical
    elseif abs(b_delta) < tol
        gwf_x = gwf_x_sph;
        gwf_y = gwf_y_sph;
        gwf_z = gwf_z_sph;
    % Case 2: Planar
    elseif abs(b_delta + 0.5) < tol
        gwf_x = gwf_x_pln;
        gwf_y = gwf_y_pln;
        gwf_z = gwf_z_pln;
    else
        error(['Gradient waveform for b_delta at position ' num2str(ii) ' not supplied.']);
    end
    % Scale gradient waveform to b-value
    [gwf_x_scaled, gwf_y_scaled, gwf_z_scaled] = scaleGradientWaveform(xps_gwf.b(ii), gradAmp, delta_t, gwf_x, gwf_y, gwf_z);    
    % Rotate gradient waveform
    [gwf_x_rot, gwf_y_rot, gwf_z_rot] = rotateGradientWaveform(orientation, gwf_x_scaled, gwf_y_scaled, gwf_z_scaled);
    % Scale to scanner gradients
    xps_gwf.gwf_x(ii, :) = gwf_x_rot * gradAmp;
    xps_gwf.gwf_y(ii, :) = gwf_y_rot * gradAmp;
    xps_gwf.gwf_z(ii, :) = gwf_z_rot * gradAmp;
end

% Calculated btomega
xps_gwf = mdm_xps_add_btomega(xps_gwf, rps);

% Save updated xps
if bSave
    [pn, fn, ext] = msf_fileparts(fn_xps);
    fn_save = fullfile(pn, [fn '_gwf' ext]);
    mdm_xps_save(xps_gwf, fn_save);
end

end


function [gwf_x_rot, gwf_y_rot, gwf_z_rot] = rotateGradientWaveform(v, gwf_x, gwf_y, gwf_z)

rot_matrix = jm_rotMatrix(v);
gwf = cat(1, gwf_y, gwf_x, gwf_z);  % gwf in Daniel's coordinates is RPS, gwf in Scanner's coordinates is PRS

gwf_rot = rot_matrix * gwf;
gwf_rot(abs(gwf_rot) < 1e-6) = 0;

gwf_x_rot = gwf_rot(2, :);
gwf_y_rot = gwf_rot(1, :);
gwf_z_rot = gwf_rot(3, :);

end


function [gwf_x_scaled, gwf_y_scaled, gwf_z_scaled] = scaleGradientWaveform(b, gradAmp, delta_t, gwf_x, gwf_y, gwf_z)

gyro_ratio = 2.6752218744e8;

Q = [cumsum(gwf_y); cumsum(gwf_x); cumsum(gwf_z)] * delta_t * gradAmp * gyro_ratio;  % q-tensor
mat_int = jm_mat_int(size(Q, 2), delta_t);  % integration matrix
B = Q * mat_int * Q';  % b-tensor
max_b = trace(B);  % maximum possible b-value in SI units

if b > max_b
    error('b-value cannot exceed maximum possible b-value.')
end

gwf_x_scaled = gwf_x * sqrt(b / max_b);
gwf_y_scaled = gwf_y * sqrt(b / max_b);
gwf_z_scaled = gwf_z * sqrt(b / max_b);

end
