% This routine performs numerical integration of the planets in the solar
% system and the Moon 100 years forward and compares the resulting relative
% positions w.r.t. Sun to outputs obtained from the JPL Horizons system

% REFERENCES: 
%  [1] Urban, Seidelmann - Explanatory Supplement to the Astronomical
%  Almanac, 3rd edition, University Science Books, 2013.
%  [2] Horizons System, available at 
%  https://ssd.jpl.nasa.gov/horizons/app.html#/
%  [3] Park, Folkner, Williams, Boggs - The JPL Planetary and Lunar 
%  Ephemerides DE440 and DE441, The Astronomical Journal, vol 161, No 3,
%  2021.

% Initial state of the Moon libration angles at JT_epoch 
% (phi, theta, psi, phi1, theta1, psi1) [1]:
LIB_initial = [
    0.00512995970515812456; 
    0.00004524704499022800;
    0.38239065587686011507; 
    -0.00000223092763198743;
    1.29414222411027863099;
    0.22994485870136698411
];

LIB_initial = [
    0.00512830031411853500;
    0.00004573724185991433;
    0.38239278420173690000;
   -0.00000218986174567295;
    1.29416700274878300000;
    0.22994486018992250000
];

% From DE421 obtained using the python package jplephem.
LIB_initial = [
    5.128132058714363e-03;
    1.165507165777481e-04;
    3.823932005230066e-01;
    1.461912823858170e-05;
    1.294168056057082e+00;
    2.298367282420818e-01
];

% Indices of the bodies.
ind_sun     = 1;
ind_mercury = 2;
ind_venus   = 3;
ind_earth   = 4;
ind_moon    = 5;
ind_mars    = 6;
ind_jupiter = 7;
ind_saturn  = 8;
ind_uranus  = 9;
ind_neptune = 10;
ind_pluto   = 11;

num_targets = 11;

% The intial conditions are obtained from JPL Horizons with origin at the
% body center for all bodies except the Sun. This is done since the actual 
% heliocentric initial condition in sources like [1] contain contributions 
% to the Solar System Barycenter from smaller non-planetary objects.

% Initial condition for the position (au).
R_initial = [
    4.501771046200178E-03,  7.667597654825399E-04,  2.662206318863734E-04; ...
    3.572602077668869E-01, -9.154904799747744E-02, -8.598103172768061E-02; ...
    6.082494317526214E-01, -3.491324458403040E-01, -1.955443448720758E-01; ...
    1.160247250568800E-01, -9.265813150908608E-01, -4.017930739420348E-01; ...
    1.152165477160266E-01, -9.285759450811601E-01, -4.028803366239009E-01; ...
   -1.146885853705725E-01, -1.328366526279570E+00, -6.061551990541251E-01;...
   -5.384209277643075E+00, -8.312483870146020E-01, -2.250951187017292E-01; ...
    7.889888161590055E+00,  4.595710989682622E+00,  1.558429793513401E+00; ...
   -1.826990605559244E+01, -1.162723764514111E+00, -2.503714950236171E-01; ...
   -1.605954001605754E+01, -2.394295935145348E+01, -9.400423444456409E+00; ...
   -3.048781548042960E+01, -8.731761132527535E-01,  8.911305385348918E+00
];

% Initial condition for the velocity (au/d).
V_initial = [
   -3.517482063742068E-07,  5.177625639825057E-06,  2.229090837927562E-06; ...
    3.367845709043980E-03,  2.488934292987171E-02,  1.294407129214571E-02; ...
    1.095242018625395E-02,  1.561250662984909E-02,  6.328876605810963E-03; ...
    1.680431652771057E-02,  1.745166161927448E-03,  7.570134278727792E-04; ...
    1.740540134429601E-02,  1.577720694761970E-03,  6.714512881328027E-04; ...
    1.448200480836478E-02,  2.372854532106865E-04, -2.837498250145917E-04; ...
    1.092364404049427E-03, -6.523294106045839E-03, -2.823012134541729E-03; ...
   -3.217204775856843E-03,  4.330632712157481E-03,  1.926417212037553E-03; ...
    2.215425014444429E-04, -3.767652400674729E-03, -1.653244046240149E-03; ...
    2.643121823195712E-03, -1.503490013808834E-03, -6.812710872439278E-04; ...
    3.225591308955338E-04, -3.148753752998151E-03, -1.080178675229524E-03
];

% Astronomical unit (km).
au = 149597870.691;

% Gauss' (gravitational) constant from Table 8.4 in [1].
k = 0.01720209895;

% Earth-Moon mass ratio from Table 8.4 in [1].
em_ratio = 81.30056;

% Standard gravitational parameters from Table 8.4 in [1].
mu = k*k ./ [
    1;
    6023600;
    408523.71;
    328900.5614 * (em_ratio + 1) / em_ratio;
    328900.5614 * (em_ratio + 1);
    3098708;
    1047.3486;
    3497.898;
    22902.98;
    19412.24;
    135200000
]';
%mu(5) =  1.093189443939597e-11
mu(5) = 1.093189450705846e-11;

% From Table 2 in [3].
GM_sun = 132712440041.279419; 
GM_moon= 4902.800118; 
GM_mer = 22031.868551;
GM_ven = 324858.592000; 
GM_ear = 398600.435507;
GM_mar = 42828.375816;
GM_jup = 126712764.100000;
GM_sat = 37940584.841800;
GM_ura = 5794556.400000;
GM_nep = 6836527.100580;
GM_plu = 975.500000;

mu(2) = k * k * GM_mer / GM_sun;
mu(3) = k * k * GM_ven / GM_sun;
mu(4) = k * k * GM_ear / GM_sun;
mu(5) = k * k * GM_moon / GM_sun;
mu(6) = k * k * GM_mar / GM_sun;
mu(7) = k * k * GM_jup / GM_sun;
mu(8) = k * k * GM_sat / GM_sun;
mu(9) = k * k * GM_ura / GM_sun;
mu(10)= k * k * GM_nep / GM_sun;
mu(11)= k * k * GM_plu / GM_sun;

% Load data for checking the accuracy of the results.
jpl_reference;

% Make coordinates consistent with the Sun expressed with respect to SSB by 
% adding it to all OSVs.
for ind_target = 2:num_targets    
    R_initial(ind_target, :) = R_initial(ind_target, :) + R_initial(1, :);
    V_initial(ind_target, :) = V_initial(ind_target, :) + V_initial(1, :);
end

% OSVs related to the DoFs during the integration.
R = R_initial;
V = V_initial;

% Compute relativistic barycenter.
[r_bary, v_bary] = barycenter(R, V, mu, true);

% Offset bodies w.r.t. so that the relativistic barycenter is approximately
% at the origin.
for ind_target = 1:num_targets
    R(ind_target, :) = R(ind_target, :) - r_bary;
    V(ind_target, :) = V(ind_target, :) - v_bary;
end

reference = {
    [0, 0, 0, 0, 0];
    mercury_exp;
    venus_exp;
    earth_exp;
    moon_exp;
    mars_exp;
    jupiter_exp;
    saturn_exp;
    uranus_exp;
    neptune_exp;
    pluto_exp
};
names = ["Sun", "Mercury", "Venus", "Earth", "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"];
names_short = ["S", "M", "V", "E", "Mo", "Ma", "J", "S", "U", "N", "P"];


% Matrix with position errors in (km) for all targets for each timestep.
ERR = [];
% Time step size in days.
JT_timestep = 0.05;
% Start date 1969-Jun-28 00:00:00.0000
JT_epoch = 2440400.50;
% Compute number of timesteps.
num_timesteps = floor(365.25 * 100 / JT_timestep + 1);
% Current Julian time.
JT = JT_epoch;
% Time after epoch.
t = 0;


[R, V] = moon_frombary(R, V, mu, 4, 5);

% We first have 6 DoFs for the libration angles and their time derivatives
% and then additional 6 DoFs (x, y, z, x1, y1, z1) for each body.
y = [LIB_initial; osv_to_dof(R, V)];
[R, V] = moon_tobary(R, V, mu, 4, 5);

y_start2 = y;
R_start2 = R;
V_start2 = V;


ERR = [];
ERR_angle = [];
LIB_out = [];
F = [];
JT = t + JT_epoch;

for timestep = 1:num_timesteps
    % Prepare output print to display during the integration.
    err_step = zeros(1, num_targets);
    err_angle= zeros(1, num_targets);
    err_str = "";
    err_str_mas = "";
    filled = false;

    for ind_target = 1:num_targets
        ref_target = reference{ind_target};
        ref_earth  = reference{ind_earth};

        results = find(abs(ref_target(:, 1) - JT) < 1/86400);
        if length(results) > 0
            year         = ref_target(results, 2);
            if ~filled
                err_str = sprintf("%d - %.10f ", year, JT);
            end

            r_target     = R(ind_target, 1:3) - R(1, 1:3);
            r_earth      = R(ind_earth, 1:3) - R(1, 1:3);

            if ind_target == 4
                vector_target =r_earth / norm(r_earth);
            else
                vector_target = (r_target - r_earth) / norm(r_target - r_earth);;
            end

            r_target_exp = ref_target(results, 3:5);
            r_earth_exp  = ref_earth(results, 3:5);
            if ind_target == 4
                vector_exp =  (r_earth_exp) / norm(r_earth_exp);
            else
                vector_exp =  (r_target_exp - r_earth) / norm(r_target_exp - r_earth);
            end
            err_target = norm(r_target - r_target_exp) * au;
            %err_angle_mas = acosd(min(1, abs(dot(vector_target, vector_exp)))) * 3600000
            err_angle_mas = vectors_angle(vector_target, vector_exp) * 3600000;


            err_step(1, ind_target) = err_target;
            err_angle(1, ind_target) = err_angle_mas;

            err_str_target     = sprintf("%s %.02f ", names_short(ind_target), err_target);
            err_str_target_mas = sprintf("%s %.01f ", names_short(ind_target), err_angle_mas);
            err_str = strcat(err_str, err_str_target);
            err_str_mas = strcat(err_str_mas, err_str_target_mas);
            filled = true;
        end        
    end
    if filled
        disp(err_str);
        disp(err_str_mas);
        ERR = [ERR; year, err_step];
        ERR_angle = [ERR_angle; year, err_angle];
        LIB_out = [LIB_out; y(1:6)'];
    end


    % Perform the actual integration.
    if timestep <= 8 
        [y, t] = runge4(@func_lib, t, y, JT_timestep, mu);
        f = func_lib(t, y, mu);
        F = [f, F];
    else 
        [y, t, F] = adams8(@func_lib, t, y, F, JT_timestep, mu);
    end
    % Convert DoFs back to OSVs.
    [R, V] = dof_to_osv(y(7:end));
    [R, V] = moon_tobary(R, V, mu, 4, 5);

    JT = t + JT_epoch;

    if filled
        % Visualize the error.
        figure(1);
        clf
        subplot(1, 2, 1);
        
        year        = ERR(:, 1);
        err_mercury = ERR(:, 3);
        err_venus   = ERR(:, 4);
        err_earth   = ERR(:, 5);
        err_moon    = ERR(:, 6);
        err_mars    = ERR(:, 7);
        err_jupiter = ERR(:, 8);
        err_saturn  = ERR(:, 9);
        err_uranus  = ERR(:, 10);
        err_neptune = ERR(:, 11);
        err_pluto   = ERR(:, 12);
        
        endyear = 1969 + length(ERR) - 1;
        loglog(year-1969, err_mercury, 'r',   'LineWidth', 2)
        hold on;
        loglog(year-1969, err_venus,   'g',   'LineWidth', 2)
        loglog(year-1969, err_earth,   'b',   'LineWidth', 2)
        loglog(year-1969, err_moon,    'bo',  'LineWidth', 2)
        loglog(year-1969, err_mars,    'm',   'LineWidth', 2)
        loglog(year-1969, err_jupiter, 'c',   'LineWidth', 2)
        loglog(year-1969, err_saturn,  'y',   'LineWidth', 2)
        loglog(year-1969, err_uranus,  'k',   'LineWidth', 2)
        loglog(year-1969, err_neptune, 'k--', 'LineWidth', 2)
        loglog(year-1969, err_pluto,   'k:',  'LineWidth', 2)
        
        legend('Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', ...
            'Neptune', 'Pluto', 'FontSize', 18, 'Location', 'NorthWest');
        xlabel('Year', 'FontSize', 18);
        ylabel('Error in Position (km)', 'FontSize', 18);
        %xlim([1969, endyear])
        xlim([1 100]);
        ylim([0.01 1000])
        grid on
        title(sprintf('Numerical Integration Error %d-%d', min(year), max(year)), 'FontSize', 18);

        subplot(1, 2, 2);
        
        year        = ERR(:, 1);
        err_mercury = ERR_angle(:, 3);
        err_venus   = ERR_angle(:, 4);
        err_sun     = ERR_angle(:, 5);
        err_moon    = ERR_angle(:, 6);
        err_mars    = ERR_angle(:, 7);
        err_jupiter = ERR_angle(:, 8);
        err_saturn  = ERR_angle(:, 9);
        err_uranus  = ERR_angle(:, 10);
        err_neptune = ERR_angle(:, 11);
        err_pluto   = ERR_angle(:, 12);
        
        endyear = 1969 + length(ERR) - 1;
        loglog(year-1969, err_mercury, 'r',   'LineWidth', 2)
        hold on;
        loglog(year-1969, err_venus,   'g',   'LineWidth', 2)
        loglog(year-1969, err_sun,     'b',   'LineWidth', 2)
        loglog(year-1969, err_moon,    'bo',  'LineWidth', 2)
        loglog(year-1969, err_mars,    'm',   'LineWidth', 2)
        loglog(year-1969, err_jupiter, 'c',   'LineWidth', 2)
        loglog(year-1969, err_saturn,  'y',   'LineWidth', 2)
        loglog(year-1969, err_uranus,  'k',   'LineWidth', 2)
        loglog(year-1969, err_neptune, 'k--', 'LineWidth', 2)
        loglog(year-1969, err_pluto,   'k:',  'LineWidth', 2)
        
        legend('Mercury', 'Venus', 'Sun', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', ...
            'Neptune', 'Pluto', 'FontSize', 18, 'Location', 'NorthWest');
        xlabel('Year', 'FontSize', 18);
        ylabel('Error in Angle from Earth (mas)', 'FontSize', 18);
        %xlim([1969, endyear])
        xlim([1 100]);
        ylim([0.01 3600])
        grid on
        title(sprintf('Numerical Integration Error %d-%d', min(year), max(year)), 'FontSize', 18);
    end
end