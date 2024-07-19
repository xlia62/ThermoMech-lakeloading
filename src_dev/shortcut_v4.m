% S_10_frequent is the one I used in early_seis_v3.docx
% Start new numbering scheme; case 1 will be whatever S_10_frequent was
s
% What to vary
% Plate age, convergence velocity?

% Plotting of melt volumes - use Copy_of_MeltVolume.m

function OutputDIR=shortcut_v4(modelrun)

% % Global variables
% global modeltime remove_time output_freq
% global v0 mantle_lith dtpotential add_plume Texcess
% global weak_seed v1 xweak fault melting

% Global variables
global microcontinent_right plate_age_old melting ridge_loc ridge_vel
global modeltime

global dYtdt d2Yt dx tstp

switch modelrun
        
        case 1
        % S_10_frequent (no need to rerun)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 100e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sib_test';
        
        case 2
        % reduce erosion rate from 1 mm/yr to 0.5 mm/yr
        % As of 4 Myr, looks better (have eroded upper crust in uplifted
        % continent, not all of crust and part of mantl lith, as had been
        % the case in S_10_frequent)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 100e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 0.5; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_2';
        
        case 3
        % reduce topography resolution from 1 km to 500 m
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 100e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 0.5e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_3';
        
        case 4
        % reduce topography resolution from 1 km to 500 m
        % reduce erosion rate from 1 mm/yr to 0.5 mm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 100e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 0.5; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 0.5e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_4';
        
        case 5
        % reduce transport length scale from 100 km to 50 km
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 100e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 50e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_5';
        
        
        case 6
        % case 1, younger (50 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 50e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_6';
        
        case 7
        % case 1, younger (75 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 75e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_7';
        
        case 8
        % case 1, younger (80 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_8';
        
        case 9
        % case 1, younger (90 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 90e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_9';
        
        case 10
        % case 1, younger (60 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 60e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_10';
        
        case 11
        % case 1, younger (70 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 70e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_11';
        
        case 12
        % case 1, younger (85 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 85e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_12';
        
        case 13
        % case 1, younger (95 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 95e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_13';
        
        case 14
        % case 1, younger (85 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 85e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_14';
        
        case 15
        % case 1, younger (65 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 65e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_15';
        
        case 16
        % case 1, older (105 Myr) plate (note: 120 Myr scenario - plate started to subduct at the wrong end)
        % 105 Myr plate also begins to subduct
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 105e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 50e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_16';
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Vary convergence velocity from 3 cm/yr
        % Use case 8 (80 Myr plate) as standard
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 17
        % 2.8 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 60e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 2.8;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_17';
        
        case 18
        % 3.2 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 60e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3.2;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_18';
        
        case 19
        % 2.6 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 60e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 2.6;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_19';
        
        case 20
        % 3.4 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 60e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3.4;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_20';
        
        case 21
        % 3.6 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3.6;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_21';
        
        case 22
        % 3.8 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 3.8;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_22';
        
        
        case 23
        % 4 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 4;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_23';
             
        
        case 24
        % 4.2 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 4.2;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_24';
        
        
        
        case 25
        % 4.4 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 4.4;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_25';
        
        
        case 26
        % 2.4 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 2.4;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_26';
        
        case 27
        % 2.2 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 2.2;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_27';
        
        
        case 28
        % 2 cm/yr
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 2;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_28';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Try thinning Sibumasu
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 29
        % same as case 28
        % in forced_subduction.m, reduced cont_mantle_lith_thick_right from
        % 79 km to 59 km
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 2;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 1; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_29';
        
        
        case 30
        % reduce erosion rate by half from case 29
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 2;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 0.5; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_30';
        
        
        case 31
        % reduce erosion rate by half from case 30
        
        microcontinent_right = 0e3;   % how far (km) is cont from right edge of model?
        plate_age_old = 80e6;          % age of downgoing plate (years)
        melting = true;
        modeltime = 40e6;
        
        ridge_loc = 500e3;  % Location of "mid-ocean ridge"
        ridge_vel = 2;      % Half spreading rate, cm/yr
        
        % topography/erosion/sedimentation
        dYtdt = 0.25; % (mm/yr); erosion rate
        d2Yt = 10e3; % (m); max elevation
        dx = 100e3; % (m); transport length scale
        tstp = 1e3; % (m); topography resolution
        
        InputGrid = 'InputModels/M0ad_grid';
        InputModel = 'InputModels/forced_subduction';
        InputBC = 'InputModels/M6_bc';
        
        OutputDIR = 'sibumasu_31';
        

end
RunThermoMec2D(InputGrid,InputModel,InputBC,OutputDIR)

% Effect of plate age
% Needs to be older plate to get the rollback and pulse of melting
% Too old - assumption of constant plate age breaks down

% Effect of convergence velocity
% Slower convergence tends to have bigger pulse

% Erosion rate
% Not sure yet