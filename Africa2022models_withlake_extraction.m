
function OutputDIR=Africa2022models_withlake(modelrun)

% Africa Models
InputGrid = 'InputModels/CrustalExtensionFineGrid300x250km';
%InputGrid = 'InputModels/CrustalExtensionFineGrid300x250km_fine';
InputBC = 'InputModels/Africaext_BC_lake_v2';

%InputModel = 'InputModels/AfricaExtension_newlem';
%InputModel = 'InputModels/AfricaExtension_FS';

% Global variables
global water_depth dtpotential tmoho mantle_lith
global modeltime output_freq add_time max_lake_level_change
global v0 lake_level_change restart step2restart timemax
add_time = false;
% output_freq = 200e3;
% restart = false; 

% if restart 
%     InputModel = 'InputModels/AfricaExtension_FS_restart_v2';
% else
%     InputModel = 'InputModels/AfricaExtension_FS_lake_v2';
% end

switch modelrun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 0.1
        % Model 1: reference model
        output_freq = 200e3;
        restart = false; 
        timemax = 1000e3;
        InputModel = 'InputModels/AfricaExtension_FS_lake_extraction_dry';
        lake_level_change = false;
        modeltime = 15.01e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/reference_71_dry_test';
        
    case 0.2
        % Model: change mantle potential temperature
        output_freq = 200e3;
        restart = false;
        timemax = 1000e3;
        InputModel = 'InputModels/AfricaExtension_FS_lake_extraction_dry';
        lake_level_change = false;
        modeltime = 15.01e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 10;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/reference_Tp10_dry';
        
    case 0.3
        % Model: change mantle material 
        output_freq = 200e3;
        restart = false; 
        InputModel = 'InputModels/AfricaExtension_FS_lake_extraction_wetmantle';
        lake_level_change = false;
        modeltime = 15e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/reference_wetmantle';
        
    case 0.4
        % Model: LAB depth 
        output_freq = 200e3;
        restart = false; 
        InputModel = 'InputModels/AfricaExtension_FS_lake_extraction_dry';
        lake_level_change = false;
        modeltime = 15.01e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 65e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/reference_LAB65_dry';

    case 0.5 
        % Model 1: reference model with no lake, just air 
        output_freq = 200e3;
        restart = false; 
        InputModel = 'InputModels/AfricaExtension_FS_air_extraction';
        lake_level_change = false;
        modeltime = 17.1e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/reference_air';
        
        
    case 0.6
        % Model 1: reference model with weak curst by a factor of 0.1 and
        % the weak zoon
        output_freq = 200e3;
        restart = false;
        timemax = 1000e3;
        InputModel = 'InputModels/AfricaExtension_FS_lake_extraction_weakcrust';
        lake_level_change = false;
        modeltime = 15.01e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/reference_weakcrust_dry';
        
        
    case 0.7
        % Model 1: reference model with weak curst by a factor of 10 also
        % incluidng the wewak zone 
        output_freq = 200e3;
        timemax = 1000e3;
        restart = false; 
        InputModel = 'InputModels/AfricaExtension_FS_lake_extraction_strongcrust';
        lake_level_change = false;
        modeltime = 15.01e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/reference_strongcrust_dry';
        
    case 0.8
        % Model 1: reference model with different velocity 
        output_freq = 200e3;
        timemax = 1000e3;
        restart = false; 
        InputModel = 'InputModels/AfricaExtension_FS_lake_extraction_dry';
        lake_level_change = false;
        modeltime = 15.01e6;
        v0 = 0.275; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/reference_v27_dry';
        
    case 0.9
        % Model 1: reference model
        output_freq = 200e3;
        restart = false; 
        timemax = 1000e3;
        InputModel = 'InputModels/AfricaExtension_FS_lake_extraction_dry';
        lake_level_change = false;
        modeltime = 17.01e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/reference_71_dry';
        
 % review: high resoluiton 
    case 0.91
        % Model 1: reference model
        output_freq = 200e3;
        restart = false; 
        timemax = 1000e3;
        InputModel = 'InputModels/AfricaExtension_FS_lake_extraction_dry';
        lake_level_change = false;
        modeltime = 15.01e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/reference_71_dry_highres';
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    case 1.1
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_dry_test';

   
        
    case 1.2
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        max_lake_level_change = 1000;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = false;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_71_dry';
 
    case 1.11
        output_freq = 5e3;
        restart = true; 
        step2restart = 81;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 16.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_81_1000m';

   
        
    case 1.21
        output_freq = 5e3;
        restart = true; 
        step2restart = 81;
        max_lake_level_change = 1000;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = false;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 16.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_81_1000m';        
        
        
    case 1.3
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        max_lake_level_change = -1000;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_1000m_rev';
        
    case 1.4
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = false;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_71_1000m_rev';
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart:lake level 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2.1  % lake level change of 600 m
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        max_lake_level_change = 600; 
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_600m_test';  
    
    case 2.2  % lake level change of 300 m

        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        max_lake_level_change = 300; 
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_300m';
        
    case 2.3  % lake level change of 600 m
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        max_lake_level_change = 1500; 
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_1500m';  
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart: extension rate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 3.1
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.5; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_1000m_halfexten_0.5';

    case 3.2
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.1; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_1000m_halfexten_0.1';

    case 3.3
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = false;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.5; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_71_1000m_halfexten_0.5';


    case 3.4
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = false;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.1; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_71_1000m_halfexten_0.1';        

    case 3.5
        output_freq = 5e3;
        timemax = 5e3;
        restart = true; 
        step2restart = 71;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_v27';
        lake_level_change = true;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.275; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_v27';        
        
    case 3.6
        output_freq = 5e3;
        timemax = 5e3;
        restart = true; 
        step2restart = 71;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_v27';
        lake_level_change = false;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.275; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_71_v27';        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart LAB depth 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4.1
        output_freq = 5e3;
        restart = true; 
        step2restart = 71; %66;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_LAB65';
        lake_level_change = true;
        % step71:13.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6; %13.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 65e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_LAB65';  
        
    case 4.2
        output_freq = 5e3;
        restart = true; 
        step2restart = 71; % 66;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_LAB65';
        lake_level_change = false;
        % step71:13.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6; %13.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 65e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_71_LAB65';  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mantle potential temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 5.1
        output_freq = 5e3;
        restart = true; 
        timemax = 5e3;
        step2restart = 71; %44;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_Tp30';
        lake_level_change = true;
        % step71:14.0419 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6; %8.95e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 10;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_1000m_Tp30';  
        
    case 5.2
        output_freq = 5e3;
        restart = true; 
        timemax = 5e3;
        step2restart = 71; %44;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_Tp30';
        lake_level_change = false;
        % step71:13.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6; %8.95e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 10;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_71_1000m_Tp30';  

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rheology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 6.1
        output_freq = 5e3;
        restart = true; 
        step2restart = 38;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_wetmantle';
        lake_level_change = true;
        % step38:7.403 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 7.75e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_38_1000m_wetmantle';  
        
    case 6.2
        output_freq = 5e3;
        restart = true; 
        step2restart = 38;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_wetmantle';
        lake_level_change = false;
        % step71:13.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 7.75e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_38_1000m_wetmantle';  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart weak or strong crust 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    case 7.1
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        timemax = 5e3;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_weakcrust';
        lake_level_change = true;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_weakcrust2';

   
        
    case 7.2
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        timemax = 5e3;
        max_lake_level_change = 1000;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_weakcrust';
        lake_level_change = false;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_71_weakcrust2';
        
     case 7.3
        output_freq = 5e3;
        restart = true; 
        timemax = 5e3;
        step2restart = 71;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_strongcrust';
        lake_level_change = true;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_strongcrust2';
     
    case 7.4
        output_freq = 5e3;
        restart = true;
        timemax = 5e3;
        step2restart = 71;
        max_lake_level_change = 1000;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_strongcrust';
        lake_level_change = false;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_71_strongcrust2';       
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart:lake level change pattern 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 8.1  % lake level change of one period
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        timemax = 5e3;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        max_lake_level_change = 1000; 
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_oneperiod';          
        
        
    case 8.2  % lake level change of one period
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        timemax = 5e3;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        max_lake_level_change = 1000; 
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_justdrop';          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart: high resolution (150 m) restart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 9.1  % lake level change of 1000 m with high resolution
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_highres';
        lake_level_change = true;
        max_lake_level_change = 1000; 
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_1000m_highres_test';
        
    case 9.2  % lake level no change but high resolution 
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_highres';
        lake_level_change = false;
        max_lake_level_change = 1000; 
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.1*2+0.075*2 =0.35
        modeltime = 14.35e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_nochange_restart_71_1000m_highre_test';  
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart with shorter period 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    case 10.1
        output_freq = 5e3;
        restart = true; 
        step2restart = 72;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_40kyr';
        lake_level_change = true;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.04*2+0.075*2 =0.23
        modeltime = 14.43e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_dry_40kyr_test';
        % the reference model would be the same: Lake_nochange_restart_71_1000m_dry       
   case 10.2
        output_freq = 5e3;
        restart = true; 
        step2restart = 72;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction_40kyr';
        lake_level_change = false;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.04*2+0.075*2 =0.23
        modeltime = 14.43e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_nochange_71_dry_40kyr';
        % the reference model would be the same: Lake_nochange_restart_71_1000m_dry       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart with 90% sediments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    case 11.1
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = true;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.04*2+0.075*2 =0.23
        modeltime = 14.4e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_dry_sed';
        % the reference model would be the same: Lake_nochange_restart_71_1000m_dry       
   case 11.2
        output_freq = 5e3;
        restart = true; 
        step2restart = 71
        max_lake_level_change = 1000; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_extraction';
        lake_level_change = false;
        % step71:14.0024 Ma 
        % period = 0.1Ma, buffer= 0.075Ma
        % total time = 0.04*2+0.075*2 =0.23
        modeltime = 14.4e6;   
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_nochange_71_dry_sed';
        % the reference model would be the same: Lake_nochange_restart_71_1000m_dry 
end

% if restart 
%     RunThermoMec2D_lake_restart(InputGrid,InputModel,InputBC,OutputDIR)
% else
RunThermoMec2D_lake_extraction(InputGrid,InputModel,InputBC,OutputDIR)
end


