
function OutputDIR=Africa2022models_withlake(modelrun)

% Africa Models
InputGrid = 'InputModels/CrustalExtensionFineGrid400x200km_quick';
InputBC = 'InputModels/Africaext_BC_lake2';

%InputModel = 'InputModels/AfricaExtension_newlem';
%InputModel = 'InputModels/AfricaExtension_FS';

% Global variables
global water_depth dtpotential tmoho mantle_lith
global modeltime output_freq add_time 
global v0 lake_level_change restart

add_time = false;
output_freq = 20e3;
restart = true; 

if restart 
    InputModel = 'InputModels/AfricaExtension_newlem_withlake_restart';
else
    InputModel = 'InputModels/AfricaExtension_newlem_withlake';
end

switch modelrun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 0
        % Model 1: change mantle potential temperature
        lake_level_change = false;
        modeltime = 30e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_restart_thermodynamic_reference';
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lake level 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0.1
        % Model 1: change water level to 1 km 
        lake_level_change = true;
        modeltime = 23e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_restart_90_300m';    
    
    case 0.11
        % Model 1: change water level to 1 km
        lake_level_change = false;
        modeltime = 23e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_restart_90_1000m_from_23k';
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 0.2  % extract
        % Model 1: change mantle potential temperature
        lake_level_change = true;
        modeltime =  23.10e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;
        OutputDIR = 'output/AfricaModels2022/Lake_drop_restart_110_extract';
       
    case 0.21  % extract
        % Model 1: change mantle potential temperature
        lake_level_change = false;
        modeltime =  23.10e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_restart_110_extract';    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extension rate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 0.5
        % Model 1: change mantle potential temperature
        lake_level_change = false;
        modeltime = 23e6;
        v0 = 0.5; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_restart_90_v_0.5_800m';

    case 0.51
        % Model 1: change mantle potential temperature
        lake_level_change = true;
        modeltime = 23e6;
        v0 = 0.5; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_restart_90_v_0.5_800m';     
    case 0.6
        % Model 1: change mantle potential temperature
        lake_level_change = false;
        modeltime = 23e6;
        v0 = 0.1; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_restart_90_v_0.1_800m';

    case 0.61
        % Model 1: change mantle potential temperature
        lake_level_change = true;
        modeltime = 23e6;
        v0 = 0.5; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_restart_90_v_0.1_800m';  
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dtp100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        % Model 1: change mantle potential temperature
        lake_level_change = true;
        modeltime = 21.1e6; %30e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 100;
        tmoho = 550; 
        InputModel = 'InputModels/AfricaExtension_newlem_withlake_restart_dtp100';
        OutputDIR = 'output/AfricaModels2022/Lake_drop_dtp100_restart_100_800m';        
        
    case 1.1
        % Model 1: change mantle potential temperature
        lake_level_change = false;
        modeltime = 21.1e6; %30e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 100;
        tmoho = 550;
        InputModel = 'InputModels/AfricaExtension_newlem_withlake_restart_dtp100';
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_dtp100_restart_100_800m';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lab depth 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 2
        % Model 2: change depth of LAB
        lake_level_change = true;
        modeltime = 24.1e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 70e3; % m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                       
        InputModel= 'output/AfricaModels2022/Lake_drop_restart_thermodynamic_reference_lab70';
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_lab70_restart_115_800m';
        
    case 2.1
        % Model 2: change depth of LAB
        lake_level_change = false;
        modeltime = 24.1e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 70e3; % m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                       
        InputModel = 'output/AfricaModels2022/Lake_drop_restart_thermodynamic_reference_lab70';
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_lab70_restart_115_800m'; 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rheology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 3     
        % Model 3: change rheology (wet mantle)
        lake_level_change = true;
        modeltime = 23.10e6;
        v0 = 0.25;
        water_depth = [];
        mantle_lith = 60e3;% m
        dtpotential = 0;
        tmoho = 550;
        InputModel = 'InputModels/AfricaExtension_newlem_withlake_restart_wetmantle';                
        OutputDIR = 'output/AfricaModels2022/Lake_drop_wetmantle_restart_110_800m';
        
    case 3.1     
        % Model 3: change rheology (wet mantle)
        lake_level_change = false;
        modeltime = 23.1e6;
        v0 = 0.25;
        water_depth = [];
        mantle_lith = 60e3;% m
        dtpotential = 0;
        tmoho = 550;
        InputModel = 'InputModels/AfricaExtension_newlem_withlake_restart_wetmantle';                
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_wetmantle_restart_110_800m';
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cohesion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4                      
        % Model 4: change rheology (low cohesion)     
        lake_level_change = true;
        modeltime = 23.1e6;
        v0 = 0.25;
        water_depth = [];
        mantle_lith = 60e3; % m
        dtpotential = 0;
        tmoho = 550;       
        InputModel = 'InputModels/AfricaExtension_newlem_withlake_restart_lowc';    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_lowc_restart_110_800m';
        
    case 4.1                      
        % Model 4: change rheology (low cohesion)     
        lake_level_change = false;
        modeltime = 23.1e6;
        v0 = 0.25;
        water_depth = [];
        mantle_lith = 60e3; % m
        dtpotential = 0;
        tmoho = 550;       
        InputModel = 'InputModels/AfricaExtension_newlem_withlake_restart_lowc';    
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_lowc_restart_110_800m';
end

% if restart 
%     RunThermoMec2D_lake_restart(InputGrid,InputModel,InputBC,OutputDIR)
% else
RunThermoMec2D_lake(InputGrid,InputModel,InputBC,OutputDIR)
end


