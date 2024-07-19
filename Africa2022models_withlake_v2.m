
function OutputDIR=Africa2022models_withlake(modelrun)

% Africa Models
InputGrid = 'InputModels/CrustalExtensionFineGrid300x250km';
%InputGrid = 'InputModels/CrustalExtensionFineGrid300x250km_quick_test';
InputBC = 'InputModels/Africaext_BC_lake_v2';

%InputModel = 'InputModels/AfricaExtension_newlem';
%InputModel = 'InputModels/AfricaExtension_FS';

% Global variables
global water_depth dtpotential tmoho mantle_lith
global modeltime output_freq add_time 
global v0 lake_level_change restart
add_time = false;
% output_freq = 200e3;
% restart = false; 

% if restart 
%     InputModel = 'InputModels/AfricaExtension_FS_restart_v2';
% else
%     InputModel = 'InputModels/AfricaExtension_FS_lake_v2';
% end
InputModel = 'InputModels/AfricaExtension_FS_lake_v2';

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
        
    case 0.1
        % Model 1: change mantle potential temperature
        output_freq = 200e3;
        restart = false; 
        InputModel = 'InputModels/AfricaExtension_FS_lake_v2';
        lake_level_change = false;
        modeltime = 16e6;
        v0 = 0.25;
        % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_model_reference_fine_test';
        
    case 0.11
        % Model 1: model with fine resolution 
        output_freq = 200e3;
        restart = false; 
        InputModel = 'InputModels/AfricaExtension_FS_lake_v2';
        lake_level_change = false;
        modeltime = 15.21e6;
        v0 = 0.25; % cm
        % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_model_reference_fine';
        
    case 0.12
        % Model 1: model with fine resolution and with changed strength
        % profile
        output_freq = 200e3;
        restart = false; 
        InputModel = 'InputModels/AfricaExtension_FS_lake_v2';
        lake_level_change = false;
        modeltime = 15.21e6;
        v0 = 0.25; % cm
        % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_model_reference_fine_strengh_profile';    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    case 1
        output_freq = 5e3;
        restart = true; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_v2';
        lake_level_change = true;
        modeltime = 13.966e6;   % 73 is at 14.4 Ma +0.3 % 70 at 13.816 +0.3 
        % 69 at 13.6166 +0.2 +0.075*2 = 13.9666 
        % for refernece fine, holes appear at 13.966 Ma
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_model_restart_69_1000m';

       
    case 1.1
        output_freq = 5e3;
        restart = true; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_v2';
        lake_level_change = false;
        modeltime = 13.966e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_model_restart_69_1000m';    
        
    case 1.2
        output_freq = 5e3;
        restart = true; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_v2';
        lake_level_change = true;
        modeltime =13.966e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_model_restart_69_600m';
        
    case 1.21
        output_freq = 5e3;
        restart = true; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_v2_300';
        lake_level_change = true;
        modeltime =13.966e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_model_restart_69_300m';   
    
    case 1.22
        output_freq = 5e3;
        restart = true; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_v2_3000';
        lake_level_change = true;
        modeltime =13.966e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_model_restart_69_3000m';   
        
	case 1.23
        output_freq = 5e3;
        restart = true; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_v2_print';
        lake_level_change = true;
        modeltime = 13.966e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_drop_model_restart_69_1000m_print';           

    case 1.24
        output_freq = 5e3;
        restart = true; 
        InputModel = 'InputModels/AfricaExtension_FS_restart_v2_print';
        lake_level_change = false;
        modeltime = 13.966e6;
        v0 = 0.25; % was 0.25
        water_depth = [];
        mantle_lith = 60e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 550;                    
        OutputDIR = 'output/AfricaModels2022/Lake_nochange_model_restart_69_1000m_print';            
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 2  % extract
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
       
    case 2.1  % extract
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
% lake level 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
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
    
    case 3.1
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
% extension rate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 4
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

    case 4.1
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
    case 4.2
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

    case 4.3
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
    case 5
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
        
    case 5.1
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

    case 6
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
        
    case 6.1
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
    case 7    
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
        
    case 7.1     
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
    case 8                      
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
        
    case 8.1                      
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
RunThermoMec2D_lake_v2(InputGrid,InputModel,InputBC,OutputDIR)
end


