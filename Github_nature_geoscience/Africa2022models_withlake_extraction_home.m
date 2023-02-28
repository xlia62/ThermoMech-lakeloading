
function OutputDIR=Africa2022models_withlake_extraction_home(modelrun)

% Africa Models
InputGrid = 'InputModels/CrustalExtensionFineGrid300x250km';

InputBC = 'InputModels/Africaext_BC_lake_v2';

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
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restart 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    case 1.1
        output_freq = 5e3;
        restart = true; 
        step2restart = 71;
        max_lake_level_change = 800; 
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
        OutputDIR = 'output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_dry';

end

% if restart 
%     RunThermoMec2D_lake_restart(InputGrid,InputModel,InputBC,OutputDIR)
% else
RunThermoMec2D_lake_extraction(InputGrid,InputModel,InputBC,OutputDIR)
end


