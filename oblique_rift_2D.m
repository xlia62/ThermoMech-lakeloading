
function OutputDIR=oblique_rift_2D(modelrun)

% Africa Models
InputGrid = 'InputModels/CrustalExtension_oblique2d400x220km';
InputBC = 'InputModels/Africaext_BC_oblique';


% Global variables
global water_depth dtpotential tmoho mantle_lith
global modeltime output_freq add_time max_lake_level_change
global v0 lake_level_change restart step2restart timemax
add_time = false;


switch modelrun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 0.1
        % Model 1: reference model
        output_freq = 200e3;
        restart = false; 
        timemax = 1000e3;
        InputModel = 'InputModels/AfricaExtension_input_oblique_2D';
        lake_level_change = false;
        modeltime = 15.01e6;
        v0 = 0.25; % cm half extension rate 
        water_depth = [];
        mantle_lith = 80e3;% m
        %water_depth = 10000;
        dtpotential = 0;
        tmoho = 576; %C                    
        OutputDIR = 'output/oblqiuerift/oblique_melt_test';
end

% if restart 
%     RunThermoMec2D_lake_restart(InputGrid,InputModel,InputBC,OutputDIR)
% else
RunThermoMec2D_oblique_rift(InputGrid,InputModel,InputBC,OutputDIR)
end


