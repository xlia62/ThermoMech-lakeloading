
function OutputDIR=Africa2021models(modelrun)

% Africa Models
InputGrid = 'InputModels/CrustalExtensionFineGrid400x200km';
InputBC = 'InputModels/Africaext_BC';
InputModel = 'InputModels/AfricaExtension_newlem2';
%InputModel = 'InputModels/AfricaExtension_FS';

% Global variables
global water_depth dtpotential tmoho
global modeltime output_freq add_time
global v0

add_time = false;
output_freq = 200e3;

switch modelrun
    
    case 1
        % Model 1: AfricaMaterials_layercorr
        
        modeltime = 5.1e6;
        v0 = 0.5;
        water_depth = [];
        dtpotential = 0;
        tmoho = 550;
                        
        OutputDIR = 'output/AfricaModels2021/fastscape_line2';

    case 2
        % Model 1: AfricaMaterials_layercorr
        modeltime = 5e6;
        v0 = 0.25;
        water_depth = [];
        dtpotential = 0;
        tmoho = 500;
                        
        OutputDIR = 'output/AfricaModels2021/noerosion';

    case 3
               
        % Model 1: AfricaMaterials_layercorr
        modeltime = 12e6;
        v0 = 0.5;
        water_depth = [];
        dtpotential = 100;
        tmoho = 600;
                        
        OutputDIR = 'output/AfricaModels2021/fastT600';
        
    case 4
                       
        % Model 1: AfricaMaterials_layercorr
        
        modeltime = 7e6;
        v0 = 0.5;
        water_depth = [];
        dtpotential = 200;
        tmoho = 600;
                        
        OutputDIR = 'output/AfricaModels2021/fastT600hot';
        
end

RunThermoMec2D(InputGrid,InputModel,InputBC,OutputDIR)
