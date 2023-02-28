function OutputDIR=Africa2019models(modelrun)

% Africa Models
InputGrid = 'InputModels/CrustalExtensionFineGrid400x200km';
InputBC = 'InputModels/Africaext_BC';

% Global variables
global water_depth dtpotential tmoho
global modeltime output_freq add_time

add_time = false;
modeltime = 16e6;
output_freq = 200e3;

switch modelrun
    
    case 1
        % Model 1: no lake, ambient mantle
        InputModel = 'InputModels/AfricaExtension_noerosion';

        water_depth = [];
        dtpotential = 0;
        tmoho = 600;
                        
        OutputDIR = 'AfricaModels2021/noerosion';
        
    case 2
        % Model 2: no lake, warmer mantle 150K
        InputModel = 'InputModels/AfricaExtension_diffusion';

        water_depth = [];
        dtpotential = 0;
        tmoho = 600;
                        
        OutputDIR = 'AfricaModels2021/diffusion';
        
    case 3
        % Model 2: no lake, warmer mantle 150K
        InputModel = 'InputModels/AfricaExtension_newlem2';

        water_depth = [];
        dtpotential = 0;
        tmoho = 550;
                        
        OutputDIR = 'AfricaModels2021/newlem4';
        
    case 4
        % Model 4: dee lake, warmer mantle 150K
        
        water_depth = 1200;
        dtpotential = 150;
        tmoho = 500;
                        
        OutputDIR = '/data1/adallenl/emplacement/AfricaExt4';

    case 5
        % Model 5: shallow lake, warmer mantle 150K
        
        water_depth = 400;
        dtpotential = 0;
        tmoho = 600;
                        
        OutputDIR = '/data1/adallenl/emplacement/AfricaExt5';

    case 6
        % Model 5: shallow lake, warmer mantle 150K
        
        water_depth = 400;
        dtpotential = 150;
        tmoho = 600;
                        
        OutputDIR = '/data1/adallenl/emplacement/AfricaExt6';
        
    case 7
        % Model 1: no lake, ambient mantle
        
        water_depth = [];
        dtpotential = 0;
        tmoho = 600;
                        
        OutputDIR = '/data1/adallenl/emplacement/AfricaExt7';

        
end

RunThermoMec2D(InputGrid,InputModel,InputBC,OutputDIR)