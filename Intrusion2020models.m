function OutputDIR=Intrusion2020models(modelrun)

% Intrusion Models
InputGrid = 'InputModels/EmplacementGrid';
InputModel = 'InputModels/IntrusionModel';
InputBC = 'InputModels/AllFreeSlipBC';

% Global variables
global water_depth dtpotential tmoho
global modeltime output_freq add_time timemax

add_time = false; % yr
modeltime = 10000; % yr total model run time % was 50ka
output_freq = 4000; %yr how often you write output files
timemax = 1000; %yr maximum time step size

switch modelrun
    
    case 1
        % Model 1: no lake, ambient mantle
        
        water_depth = [];
        dtpotential = 0;  % Change of mantle potential temperature
        tmoho = 627;   % degrees C
        
        OutputDIR = 'emplacement/EMP1700_25_14';
    case 2
        % Model 1: no lake, ambient mantle
        
        water_depth = [];
        dtpotential = 0;  % Change of mantle potential temperature
        tmoho = 627;   % degrees C
        
        OutputDIR = 'emplacement/EMP1600_0_16_test';
end

RunThermoMec2D(InputGrid,InputModel,InputBC,OutputDIR)