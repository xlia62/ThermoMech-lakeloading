function [XMELT,HLAT] = Melt_fraction(MPR,MTK,MI)

% Function [xmelt,hlat]=Melt_fraction(ppa,mtk,rock)
%
% This function compute melt fraction (xmelt) and latent heat (hlat)
% at given pressure (ppa), temperature (mtk)and rock type (rock)
% Function returns solution for melt fraction (xmelt) and
% respective latent heat increment (hlat)

% Calculate melt fraction using marker rock type
% Standard 19
% MatProp{1} = 'Sticky Air';
% MatProp{2} = 'Sticky Water';
% MatProp{3} = 'Sediments';
% MatProp{4} = 'Upper Oceanic Crust';
% MatProp{5} = 'Lower Oceanic Crust';
% MatProp{6} = 'Upper Continental Crust';
% MatProp{7} = 'Lower Continental Crust';
% MatProp{8} = 'Oceanic Lithospheric Mantle';
% MatProp{9} = 'Continental Lithospheric Mantle';
% MatProp{10} = 'Upper Mantle';
% MatProp{11} = 'Hydrated Mantle';
% MatProp{12} = 'Plume Material';
% MatProp{13} = 'Eclogite';
% MatProp{14} = 'Lower Mantle'; % Perovskite
% MatProp{15} = 'PvLith';       % Equivalent to ower mantle used to track lithosphere.
% MatProp{16} = 'Upper Oceanic Plateau';
% MatProp{17} = 'Lower Oceanic Plateau';
% MatProp{18} = 'Upper Continental Crust Stripe';
% MatProp{19} = 'Lower Continental Crust Stripe';

marknum = length(MI);
XMELT = zeros(size(MI));
HLAT = XMELT;

for i = 1:marknum
    
    tl=0; % Liquidus temperature
    mtk = MTK(i);
    rock_type = MI(i);
    P = MPR(i)*1e-6; % MPa
    
    switch rock_type
                    
        % 3 = Sediments: latent heat 300 kJ/kg
        case 3
            % Solidus Temperature
            if (P<1200)
                ts=889+17900/(P+54)+20200/(P+54)^2;
            else
                ts=831+0.06*P;
            end
            % Liquidus temperature
            tl=1262+0.09*P;
            % Latent heat
            HL=300000;
            
        % 4, 16 = Basalt: latent heat 380 kJ/kg
        case {4,16}
            % Solidus Temperature
            if (P<1600)
                ts=973-70400/(P+354)+77800000/(P+354)^2;
            else
                ts=935+0.0035*P+0.0000062*P^2;
            end
            % Liquidus temperature
            tl=1423+0.105*P;
            % Latent heat
            HL=380000;
            
        % 5,17 = Gabbro: latent heat 380 kJ/kg
        case {5,17}
            % Solidus Temperature
            if (P<1600)
                ts=973-70400/(P+354)+77800000/(P+354)^2;
            else
                ts=935+0.0035*P+0.0000062*P^2;
            end
            % Liquidus temperature
            tl=1423+0.105*P;
            % Latent heat
            HL=380000;
            
        % 6,18 = Granodiorite, Upper continental crust: latent heat 300 kJ/kg
        case {6,18}
            % Solidus Temperature
            if (P<1200)
                ts=889+17900/(P+54)+20200/(P+54)^2;
            else
                ts=831+0.06*P;
            end
            % Liquidus temperature
            tl=1262+0.09*P;
            % Latent heat
            HL=300000;
            
        % 7,19 = Diorite, Lower continental crust: latent heat 380 kJ/kg
        case {7,19}
            % Solidus Temperature
            if (P<1600)
                ts=973-70400/(P+354)+77800000/(P+354)^2;
            else
                ts=935+0.0035*P+0.0000062*P^2;
            end
            % Liquidus temperature
            tl=1423+0.105*P;
            % Latent heat
            HL=380000;
                
        % 8,9 = Lithospheric mantle (dry): latent heat 400 kJ/kg
        case {8,9}
            % Solidus Temperature
            if (P<10000)
                ts=1394+0.132899*P-0.000005104*P^2;
            else
                ts=2212+0.030819*(P-10000);
            end
            % Liquidus temperature
            tl=2073+0.114*P;
            % Latent heat
            HL=400000;
            
        % 10 = Asthenospheric mantle (dry): latent heat 400 kJ/kg
        case 10
            % Solidus Temperature
            if (P<10000)
                ts=1394+0.132899*P-0.000005104*P^2;
            else
                ts=2212+0.030819*(P-10000);
            end
            % Liquidus temperature
            tl=2073+0.114*P;
            % Latent heat
            HL=400000;
            
        % 11 = Hydrated mantle (wet): latent heat 400 kJ/kg
        case 11
            % Solidus Temperature
            if (P<2400)
                ts=1240+49800/(P+323);
            else
                ts=1266-0.0118*P+0.0000035*P^2;
            end
            % Liquidus temperature
            tl=2073+0.114*P;
            % Latent heat
            HL=400000;
            
        % 12 = Plume Material, Asthenospheric mantle (dry): latent heat 400 kJ/kg -- plume
        case 12
            % Solidus Temperature
            if (P<10000)
                ts=1394+0.132899*P-0.000005104*P^2;
            else
                ts=2212+0.030819*(P-10000);
            end
            % Liquidus temperature
            tl=2073+0.114*P;
            % Latent heat
            HL=400000;
    end
    
    % Melt fraction and latent heat calc, check
    if (tl>0)
        % Solidus and liquidus must not entersect
        % in the extrapolation region
        if (ts>tl-100)
            ts=tl-100;
        end
        % Melt fraction
        XMELT(i) = (mtk-ts)/(tl-ts);
        if (XMELT(i) < 0)
            XMELT(i) = 0;
        end
        if (XMELT(i) > 1)
            XMELT(i) = 1;
        end
        % Latent heat calc
        HLAT(i)=HL*XMELT(i);
    end
    
end

