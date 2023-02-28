
% Build to run multithreaded
%mex FFLAGS='$FFLAGS -fopenmp' FOPTIMFLAGS='-O2' LDFLAGS='$LDFLAGS,-lgomp' mexmarker_rheology.F90 -compatibleArrayDims
mex FFLAGS='$FFLAGS -fopenmp' FOPTIMFLAGS='-O2' mexmarker_rheology.F90 -compatibleArrayDims
% no multithread
%mex FOPTIMFLAGS='-O2' mexmarker_rheology.F90 -compatibleArrayDims
 
%The rest the F90 files you can compile with
mex FOPTIMFLAGS='-O2' mexMETA2grid.F90 -compatibleArrayDims		
mex FOPTIMFLAGS='-O2' mexmovemarkersFast.F90 -compatibleArrayDims
mex FOPTIMFLAGS='-O2' mexMarkers2grid.F90 -compatibleArrayDims
mex FOPTIMFLAGS='-O2' mexmarker_pressure_strainrate.F90 -compatibleArrayDims		
mex FOPTIMFLAGS='-O2' mexsubgridS_diffusion.F90 -compatibleArrayDims
mex FOPTIMFLAGS='-O2' mexsubgridT_diffusion.F90 -compatibleArrayDims

% !cp mexmarker_rheology.mexmaci64 ..
% !cp mexmarker_pressure_strainrate.mexmaci64 ..
% !cp mexMETA2grid.mexmaci64 ..
% !cp mexmovemarkersFast.mexmaci64 ..
% !cp mexMarkers2grid.mexmaci64 ..
% !cp mexsubgridS_diffusion.mexmaci64 ..
% !cp mexsubgridT_diffusion.mexmaci64 ..