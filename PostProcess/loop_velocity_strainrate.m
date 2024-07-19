function loop_velocity(inpath,nstart,dt,nend,xlims,ylims)
%% 

MPROP = 'MEII';
k = [nstart:dt:nend];
for step = k
    disp(step);
    plot_contour_extraction_velocity_strainrate_onlymat(inpath,step,xlims,ylims, MPROP);
end

end

% function loop_velocity(inpath,nstart,dt,nend)
% %% 
% 
% MPROP = 'MEII';
% k = [nstart:dt:nend];
% for step = k
%     disp(step);
%     %plot_contour_extraction_velocity_strainrate(inpath,step,xlims,ylims, MPROP);
%     plot_particles_plastic(inpath, step, 5e-16, 1e-13, 1,2)
% end
% 
% end