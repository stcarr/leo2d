function plot_CCCM(E,O,figtitle,ztitle,type)
%
% Plot current-current correlation measure type observable
% E - energy mesh
% O - Observable
% figtitle - string for the figure title
% ztitle - value plotted on z axis
% type   - figure type, type = 1 ... surf()

 if type == 1
    surf(E,E,O,'EdgeColor','none'); hold on; % simple plot
    view(2); axis square;
    xlabel('Energy E_1'); ylabel('Energy E_2');zlabel(ztitle);
    title(figtitle);
 end
 
end

