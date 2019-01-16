function plot_idos(E,M);
     %
     % Integrate DoS in M against energy scale in E
     % And plot
     %
     E_h = flipud(E);
     dos_h = flipud(M);

     idos(1) = 0;

     for idx = 2:length(E_h)
        idos(idx) = trapz(E_h(1:idx),dos_h(1:idx));
     end

     clf
     hold on

     uc_area = 5.7; % unit cell area in A^-2, roughly (need exact #!)

     idos_scaled = idos*11; % int units [per unit-cell]
     dos_scaled = dos_h*11/uc_area;    % in units of eV^-1 A^-1

     clf

     subplot(2,1,1)
     plot(E_h,idos_scaled);
     xlabel('Energy [eV]')
     ylabel('IDoS [per unit-cell]')
     axis([E_h(1) E_h(end) 0 idos_scaled(end)])

     set(groot,'DefaultTextInterpreter','Latex')

     subplot(2,1,2)
     plot(E_h,dos_scaled);
     xlabel('Energy [eV]');
     ylabel('DoS [eV$^{-1}$ \AA $^{-2}$]')

end
