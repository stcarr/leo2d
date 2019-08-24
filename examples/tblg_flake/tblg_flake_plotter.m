clear all;
dos = dlmread('tblg_flake.cheb');
pos = dlmread('tblg_flake_pos.dat');
%kpts = dlmread('tblg_flake.kpts');

ax_m = 1.5;

clf
subplot(3,1,1)
hold on
box on
title('twisted bilayer graphene (flake)')
plot(dos(1,:),dos(2,:)+dos(3,:))
% 2,3 correspond to A/B orbitals at center unit-cell 

ylabel('DoS')
xlabel('Energy (eV)')
%plot([bands; bands(1,:)],'k')

subplot(3,1,[2 3])
hold on
box on
title('orbital positions')
l1_pos = pos( [pos(:,4) == 0],:);
plot3(l1_pos(:,1),l1_pos(:,2),l1_pos(:,3),'.')
l2_pos = pos( [pos(:,4) == 1],:);
plot3(l2_pos(:,1),l2_pos(:,2),l2_pos(:,3),'.')
axis equal

%{
subplot(2,2,4)
hold on
box on
title('k point sampling')
plot(kpts(:,1),kpts(:,2),'k')
axis equal
%}
