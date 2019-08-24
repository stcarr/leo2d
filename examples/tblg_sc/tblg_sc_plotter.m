clear all;
bands = dlmread('tblg_sc.cheb');
pos = dlmread('tblg_sc_pos.dat');
kpts = dlmread('tblg_sc.kpts');

ax_m = 1.5;

clf
subplot(2,2,[1 3])
hold on
box on
title('twisted bilayer graphene')
plot([bands; bands(1,:)],'k')
axis([1 size(bands,1)+1 -ax_m ax_m])

subplot(2,2,2)
hold on
box on
title('orbital positions')
l1_pos = pos( [pos(:,4) == 0],:);
plot3(l1_pos(:,1),l1_pos(:,2),l1_pos(:,3),'.')
l2_pos = pos( [pos(:,4) == 1],:);
plot3(l2_pos(:,1),l2_pos(:,2),l2_pos(:,3),'.')
axis equal

subplot(2,2,4)
hold on
box on
title('k point sampling')
plot(kpts(:,1),kpts(:,2),'k')
axis equal
