clear all;
bands = dlmread('mlg_nostrain.cheb');
pos = dlmread('mlg_nostrain_pos.dat');
kpts = dlmread('mlg_nostrain.kpts');
sc = dlmread('mlg_nostrain_supercell.dat');

strain_pos = dlmread('mlg_strain_pos.dat');
strain_bands = dlmread('mlg_strain.cheb');
strain_kpts = dlmread('mlg_strain.kpts');
strain_sc = dlmread('mlg_strain_supercell.dat');

ax_m = 11;

clf
subplot(2,2,[1 3])
hold on
box on
title('strained monolayer graphene')
plot([bands; bands(1,:)],'k')
plot([strain_bands; strain_bands(1,:)],'r')

axis([1 size(bands,1)+1 -ax_m+2 ax_m+2])

subplot(2,2,2)
hold on
box on
title('orbital positions')
for dx = -3:3
    for dy = -3:3
        dr1 = dx*sc(1,4:5) + dy*sc(2,4:5);
        dr2 = dx*strain_sc(1,4:5) + dy*strain_sc(2,4:5);

        plot3(dr1(1)+pos(:,1),dr1(2)+pos(:,2),pos(:,3),'.k')
        plot3(dr2(1)+strain_pos(:,1),dr2(2)+strain_pos(:,2),strain_pos(:,3),'.r')
        axis equal
    end
end

subplot(2,2,4)
hold on
box on
title('k point sampling')
plot(kpts(:,1),kpts(:,2),'k')
plot(strain_kpts(:,1),strain_kpts(:,2),'r')

axis equal
