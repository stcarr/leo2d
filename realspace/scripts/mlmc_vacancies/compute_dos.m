%
% Compute DoS approximation from eigenvalues
%
function [dos,eng] = compute_dos(eigvals,de,regeps);
emin = min(eigvals);
emax = max(eigvals);
eng = [emin:de:emax];

ne = length(eng);
for i=1:ne
    
    dos(i) = -1/(ne*pi)*imag(sum(1./(eng(i) - eigvals + 1i*regeps)));
    
end

plot(eng,dos); hold on;
xlabel('E');
ylabel('DoS');
title('Density of states');
