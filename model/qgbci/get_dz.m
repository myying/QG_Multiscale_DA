function dz = get_dz(z);

%   dz = get_dz(z)
%
%   This version gives dz with length(z)

nz = length(z);
Dz = z(1:end-1)-z(2:end);
dz(1) = Dz(1);
dz(2:nz-1) = (Dz(1:end-1)+Dz(2:end))/2;
dz(nz) = Dz(end);

