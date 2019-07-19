function G = pv_stretch_opz(z,rho,F,dim,Fe)

%  G = pv_stretch_opz(z,rho,F,dim,Fe)
%
%     Create tridiagonal matrix G that opearates on streamfunction psi to
%     get vortex stretching term in QG potential vorticity.  i.e.
%
%       q_s = G psi = d_z S d_z psi, 
%
%     where 
% 
%       S = -(rho0 f^2/g) (d_z rho)^{-1} = f^2/N^2
%  
%     Inputs:
%
%       z    = z coordinates for values in vector rho
%       rho  = background potential density profile 
%       F    = f^2*rho0/g [L^2/(H*drho)]  
%              (set []=1 above for dimensional, dim=1 below)
%
%     Optional inputs: 
%
%       dim  = 0: rho and z normalized as in SQG model (default)
%              1: rho and z not normalized 
%       Fe   = F*(rho(2)-rho(1))/(rho(1)-rho_top) (default=0)
%              (for external deformation scale at upper surface)
%             
%     where
%
%       rho0 = average of rho,
%       drho = (rho(nz) - rho(1))/(nz-1), 
%       nz = length(rho), 
%       H = total depth
%       f = Coriolis parameter
%       g = 9.81 m/s^2.  
%

switch nargin
  case 5
  case 4, Fe=0;
  case 3, Fe=0; dim=0;
  otherwise, error('need 3, 4, or 5 arguments')
end

nz = length(z);

rho = rho(:);  
z = z(:); 

dz = get_dz(z);
drho = rho(2:end) - rho(1:end-1); % length nz-1

if dim==0
  drho0 = (rho(end)-rho(1))/(nz-1);
  drho = drho/drho0;  
  dz = dz/sum(dz);
end

% Make subdiagonal, diagonal and superdiagonal for tridiagonal
% stretching operator

sub = zeros(nz-1,1);
mid = zeros(nz,1);
sup = zeros(nz-1,1);

sub   = F./(dz(2:end).*drho);
sup   = F./(dz(1:end-1).*drho);

mid(1)      = -F/(dz(1)*drho(1)) - Fe/dz(1);
sup(1)      =  F/(dz(1)*drho(1));
mid(end)    = -F/(dz(end)*drho(end));
sub(end)    =  F/(dz(end)*drho(end));

if nz>2 
   mid(2:end-1) = -F./dz(2:end-1) .* (1./drho(1:end-1) + 1./drho(2:end));
end

G = diag(sub,-1) + diag(mid,0) + diag(sup,1);

