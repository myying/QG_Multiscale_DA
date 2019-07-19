function [wi_max,wr_max,psiv] = qggrz(z,rho,U,V,F,betax,betay,kvec,lvec,dim)

%  [wi_max,wr_max,psiv] = qggrz(z,rho,U,V,F,betax,betay,kvec,lvec,dim)
%     Calculates the growth rate for the stratified QG model 
%     with a density profile, rho(:), mean velocity components
%     U(:) and V(:), each with values specified at vertical
%     postions z(:).  
%
%     Inputs:
%
%       z     = z coordinates for values in vector rho
%       rho   = background potential density profile
%       U     = mean zonal velocity
%       V     = mean meridional velocity
%       F     = f^2*rho0/g * [L0^2/(H0*drho0)]  
%       betax = beta0      * [L0^2/U0]
%       betay = beta0      * [L0^2/U0]
%
%       kvec  = vector of zonal wavenumbers
%       lvec  = vector of meridional wavenumbers
%
%     Optional inputs: 
%
%       dim  = 0: rho and z normalized as in spectral QG model (default)
%              1: rho and z not normalized (values in braces []
%                 above set to 1)
%             
%     where
%
%       rho0  = average of rho,
%       drho0 = (rho(nz) - rho(1))/(nz-1), 
%       nz    = length(rho), 
%       H0    = total depth
%       f0    = Coriolis parameter
%       g     = 9.81 m/s^2.  
%
%     Outputs:   
%
%        wi_max(length(kvec),length(lvec)):  
%               array of maximum growth rates (wrt vertical
%               wavenumber) for each wavenumber pair (k,l)
%        wr_max(length(kvec),length(lvec)):   
%               array of real part of frequency corresponding to
%               fastest growing mode for each (k,l).  If (k,l) is
%               stable, array contains largest frequency.
%        psiv(length(kvec),length(lvec),nz):
%               array of eigenfunctions that correspond to wi_max
%               for each (k,l)
%
%     See PV_STRETCH_OPZ

switch nargin
  case 10
  case 9, dim=0;
  otherwise, error('need 9 or 10 arguments')
end

% Grid sizes
nkx = length(kvec);  nky = length(lvec); nz = max(length(U),length(V));

% Columnize
z = z(:); rho = rho(:); U = U(:); V = V(:); 

% Find layer thicknesses dz
dz = get_dz(z);

% Stretching operator
G = pv_stretch_opz(z,rho,F,dim);

% Mean PV gradients
Q_x    = betax + G*V;
Q_y    = betay - G*U;

% Output arrays
wi_max = zeros(nkx,nky);  
wr_max = zeros(nkx,nky);

vec = 0;
if nargout>2
  psiv = zeros(nkx,nky,nz);
  vec  = 1;
end
  
kc = 1;
for k = kvec
  lc = 1;
  for l = lvec
    
    % Calculate K-dependent terms
    K2 = k^2 + l^2;          % Total wavenumber
    KdotU = k*U + l*V;       % Linear advection..
    KxdelQ = k*Q_y - l*Q_x;
    
    % Matrices:  w B psi = A psi
    B = G - K2*eye(nz);
    
    % To remove vorticity from top and bottom levels, uncomment the
    % following lines
    %B(1,:) = G(1,:);
    %B(end,:) = G(end,:);
    
    A = diag(KxdelQ) + diag(KdotU) * B;

    % Solve eigenvalue problem
    if vec==1
      [evec,D] = eig(A,B);
      w = diag(D);
    else
      w = eig(A,B);
    end
    if max(abs(imag(w))) > eps
      [wi_max(kc,lc),ind] = max(imag(w));
      wr_max(kc,lc) = real(w(ind));
    else
      [wr_max(kc,lc),ind] = max(real(w));
      wi_max(kc,lc) = imag(w(ind));
    end
    if vec==1, 
      psiv(kc,lc,:) = evec(:,ind); 
    end

    % Perhaps you'd rather have the amplitude and phase as outputs...
    % amp(kc,lc,:)=sqrt(real(evec(:,ind)).^2+imag(evec(:,ind)).^2);
    % phase(kc,lc,:)=atan2(imag(evec(:,ind)),real(evec(:,ind)));
    
    lc = lc+1;
  end
  
  kc = kc+1;
end

