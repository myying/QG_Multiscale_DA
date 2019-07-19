function [kd,pm] = pmodesz(G,z);

% [kd,pm] = pmodesz(G,z)  
%
%     Calculates the deformation wavenmbers and vertical modes psi
%     solving
%
%     d/dz(f^2/N^2 dpsi/dz) = -kd^2 psi
%
%     where N^2 is tbe BV frequency and f is Coriolis.  Input is
%     matrix A such that, in discrete form, problem is
%
%     G psi = -kd^2 psi
%
%     where G can be obtained from pv_stretch_opz.  The vector z
%     contains the coordinates of the density rho needed to make G.
%     Vertical modes psi are given in matrix pm organized as pm(z,mode).
%     These eigenvectors psi_m(z) are normalized so that:
%     Integral(psi_m(z)*psi_n(z)))dz = Sum(dz*psi_m(z)*psi_n(z)) = delta_mn.
%
%     See also PV_STRETCH_OPZ

[V,D] = eig(G);

[kd,ri]=sort(sqrt(-diag(D)));
%kd(1)=0; 
pm(:,:)=V(:,ri);

% Now normalize pm so that <pm_n,pm_m> = delta_mn

nz = length(z);
dz = get_dz(z);
dz = dz(:)/sum(dz);
for m=1:nz
   alpha=(pm(:,m).*dz)'*pm(:,m);
   pm(:,m)=pm(:,m)/sqrt(alpha);
   if (pm(1,m)<0)
      pm(:,m) = -pm(:,m);
   end
end

