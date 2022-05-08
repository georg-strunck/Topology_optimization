% Function: Infill optimization
% Author: Jun Wu (j.wu-1@tudelft.nl)
% Version: 2017-06-19

% Example Infill(200,100,[2],600) or Infill(200,100,[1,2],600),   Infill(100,50,[1,2],300)

function Infill(nelx,nely,mdof,nloop)
%mdof[1,2]
% 1 total volume
% 2 upper bound
choice=[1,2];
for y = choice
    if y==1
        mdof=[1]
    elseif y==2
        mdof=[1,2]
    end

close all;
mkdir('images');
volfrac = 0.2;  % total volume constraint --> for normal topopt and in general how much of designspace is filled, constant (no max or min, exact)
vol_max = 0.3;  % local volume constraint --> how much of circle with radius r_hat is filled MAX
penal = 1;      % stiffness penalty       --> how important is stiffness? --> 3=solid, 1=lattice/varying density
p = 16;         % pNorm                   --> used for constraints differentiation approximation
r_hat = 15;      % pNorm radius            --> radius of local volume
rmin = 1.6;     % density filter radius
move = 0.01;    % limited move for the design variables
beta = 1;       % beta continuation
eta = 0.5;      % projection threshold, fixed at 0.5
vol_max_pNorm = (nelx*nely*vol_max^p)^(1/p);

%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
iLoad = 5;
Fsparse = sparse(2*(nely+1)*(nelx+1),1);
if iLoad == 0
    Fsparse(2*nelx*(nely+1)+2*(nely/2)+2,1) = -1;
    fixeddofs   = union([1:1:2],[2*(nely):1:2*(nely+1)]); %#ok<*NBRAK>
elseif iLoad == 1
    Fsparse(2,1) = -1;
    fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
elseif iLoad == 2
    Fsparse(2*nelx*(nely+1)+2*(nely/2)+2,1) = -1;
    fixeddofs   = union([2*(nely):1:2*(nely+1)],[2*(nelx+1)*(nely+1)-2:1:2*(nelx+1)*(nely+1)]);
elseif iLoad == 3
    Fsparse(2*(nely+1)*(nelx+1),1) = -1
    Fsparse(2*(nely+1)*(nelx)+nely+2,1) = -1;
    fixeddofs = union([1:1:2*(nely+1)],[1]);
elseif iLoad == 4                                                 % standing up loadcase
    %Fsparse(2*(nelx*7/10+1)*(nely+1)+2*(nely/4+1)+1) = -0.6;                   %seat in air
    Fsparse(2*(nelx*3/10+1)*(nely+1)+2*(nely/4+1)+1) = -0.2;                   %steering bar at top
    Fsparse(2*(nelx*6/10+1)*(nely+1)) = -0.8;                                  %pedal on ground
    fixeddofs = union([2*(nelx*2/10+1)*(nely+1)], [2*(nelx*9/10+1)*(nely+1)]); %wheel axes on ground
elseif iLoad == 5
    Fsparse(2*(nelx*7/10+1)*(nely+1)+2*(nely/4+1)+1) = -0.6;                   %seat in air
    Fsparse(2*(nelx*3/10+1)*(nely+1)+2*(nely/4+1)+1) = -0.1;                   %steering bar at top
    Fsparse(2*(nelx*6/10+1)*(nely+1)) = -0.3;                                  %pedal on ground
    fixeddofs = union([2*(nelx*2/10+1)*(nely+1)+2*(nely*9/10+1)], [2*(nelx*9/10+1)*(nely+1)+2*(nely*9/10+1)]); %wheel axes in air
end
U = zeros(2*(nely+1)*(nelx+1),1);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

%% PREPARE PDE FILTER
edofVecF = reshape(nodenrs(1:end-1,1:end-1),nelx*nely,1);
edofMatF = repmat(edofVecF,1,4)+repmat([0 nely+[1:2] 1],nelx*nely,1);
iKF = reshape(kron(edofMatF,ones(4,1))',16*nelx*nely,1);
jKF = reshape(kron(edofMatF,ones(1,4))',16*nelx*nely,1);
iTF = reshape(edofMatF,4*nelx*nely,1);
jTF = reshape(repmat([1:nelx*nely],4,1)',4*nelx*nely,1);
sTF = repmat(1/4,4*nelx*nely,1);
TF = sparse(iTF,jTF,sTF);

Rmin = r_hat/2/sqrt(3);
KEF = Rmin^2*[4 -1 -2 -1; -1  4 -1 -2; -2 -1  4 -1; -1 -2 -1  4]/6 + ...
             [4  2  1  2;  2  4  2  1;  1  2  4  2;  2  1  2  4]/36;
sKF = reshape(KEF(:)*ones(1,nelx*nely),16*nelx*nely,1);
KF = sparse(iKF,jKF,sKF);
LF = chol(KF,'lower');

%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xTilde = x;
xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
xold1 = reshape(x,[nely*nelx,1]);
xold2 = reshape(x,[nely*nelx,1]);
low = 0;
upp = 0;

loopbeta = 0;
loop = 0;
change = 1;
%% START ITERATION

% store results
c_hist = zeros(nloop,1);        % compliance
vol_hist = zeros(nloop,1);      % volume
change_hist = zeros(nloop,1);   % maximum design change
sharp_hist = zeros(nloop,1);    % sharpness
cons_hist = zeros(nloop,2);     % constraints


while change > 0.0001 && loop < nloop
    loopbeta = loopbeta+1;
    loop = loop+1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    % K(spring constant)*U(deformation)=F(force) --> K\F solves set of lin equ for U
    U(freedofs) = K(freedofs,freedofs)\Fsparse(freedofs);               

    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    ce_disp=sum((Emin+xPhys.^penal*(E0-Emin)).*ce);          % compliance per column
    c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;

    dv = ones(nely,nelx);

    x_pde_hat = (TF'*(LF'\(LF\(TF*xPhys(:)))));
    dfdx_pde = (sum(x_pde_hat.^p))^(1/p-1) * x_pde_hat.^(p-1);
  
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    dx = beta * (1-tanh(beta*(xTilde-eta)).*tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    dc(:) = H*(dc(:).*dx(:)./Hs);
    dv(:) = H*(dv(:).*dx(:)./Hs);

    %% UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    % METHOD OF MOVING ASYMPTOTES (MMA)
    m = size(mdof,2);
    n = nelx*nely;  
    
    df0dx = reshape(dc,[nelx*nely,1]);
    df0dx2 = 0*df0dx;
    dfdx = zeros(3,nelx*nely);
    dfdx(1,1:nelx*nely) = reshape(dv,[1,nelx*nely])/(nelx*nely*volfrac);
    dfdx(2,1:nelx*nely) = TF'*(LF'\(LF\(TF*dfdx_pde(:))));
    
    ic = 2;
    tmp = reshape(dfdx(ic,:),[nely,nelx]);
    dfdx(ic,:) = reshape(H*(tmp(:).*dx(:)./Hs),[1,nelx*nely]);   
    
    dfdx2 = 0*dfdx;

    iter = loopbeta;
    xval = reshape(x,[nelx*nely,1]);
    xmin=max(0.0,xval-move);
    xmax=min(1,xval+move);

    f0val = c;
    fval = zeros(2,1);
    fval(1,1) = sum(sum(xPhys)) / (nelx*nely*volfrac) - 1;
    fval(2,1) = (sum(x_pde_hat.^p))^(1/p)- vol_max_pNorm;
    
    a0 = 1;
    a = zeros(m,1);     
    c_ = ones(m,1)*1000;
    d = zeros(m,1);
    [xmma,ymma,zmma,lam,xsi,eta_,mu,zet,s,low,upp] = ...
        mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,...
        f0val,df0dx,df0dx2,fval(mdof),dfdx(mdof,:),dfdx2(mdof,:),low,upp,a0,a,c_,d);
    xnew = reshape(xmma,[nely,nelx]);
    xold2 = xold1;
    xold1 = xval;
    
    xTilde(:) = (H*xnew(:))./Hs;
    xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    
    xPhys(1:nely*2/10-2, 1:nelx*7/10+2) = 0;                            % vision and bottom + torso free of structure
    i = 0;                                                              % allow space for the front wheel to turn
    for m = nelx*3/10:nelx*5/10
        i = i+1;
        j = 0;
        for n = nely*6/10:nely
            j = j+1;
            if (j<=i)
                xPhys(m,n) = 0;
            end
        end
    end
%       [ii,jj] = find(Fsparse);
%       ii
%       column = floor(fix((ii)/((nelx+1))))
%       row = round(rem((ii),((nely+1)))/2)
%       size(xPhys)
%       for s=1:length(ii)
%           cross = 3;   % has to be uneven! (for middle of cross)
%           crosstep = fix(cross/2)
%           if row(s)==0
%               row(s)=nely-(crosstep-1);
%           end
%           if column(s)==0
%               column(s)=nelx-(crosstep-1);
%           end
%           if row(s)<cross
%               row(s)=cross-crosstep+1;
%           end
%           if column(s)<cross
%               column(s)=cross-crosstep+1;
%           end
%           row(s), column(s)
%           xPhys(row(s)-1-crosstep:row(s)-1+crosstep, column(s)-1)=1;
%           xPhys(row(s)-1, column(s)-1-crosstep:column(s)-1+crosstep)=1;
%           size(xPhys)
%       end


    change = max(abs(xnew(:)-x(:)));
    x = xnew;

    %% PRINT RESULTS
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
        ' Vol.: ' sprintf('%6.3f',sum(sum(xPhys))/(nelx*nely)) ...
        ' Ch.: ' sprintf('%6.3f',change) ...
        ' Cons.: ' sprintf('%6.3f',fval)]);

    %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if beta < 100 && (loopbeta >= 40 || change <= 0.001)
        beta = 2*beta;
        loopbeta = 0;
        change = 1;
        penal = min(3, penal + 0.5);
        fprintf('Parameter beta increased to %g.\n',beta);
        fprintf('Parameter penal increased to %g.\n',penal);
    end
    
    %% Store current values
    c_hist(loop,1) = c;
    vol_hist(loop,1) = sum(sum(xPhys))/(nelx*nely);
    change_hist(loop,1) = change;
    cons_hist(loop,:) = fval;
    sharp_hist(loop,1) = 4*sum(sum(xPhys.*(ones(nely, nelx)-xPhys))) / (nely*nelx);

    %% PLOT DENSITIES
    figure(1);
    set(1, 'Position', [100, 450, 540, min(100+540*nely/nelx,540)]);
    colormap(gray); imagesc(-xPhys, [-1 0]); axis equal; axis tight; axis off; drawnow;
    %title('\rho');

%     filename1 = sprintf('images\\rho-It%i.png',loop);
%     saveas(1,filename1,'png');
end

%%Damage calculation & plot

xPhys_dam = xPhys;
xPhys_dam(1:nely*4/10, nelx*5/10-5:nelx*5/10+5) = 0;
sK_dam = reshape(KE(:)*(Emin+xPhys_dam(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K_dam = sparse(iK,jK,sK_dam); K_dam = (K_dam+K_dam')/2;
% K(spring constant)*U(deformation)=F(force) --> K\F solves set of lin equ for U
U_dam(freedofs) = K_dam(freedofs,freedofs)\Fsparse(freedofs);  
ce_dam = reshape(sum((U_dam(edofMat)*KE).*U_dam(edofMat),2),nely,nelx);
ce_disp_dam=sum((Emin+xPhys_dam.^penal*(E0-Emin)).*ce_dam);
c_dam = sum(sum((Emin+xPhys_dam.^penal*(E0-Emin)).*ce_dam));

if y==1
        ce_disp_reg = ce_disp;
        ce_disp_dam_reg = ce_disp_dam;
        xPhys_reg = xPhys;
        xPhys_dam_reg = xPhys_dam;
elseif y==2
        ce_disp_por = ce_disp;
        ce_disp_dam_por = ce_disp_dam;
        xPhys_por = xPhys;
        xPhys_dam_por = xPhys_dam;
end

end

%% Grid plotting

xPhys_grid = xPhys;
xPhys_grid(:) = 0;
for u = round(linspace(0,nelx-1,(nelx*volfrac/2)))
    xPhys_grid(:,u+1)=1;
end
for q = round(linspace(0,nely-1,(nely*volfrac/2)))
    xPhys_grid(q+1,:)=1;
end
% commented out code below adds extra material (cross-shape) close to where the forces
% are applied, not needed for paper
% [ii,jj] = find(Fsparse);
% ii
% column = floor(fix((ii)/((nelx+1))))
% row = round(rem((ii),((nely+1)))/2)
% for s=1:length(ii)
%   cross = 3;   % has to be uneven! (for middle of cross)
%   crosstep = fix(cross/2);
%   if row(s)==0
%       row(s)=nely-(crosstep-1);
%   end
%   if column(s)==0
%       column(s)=nelx-(crosstep-1);
%   end
%   if row(s)<cross
%       row(s)=cross-crosstep+1;
%   end
%   if column(s)<cross
%       column(s)=cross-crosstep+1;
%   end
%   xPhys_grid(row(s)-1-crosstep:row(s)-1+crosstep, column(s)-1)=1;
%   xPhys_grid(row(s)-1, column(s)-1-crosstep:column(s)-1+crosstep)=1;
% end
xPhys_grid(nely*8/10+5:nely-5, nelx*1/10:nelx*3/10-5)=1;
xPhys_grid(nely*8/10+5:nely-5, nelx*8/10-5:nelx-5)=1;
sK_grid = reshape(KE(:)*(Emin+xPhys_grid(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K_grid = sparse(iK,jK,sK_grid); K_grid = (K_grid+K_grid')/2;
% K(spring constant)*U(deformation)=F(force) --> K\F solves set of lin equ for U
U_grid(freedofs) = K_grid(freedofs,freedofs)\Fsparse(freedofs);  
ce_grid = reshape(sum((U_grid(edofMat)*KE).*U_grid(edofMat),2),nely,nelx);
ce_disp_grid=sum((Emin+xPhys_grid.^penal*(E0-Emin)).*ce_grid);%*10^(-6);
c_grid = sum(sum((Emin+xPhys_grid.^penal*(E0-Emin)).*ce_grid));

%%Grid with damage
xPhys_dam_grid = xPhys_grid;
xPhys_dam_grid(1:nely*4/10, nelx*5/10-5:nelx*5/10+5) = 0;
sK_dam_grid = reshape(KE(:)*(Emin+xPhys_dam_grid(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K_dam_grid = sparse(iK,jK,sK_dam_grid); K_dam_grid = (K_dam_grid+K_dam_grid')/2;
% K(spring constant)*U(deformation)=F(force) --> K\F solves set of lin equ for U
U_dam_grid(freedofs) = K_dam_grid(freedofs,freedofs)\Fsparse(freedofs);  
ce_dam_grid = reshape(sum((U_dam_grid(edofMat)*KE).*U_dam_grid(edofMat),2),nely,nelx);
ce_disp_dam_grid=sum((Emin+xPhys_dam_grid.^penal*(E0-Emin)).*ce_dam_grid);%*10^(-6);
c_dam_grid = sum(sum((Emin+xPhys_dam_grid.^penal*(E0-Emin)).*ce_dam_grid));

figure(1);
set(1, 'Position', [100, 450, 540, min(100+540*nely/nelx,540)]);
colormap(gray); imagesc(-xPhys_reg, [-1 0]); axis equal; axis tight; axis off; drawnow;
title('Regular TopOpt');
figure(2);
set(1, 'Position', [100, 450, 540, min(100+540*nely/nelx,540)]);
colormap(gray); imagesc(-xPhys_dam_reg, [-1 0]); axis equal; axis tight; axis off; drawnow;
title('Regular TopOpt Damage');
figure(3);
set(1, 'Position', [100, 450, 540, min(100+540*nely/nelx,540)]);
colormap(gray); imagesc(-xPhys_por, [-1 0]); axis equal; axis tight; axis off; drawnow;
title('Porous TopOpt');
figure(4);
set(1, 'Position', [100, 450, 540, min(100+540*nely/nelx,540)]);
colormap(gray); imagesc(-xPhys_dam_por, [-1 0]); axis equal; axis tight; axis off; drawnow;
title('Porous TopOpt Damage');
figure(5);
set(1, 'Position', [100, 450, 540, min(100+540*nely/nelx,540)]);
colormap(gray); imagesc(-xPhys_grid, [-1 0]); axis equal; axis tight; axis off; drawnow;
title('Manual Grid');
figure(6);
set(1, 'Position', [100, 450, 540, min(100+540*nely/nelx,540)]);
colormap(gray); imagesc(-xPhys_dam_grid, [-1 0]); axis equal; axis tight; axis off; drawnow;
title('Manual Grid Damage');

figure(7)
x = 1:nelx;
plot(x, ce_disp_dam_reg, x, ce_disp_dam_por, x, ce_disp_dam_grid)
legend('Regular TopOpt with damage', 'Porous TopOpt with damage', 'Manual Grid with damage')
figure(8)
x = 1:nelx;
plot(x, ce_disp_reg, x, ce_disp_por, x, ce_disp_grid)
legend('Regular TopOpt', 'Porous TopOpt', 'Manual Grid')

for t=1:8
    filename1 = sprintf('images\\notitle%i.png',t);
    saveas(t,filename1,'png');
end

% Publication
% Jun Wu, Niels Aage, Ruediger Westermann, Ole Sigmund, 
% Infill Optimization for Additive Manufacturing -- Approaching Bone-like Porous Structures
% IEEE Trans. on Visualization and Computer Graphics, 2017

% The code was developed based on the 110 line topology optiziation code, by E. Andreassen etc, 2011