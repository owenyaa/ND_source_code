function [reg_c,rho_c,eta_c] = l_corner(rho,eta,reg_param,U,s,b,method,M)
%L_CORNER Locate the "corner" of the L-curve.
%
% [reg_c,rho_c,eta_c] =
%        l_corner(rho,eta,reg_param)
%        l_corner(rho,eta,reg_param,U,s,b,method,M)
%        l_corner(rho,eta,reg_param,U,sm,b,method,M) ,  sm = [sigma,mu]
%
% Locates the "corner" of the L-curve in log-log scale.
%
% It is assumed that corresponding values of || A x - b ||, || L x ||,
% and the regularization parameter are stored in the arrays rho, eta,
% and reg_param, respectively (such as the output from routine l_curve).
%
% If nargin = 3, then no particular method is assumed, and if
% nargin = 2 then it is issumed that reg_param = 1:length(rho).
%
% If nargin >= 6, then the following methods are allowed:
%    method = 'Tikh'  : Tikhonov regularization
%    method = 'tsvd'  : truncated SVD or GSVD
%    method = 'dsvd'  : damped SVD or GSVD
%    method = 'mtsvd' : modified TSVD,
% and if no method is specified, 'Tikh' is default.  If the Spline Toolbox
% is not available, then only 'Tikh' and 'dsvd' can be used.
%
% An eighth argument M specifies an upper bound for eta, below which
% the corner should be found.

% The following functions from the Spline Toolbox are needed if 当所用方法不是tikh
% 或者 dsvd 的时候样条函数工具盒才会用到
% method differs from 'Tikh' or 'dsvd':
% fnder, ppbrk, ppmak, ppual, sp2pp, sorted, spbrk, spmak, sprpp.

% Per Christian Hansen, IMM, 04/16/98.

% Set default regularization method.
if (nargin <= 3)
  method = 'none';
  if (nargin==2), reg_param = [1:length(rho)]'; end
else
  if (nargin==6), method = 'Tikh'; end
end

% Set threshold for skipping very small singular values in the
% analysis of a discrete L-curve.
s_thr = eps;  % Neglect singular values less than s_thr.

% Set default parameters for treatment of discrete L-curve.
deg   = 2;  % Degree of local smooting polynomial.
q     = 2;  % Half-width of local smoothing interval.
order = 4;  % Order of fitting 2-D spline curve.                            %%%%%

% Initialization.
if (length(rho) < order)
  error('Too few data points for L-curve analysis')
end
if (nargin > 3)
  [p,ps] = size(s); [m,n] = size(U);
  if (ps==2), s = s(p:-1:1,1)./s(p:-1:1,2); U = U(:,p:-1:1); end
  beta = U'*b; xi = beta./s;
end

% Restrict the analysis of the L-curve according to M (if specified).
if (nargin==8)
  index = find(eta < M);
  rho = rho(index); eta = eta(index); reg_param = reg_param(index);
end

if (strncmp(method,'Tikh',4) | strncmp(method,'tikh',4))                   %%%%%默认的正则化方法

  % The L-curve is differentiable; computation of curvature in
  % log-log scale is easy.
  
  % Compute g = - curvature of L-curve.
  g = lcfun(reg_param,s,beta,xi);
  
  % Locate the corner.  If the curvature is negative everywhere,
  % then define the leftmost point of the L-curve as the corner.
  [gmin,gi] = min(g);
  reg_c = fminbnd('lcfun',...
    reg_param(min(gi+1,length(g))),reg_param(max(gi-1,1)),...
    [],s,beta,xi); % Minimizer.                                            %%%%% fminbnd去最小，后面的[],s,beta,xi看不懂
  kappa_max = - lcfun(reg_c,s,beta,xi); % Maximum curvature.

  if (kappa_max < 0)
    lr = length(rho);
    reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
  else
    f = (s.^2)./(s.^2 + reg_c^2);
    eta_c = norm(f.*xi);
    rho_c = norm((1-f).*beta);
    if (m>n), rho_c = sqrt(rho_c^2 + norm(b - U*beta)^2); end               %%%%这是几个意思
  end
                                                                           %%%% 默认方法就到这
elseif (strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4) | ...
        strncmp(method,'mtsv',4) | strncmp(method,'none',4))

  % The L-curve is discrete and may include unwanted fine-grained
  % corners.  Use local smoothing, followed by fitting a 2-D spline
  % curve to the smoothed discrete L-curve.

  % Check if the Spline Toolbox exists, otherwise return.
  if (exist('spdemos')~=2)
    error('The Spline Toolbox in not available so l_corner cannot be used')
  end

  % For TSVD, TGSVD, and MTSVD, restrict the analysis of the L-curve
  % according to s_thr.
  if (nargin > 3)
    index = find(s > s_thr);
    rho = rho(index); eta = eta(index); reg_param = reg_param(index);
    s = s(index); beta = beta(index); xi = xi(index);
  end

  % Convert to logarithms.
  lr = length(rho);
  lrho = log(rho); leta = log(eta); slrho = lrho; sleta = leta;

  % For all interior points k = q+1:length(rho)-q-1 on the discrete
  % L-curve, perform local smoothing with a polynomial of degree deg
  % to the points k-q:k+q.
  v = [-q:q]'; A = zeros(2*q+1,deg+1); A(:,1) = ones(length(v),1);
  for j = 2:deg+1, A(:,j) = A(:,j-1).*v; end
  for k = q+1:lr-q-1
    cr = A\lrho(k+v); slrho(k) = cr(1);
    ce = A\leta(k+v); sleta(k) = ce(1);
  end                                                                      %%%%% 优化？

  % Fit a 2-D spline curve to the smoothed discrete L-curve.
  sp = spmak([1:lr+order],[slrho';sleta']);
  pp = ppbrk(sp2pp(sp),[4,lr+1]);

  % Extract abscissa and ordinate splines and differentiate them.
  % Compute as many function values as default in spleval.
  P     = spleval(pp);  dpp   = fnder(pp);                                 % 对样条函数求导
  D     = spleval(dpp); ddpp  = fnder(pp,2);
  DD    = spleval(ddpp);
  ppx   = P(1,:);       ppy   = P(2,:);
  dppx  = D(1,:);       dppy  = D(2,:);
  ddppx = DD(1,:);      ddppy = DD(2,:);

  % Compute the corner of the discretized .spline curve via max. curvature.
  % No need to refine this corner, since the final regularization
  % parameter is discrete anyway.
  % Define curvature = 0 where both dppx and dppy are zero.
  k1    = dppx.*ddppy - ddppx.*dppy;
  k2    = (dppx.^2 + dppy.^2).^(1.5);
  I_nz  = find(k2 ~= 0);
  kappa = zeros(1,length(dppx));
  kappa(I_nz) = -k1(I_nz)./k2(I_nz);
  [kmax,ikmax] = max(kappa);
  x_corner = ppx(ikmax); y_corner = ppy(ikmax);

  % Locate the point on the discrete L-curve which is closest to the
  % corner of the spline curve.  Prefer a point below and to the
  % left of the corner.  If the curvature is negative everywhere,
  % then define the leftmost point of the L-curve as the corner.
  if (kmax < 0)
    reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
  else
    index = find(lrho < x_corner & leta < y_corner);
    if (length(index) > 0)
      [dummy,rpi] = min((lrho(index)-x_corner).^2 + (leta(index)-y_corner).^2);
      rpi = index(rpi);
    else
      [dummy,rpi] = min((lrho-x_corner).^2 + (leta-y_corner).^2);
    end
    reg_c = reg_param(rpi); rho_c = rho(rpi); eta_c = eta(rpi);
  end

elseif (strncmp(method,'dsvd',4) | strncmp(method,'dgsv',4))               

  % The L-curve is differentiable; computation of curvature in
  % log-log scale is easy.

  % Compute g = - curvature of L-curve.
  g = lcfun(reg_param,s,beta,xi,1);

  % Locate the corner.  If the curvature is negative everywhere,
  % then define the leftmost point of the L-curve as the corner.
  [gmin,gi] = min(g);
  reg_c = fminbnd('lcfun',...
    reg_param(min(gi+1,length(g))),reg_param(max(gi-1,1)),...
    [],s,beta,xi,1); % Minimizer.
  kappa_max = - lcfun(reg_c,s,beta,xi,1); % Maximum curvature.

  if (kappa_max < 0)
    lr = length(rho);
    reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
  else
    f = s./(s + reg_c);
    eta_c = norm(f.*xi);
    rho_c = norm((1-f).*beta);
    if (m>n), rho_c = sqrt(rho_c^2 + norm(b - U*beta)^2); end
  end

else, error('Illegal method'), end