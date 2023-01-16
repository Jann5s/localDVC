clear; close all


% Compiling
% -------------------------------

if true
    % debug = '-g';
    debug = '';

    cppfile = 'localDVC_kernel.cpp';
    mex(cppfile,debug,'-largeArrayDims','-DNOMINMAX')

end

% Reading some test images
% -------------------------------

% Define the path to the images
imagepath = 'virtimage3D';
imagefiles = dir(fullfile(imagepath, '*.raw'));
imagefiles = {imagefiles.name}';
images = strcat(imagepath, filesep, imagefiles);

refimage = images{1};
defimage = images{2};

siz = [160, 160, 160];

fid = fopen(refimage, 'r');
f = fread(fid, inf, '*uint8', 0);
fclose(fid);
fid = fopen(defimage, 'r');
g = fread(fid, inf, '*uint8', 0);
fclose(fid);

f = reshape(f, siz);
g = reshape(g, siz);

% convert to floating point
f = single(f);
g = single(g);

f = f ./ 255;
g = g ./ 255;

% Subset definition
% -------------------------------

C = round(siz./2) + [3,4,5];
L = 15 * [1, 1, 1];
threads = 0;

% number of dof (1: order0, 4: order1, 10: order2)
Ndof = 10;

% initial guess
a = zeros( 3 * Ndof,1);

% The equivalent code in Matlab (to verify)
% -------------------------------
tic

Ix = 1:3:3*Ndof;
Iy = 2:3:3*Ndof;
Iz = 3:3:3*Ndof;

x = (C(1) - ((L(1) - 1) / 2) + (0:L(1)-1));
y = (C(2) - ((L(2) - 1) / 2) + (0:L(2)-1));
z = (C(3) - ((L(3) - 1) / 2) + (0:L(3)-1));

x = x( x > 1 & x < siz(2) );
y = y( y > 1 & y < siz(1) );
z = z( z > 1 & z < siz(3) );

[X, Y, Z] = meshgrid(1:siz(2),1:siz(1),1:siz(3));

Xn = 2 * (X - C(1)) ./ L(1);
Yn = 2 * (Y - C(2)) ./ L(2);
Zn = 2 * (Z - C(3)) ./ L(3);

roi = false(siz);
roi(y,x,z) = true;

[fx, fy, fz] = gradient(f);

if Ndof == 1
    phi = ones( prod(siz), 1);
elseif Ndof == 4
    phi = ones( prod(siz), 4);
    phi(:,2) = Xn(:);
    phi(:,3) = Yn(:);
    phi(:,4) = Zn(:);
elseif Ndof == 10
    phi = ones( prod(siz), 10);
    phi(:, 2) = Xn(:);
    phi(:, 3) = Yn(:);
    phi(:, 4) = Zn(:);

    phi(:, 5) = Xn(:) .* Xn(:);
    phi(:, 6) = Yn(:) .* Yn(:);
    phi(:, 7) = Zn(:) .* Zn(:);

    phi(:, 8) = Yn(:) .* Zn(:);
    phi(:, 9) = Xn(:) .* Zn(:);
    phi(:,10) = Xn(:) .* Yn(:);
else
    error('Ndof %d',Ndof);
end

Ux = phi(roi,:) * a(Ix);
Uy = phi(roi,:) * a(Iy);
Uz = phi(roi,:) * a(Iz);

gt = interp3(X, Y, Z, g, X(roi) + Ux, Y(roi) + Uy, Z(roi) + Uz, 'cubic');

res = f(roi) - gt;
R = nan(siz);
R(roi) = res;

Lx = fx(roi) .* phi(roi,:);
Ly = fy(roi) .* phi(roi,:);
Lz = fz(roi) .* phi(roi,:);

b1 = zeros(3 * Ndof, 1);
b1(Ix) = Lx.' * res;
b1(Iy) = Ly.' * res;
b1(Iz) = Lz.' * res;

M1 = zeros(3 * Ndof, 3 * Ndof);
M1(Ix,Ix) = Lx.' * Lx;
M1(Ix,Iy) = Lx.' * Ly;
M1(Ix,Iz) = Lx.' * Lz;

M1(Iy,Ix) = Ly.' * Lx;
M1(Iy,Iy) = Ly.' * Ly;
M1(Iy,Iz) = Ly.' * Lz;

M1(Iz,Ix) = Lz.' * Lx;
M1(Iz,Iy) = Lz.' * Ly;
M1(Iz,Iz) = Lz.' * Lz;
toc

% Speed testing the C++ code
% -------------------------------

% first run of a newly compiled code is not cached yet
[M2, b2, r] = localDVC_kernel(f, g, a, C, L, threads);

% runtime without computing the residual image
tic
[M2, b2, r] = localDVC_kernel(f, g, a, C, L, threads);
toc

% runtime with computing the residual image
tic
[M2, b2, r, R] = localDVC_kernel(f, g, a, C, L, threads);
toc

% Example DVC code
% -------------------------------

maxit = 10;
convcrit = 1e-4;

rm = Inf;
for it = 1:maxit
    % update M and b for current a
    [M, b, r] = localDVC_kernel(f, g, a, C, L, threads);

    % compute the update in a
    da = M \ b;

    % update the degrees of freedom
    a = a + da;

    % some output to the screen
    dr = r - rm;
    rm = r;
    fprintf('it:%3d, r:%10.3e, dr:%10.3e, |da|:%10.3e\n',it,r,dr,rms(da));

    % convergence test
    if rms(da) < convcrit
        fprintf('converged\n\n');
        break
    end
end

disp(a(1:3))
