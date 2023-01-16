clear; close all

if true
    % debug = '-g';
    debug = '';

    cppfile = 'localDVC_kernel.cpp';
    mex(cppfile,debug,'-largeArrayDims','-DNOMINMAX')

end


addpath(correli_pathdef)

% Define the path to the images
imagepath = fullfile('virtimage3D');
imagefiles = dir(fullfile(imagepath, '*.raw'));
imagefiles = {imagefiles.name}';
images = strcat(imagepath, filesep, imagefiles);

refimage = images{1};
defimage = images{2};

if true
    % Load the images
    f = raw_read(refimage);
    g = raw_read(defimage);

    % convert to floating point
    f = single(f);
    g = single(g);

else
    siz = [160, 160, 160];

    f = zeros(siz,'single');
    g = ones(siz,'single');

%     x = single(linspace(-1,1,siz(2)));
%     f = repmat(x, siz(1), 1, siz(3));
%     y = single(linspace(-1,1,siz(1)));
%     f = repmat(y(:), 1, siz(2), siz(3));
    z = single(linspace(-1,1,siz(3)));
    f = repmat(permute(z, [3,1,2]), siz(1), siz(2), 1);
    g = f + 10;
end


siz = size(f);

% f = f ./ 255;
% g = g ./ 255;

% -------------------------------
mar = 10;
x = mesh_linspace(1+mar, siz(2)-mar, 2);
y = mesh_linspace(1+mar, siz(1)-mar, 2);
z = mesh_linspace(1+mar, siz(3)-mar, 2);
Mesh = mesh_gen_structured(x,y,z,'T4');

tic
init = zeros(8,3);
[cor.a, cor.r, cor.R, cor.stat] = correlate(f,g,Mesh,init);
toc


% -------------------------------

C = round(siz./2) + [3,4,5];
L = 51 * [1, 1, 1];
threads = 0;

Ndof = 4;
a = zeros( 3 * Ndof,1);
% a(1:3) = mean(cor.a);

Ix = 1:3:3*Ndof;
Iy = 2:3:3*Ndof;
Iz = 3:3:3*Ndof;

% -------------------------------
tic
x = (C(1) - ((L(1) - 1) / 2) + (0:L(1)-1));
y = (C(2) - ((L(2) - 1) / 2) + (0:L(2)-1));
z = (C(3) - ((L(3) - 1) / 2) + (0:L(3)-1));

[X, Y, Z] = meshgrid(1:siz(2),1:siz(1),1:siz(3));

% Xn = X - C(1);
% Yn = Y - C(2);
% Zn = Z - C(3);
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
ortho_view(R);

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

disp([ M1, b1] );
disp(rms(res));
figure;
imagesc(log10(abs(M1)))
colorbar

% fprintf('%12d, [%3d,%3d,%3d], %6.1f, %6.1f, [%6.1f, %6.1f, %6.1f]\n',[find(roi)-1, X(roi)-1, Y(roi)-1, Z(roi)-1, f(roi), R(roi), fx(roi), fy(roi), fz(roi)].');

%% -------------------------------



tic
[M2, b2, r] = localDVC_kernel(f, g, a, C, L, threads);
toc

tic
[M2, b2, r] = localDVC_kernel(f, g, a, C, L, threads);
toc

tic
[M2, b2, r, R] = localDVC_kernel(f, g, a, C, L, threads);
toc

ortho_view(R)

disp([M2, b2]);
disp(r);
figure;
imagesc(log10(abs(M2)))
colorbar

%% -----------------------------------


maxit = 10;
convcrit = 1e-4;

rm = Inf;
for it = 1:maxit

    [M, b, r] = localDVC_kernel(f, g, a, C, L, threads);

    da = M \ b;

    a = a + da;

    dr = r - rm;
    rm = r;

    fprintf('it:%3d, r:%10.3e, dr:%10.3e, |da|:%10.3e\n',it,r,dr,rms(da));


    if rms(da) < convcrit
        fprintf('converged\n\n');
        break
    end


end


[M, b, r, R] = localDVC_kernel(f, g, a, C, L, threads);
ortho_view(R)
