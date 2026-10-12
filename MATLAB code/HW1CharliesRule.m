%% Numerical PDEs HW1
% Problem 3: Charlies Rule
clear;clc;close all;

% manufactured solns
U_exact = @(x,y) sin(2*pi*y).*sin(2*pi*x);
V_exact = @(x,y) -3.5 + cos(2*pi*x).*(cos(2*pi*y)-1);
P_exact = @(x,y) sin(2*pi*y).*cos(2*pi*x);

% RHS
f = @(x,y) (2*pi - 8*pi*pi)*sin(2*pi*y)*sin(2*pi*x);
g = @(x,y) 4*pi*pi*cos(2*pi*x) - (8*pi*pi + 2*pi)*cos(2*pi*x)*cos(2*pi*y);

% params
a = 0;
b = 1;

%vector of different step sizes or grid points
NV = 2.^(5:5); 
hV = zeros(size(NV)); 

bslashTime = zeros(size(NV));
gmresTime = zeros(size(NV));

for l=1:length(NV)
    N=NV(l);
    h = (b-a)/N;
    hV(l) = h;
    x = linspace(a,b,2*N+1);
    y = linspace(a,b,2*N+1);
    z = a + (0:N-1)*h;          % Integer points: 0, h, ..., 1-h
    half = a + ((0:N-1)+0.5)*h; % Half points: h/2, 3h/2, ..., 1-h/2
    yVj = a + (1:N-1)*h; 


% Number of Unknowns
n = 3*(N^2) - N;
nU = N*N; 
nV = N*(N-1); 
nP = nU;

% matrix index to flat index
Uidx = @(i,j) i + (j-1)*N;
Vidx = @(i,j) nU + i + (j-1)*N;
Pidx = @(i,j) nU + nV + i + (j-1)*N;

% sparse matrix set up
row = [];
col = [];
val = [];

% U parts at i, j+1/2
for i=1:N
    for j=1:N
        k=Uidx(i,j);

        %center
        row(end+1) = k;
        col(end+1) = k;
        if j==1 || j==N
            val(end+1) = -5;
        else
            val(end+1) = -4;
        end
        

        % left neighbor
        if i > 1
            row(end+1) = k;
            col(end+1) = Uidx(i-1,j);
            val(end+1) = 1;
        else
            row(end+1) = k;
            col(end+1) = Uidx(N,j);
            val(end+1) = 1;
        end

        % right neighbor
        if i < N
            row(end+1) = k;
            col(end+1) = Uidx(i+1,j);
            val(end+1) = 1;
        else
            row(end+1) = k;
            col(end+1) = Uidx(1,j);
            val(end+1) = 1;
        end

        % below neighor
        if j > 1
            row(end+1) = k;
            col(end+1) = Uidx(i,j-1);
            val(end+1) = 1;
        end

        % above neighbor
        if j < N
            row(end+1) = k;
            col(end+1) = Uidx(i,j+1);
            val(end+1) = 1;
        end

        %Pressure
            row(end+1) = k;
            col(end+1) = Pidx(i,j);
            val(end+1) = -h;
       
        if i > 1
            row(end+1) = k;
            col(end+1) = Pidx(i-1,j);
            val(end+1) = h;
        else 
            row(end+1) = k;
            col(end+1) = Pidx(N,j);
            val(end+1) = h;
        end
    end
end

% V parts at i+1/2, j
for i=1:N
    for j=1:N-1
        k=Vidx(i,j);

        %center
        row(end+1) = k;
        col(end+1) = k;
        val(end+1) = -4;

        % left neighbor
        if i > 1
            row(end+1) = k;
            col(end+1) = Vidx(i-1,j);
            val(end+1) = 1;
        else
            row(end+1) = k;
            col(end+1) = Vidx(N,j);
            val(end+1) = 1;
        end

        % right neighbor
        if i < N
            row(end+1) = k;
            col(end+1) = Vidx(i+1,j);
            val(end+1) = 1;
        else
            row(end+1) = k;
            col(end+1) = Vidx(1,j);
            val(end+1) = 1;
        end

        % below neighbor
        if j > 1
            row(end+1) = k;
            col(end+1) = Vidx(i,j-1);
            val(end+1) = 1;
        end

        % above neighbor
        if j < N-1
            row(end+1) = k;
            col(end+1) = Vidx(i,j+1);
            val(end+1) = 1;
        end

        %pressure
       row(end+1)=k; col(end+1)=Pidx(i,j+1); val(end+1)= -h;
       row(end+1)=k; col(end+1)=Pidx(i,j);   val(end+1)=h;
    end
end

% RHS for U: integer x, half-step y
RHS = zeros(n,1);
for i = 1:N
    for j = 1:N
        RHS(Uidx(i,j)) = h^2 * f(x(2*i-1), y(2*j));
    end
end

% RHS for V: half-step x, integer y
for i = 1:N
    for j = 1:N-1
        RHS(Vidx(i,j)) = h^2 * g(x(2*i), y(2*j+1));
    end
end

% Divergence equation part
for i = 1:N
    for j = 1:N

        k = Pidx(i,j);

        % +U(i+1,j), periodic wrap
        row(end+1) = k;
        if i < N 
            col(end+1) = Uidx(i+1,j); 
        else
            col(end+1) = Uidx(1,j); 
        end
        val(end+1) = 1;


        % -U(i,j)
        row(end+1) = k; 
        col(end+1) = Uidx(i,j); 
        val(end+1) = -1;

        % V(i,j)
        if j < N
            row(end+1) = k;
            col(end+1) = Vidx(i,j);
            val(end+1) = 1;
        end

        % -V(i,j-1)
        if j > 1
            row(end+1) = k;
            col(end+1) = Vidx(i,j-1);
            val(end+1) = -1;
        end

    end
end
       
A = sparse(row,col,val,n,n);

%% Boundary Conditions

% U = 0 on the top and bottom: Do Nothing

% V = -3.5 on top and bottom
for i = 1:N
    RHS(Vidx(i,1)) = RHS(Vidx(i,1)) + 3.5;
    RHS(Vidx(i,N-1)) = RHS(Vidx(i,N-1)) + 3.5;
end

% Divergence boundary contributions

for i = 1:N

    % bottom: V(i,1) = -3.5
    RHS(Pidx(i,1)) = RHS(Pidx(i,1)) - 3.5;

    % top: V(i,N+1) = -3.5
    RHS(Pidx(i,N)) = RHS(Pidx(i,N)) + 3.5;

end

kp = Pidx(1,1);
A(kp,:) = 0;  A(kp,kp) = 1;  RHS(kp) = 0;   % pin P(1,1)


% %-------(section 1)-------------------------------------------------------------------
% %% Timing
% tic
% u = A\RHS;
% 
% elapsedTime = toc;
% bslashTime(l) = elapsedTime;
% 
% %% GMRES
% % set of A has h^2 error and h is abt .01 
% % tol should be .1*.01*.01 
% tic
% tol = .00001;
% max_its = n;
% ugm = gmres(A,RHS,[],tol,max_its);
% 
% elaspedTime = toc;
% gmresTime(l) = elaspedTime;
% 
% end
% 
% %% Timing Plot
% figure;
% plot(NV,bslashTime, '-o', 'Color', [.267,.467,.667], 'LineWidth', 2.5, 'MarkerFaceColor', [.267,.467,.667]);
% hold on;
% plot(NV,gmresTime, '-s', 'Color', [.933,.4,.467], 'LineWidth', 2.5, 'MarkerFaceColor', [.933,.4,.467]);
% hold off;
% 
% xlabel('Gridspaces (N)');
% ylabel('Time (sec)');
% legend('Backslash', 'gmres');
% grid on

%------------(section 2)------------------------------------------------------------
%% Changing Tol

%% GMRES
% set of A has h^2 error and h is abt .01 
% tol should be .1*.01*.01 
tic
tol = .00001;
badtol = .1;
max_its = n;
u = gmres(A,RHS,[],tol,max_its);

%% Extract numerical U and V

% U has N x N unknowns
U = zeros(N,N);

for i = 1:N
    for j = 1:N
        U(i,j) = u(Uidx(i,j));
    end
end

% V has N x (N-1) unknowns
V = zeros(N,N+1);

for i = 1:N
    for j = 1:N-1
        V(i,j+1) = u(Vidx(i,j)); 
    end
end

% Boundary conditions for V
V(:,1) = -3.5;       % bottom
V(:,N+1) = -3.5;     % top


%% Interpolate U and V to cell centers

Uc = zeros(N,N);
Vc = zeros(N,N);

for i = 1:N
    for j = 1:N

        % U is staggered in x
        if i < N
            Uc(i,j) = 0.5*(U(i,j) + U(i+1,j));
        else
            % periodic x boundary
            Uc(i,j) = 0.5*(U(i,j) + U(1,j));
        end

        % V is staggered in y
        Vc(i,j) = 0.5*(V(i,j) + V(i,j+1));

    end
end

%% Velocity magnitude
umag = sqrt(Uc.^2 + Vc.^2);

%% Plot velocity magnitude
xc = a + h/2 : h : b-h/2;
yc = a + h/2 : h : b-h/2;
[Xc,Yc] = meshgrid(xc,yc);

Mplot = umag';
Uplot = Uc';
Vplot = Vc';

figure;
surf(Xc, Yc, zeros(size(Xc)), Mplot);   % flat surface, colored by |u|
shading interp;
view(2);
axis equal tight;
colorbar;
colormap summer
hold on;

hs = streamslice(Xc, Yc, Uplot, Vplot);
set(hs, 'Color', 'k', 'LineWidth', 0.8);
for k = 1:numel(hs)
    hs(k).ZData = ones(size(hs(k).XData));   % lift above the surface at z=0
end

xlabel('x');
ylabel('y');
title('||u||');
hold off;

end

