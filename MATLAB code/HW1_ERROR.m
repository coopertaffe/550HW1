%% Numerical PDEs HW1
% Problem 2: Error Analysis
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
NV = 2.^(3:7); 
hV = zeros(size(NV)); 

Err_U = zeros(size(NV));
Err_V = zeros(size(NV));
Err_P = zeros(size(NV));

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


%% EXACT
[Xu,Yu] = ndgrid(z,half);
[Xv,Yv] = ndgrid(half,yVj);
[Xp,Yp] = ndgrid(half,half);
exactU = U_exact(Xu,Yu);
exactV = V_exact(Xv,Yv);
exactP = P_exact(Xp,Yp);
       
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

u = A\RHS;

%% Extract numerical U and V
U = u(1:nU); %size = NxN
V = u(nU+1:nU+nV); %size = NxN-1
P = u(nU+nV+1:end); %size = NXN 
P = P - mean(P);

%% EXACT
[Xu,Yu] = ndgrid(z,half);
[Xv,Yv] = ndgrid(half,yVj);
[Xp,Yp] = ndgrid(half,half);
exactU = U_exact(Xu,Yu);
exactV = V_exact(Xv,Yv);
exactP = P_exact(Xp,Yp);

Udiff = exactU(:) - U;
Vdiff = exactV(:) - V;
Pdiff = exactP(:) - P;

Err_U(l) = norm(Udiff)/norm(exactU(:));
Err_V(l) = norm(Vdiff)/norm(exactV(:));
Err_P(l) = norm(Pdiff)/norm(exactP(:));


end

%plotting
figure
loglog(NV, Err_U, '-o', 'Color', [.267,.467,.667], 'LineWidth', 2.5, 'MarkerFaceColor', [.267,.467,.667]) %blue
hold on
loglog(NV, Err_V, '-s', 'Color', [.933,.4,.467], 'LineWidth', 2.5, 'MarkerFaceColor', [.933,.4,.467]) %pink
loglog(NV, Err_P, '-^', 'Color', [.8,.733,.267], 'LineWidth', 2.5, 'MarkerFaceColor', [.8,.733,.267]) %yellowgreen
shiftFactor = 0.5;   % how far below the curve to place the line (smaller = lower)

% Reference slope lines
ref1 = NV.^(-1);
ref2 = NV.^(-2);

% Scale reference lines to appear near the error curves
ref1 = 1.5 * ref1 * Err_U(1)/ref1(1);
ref2 = 1.5 * ref2 * Err_U(1)/ref2(1);

loglog(NV, ref1, 'k:', 'LineWidth', 2)
loglog(NV, ref2, 'k-.', 'LineWidth', 2)

hold off
xlabel('Grid Size (N)')
ylabel('Relative Error')
legend('Err_U','Err_V','Err_P','Ref Slope: -1', 'Ref Slope: -2')
grid on