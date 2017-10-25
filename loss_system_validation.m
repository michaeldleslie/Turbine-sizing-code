%% Initialization
clearvars;
close all;
% Add lib to path
addpath(genpath(pwd));

%% Validate Y_p1
clearvars

f(1) = figure(1);
sc = 0.1:0.01:1.2;
alpha = [10 15 20 25 30 40 50];
for i = 1:length(alpha)
    alpha2 = alpha(i);
    yp1 = Y_p1_nozzle_profile_loss(alpha2, sc);
    plot(sc, yp1, sc, yp1, 'o');
    hold on;
end
hold off;
xlabel('s/c');
ylabel('Yp1');
grid;

%% Validate Y_p2
clearvars

f(2) = figure(2);
sc = 0.1:0.01:1.1;
alf = [10 15 20 25 30 35 40 50];
for i = 1:length(alf)
    alfa2 = alf(i);
    yp2 = Y_p2_impulse_airfoil_profile_loss(alfa2, sc);
    plot(sc, yp2, sc, yp2, 'o');
    hold on;
end
hold off;
xlabel('s/c');
ylabel('Yp2');
grid;

%% Validate K_p
clearvars

f(3) = figure(3);
M1 = 0.1:0.01:1.2;
M2 = 0.1:0.01:1.2;
k = 1;
KP = zeros(length(M1), length(M2));
for i1 = 1:length(M1)
    for i2 = 1:length(M2)
        KP(i1,i2) = K_p_compressibility_correction(M1(i1),M2(i2));
    end
end
colormap jet
pcolor(M1, M2, KP);
shading('interp');
colorbar;
xlabel('M1');
ylabel('M2'); 

%% Validate K_Re
clearvars

f(4) = figure(4);
roughness = 5E-6;
rho = 1.225;
c =[0.01 .1 1 10];
niu = 1.9e-6/1.225;
Re = 10.^(4:.05:9);
k = 1;
for i1 = 1:length(c)
    C2 = Re*niu./c(i1);
    for i2 = 1:length(C2)
        KRe(k) = K_Re_reynolds_correction(rho, C2(i2), c(i1), niu*rho, roughness);
        Re_c(k) = C2(i2)*c(i1)/niu;
        k = k+1;
    end
end

semilogx(Re_c,KRe,'.');
xlabel('Re_c');
ylabel('KRe');
grid;

%% Validate Y_S
clearvars

%figure(100)
alfa1 = (90-(20:5:40));
alfa2 = -(90-(30:5:50));
beta1 = alfa1;
Lc = 1:0.25:5;
k = 1;
ys = zeros(length(alfa1)*length(alfa2)*length(beta1)*length(Lc), 1);
hic = zeros(length(alfa1)*length(alfa2)*length(beta1)*length(Lc), 1);
for i1 = 1:length(alfa1)
    for i2 = 1:length(alfa2)
        for i3 = 1:length(beta1)
            %for i4 = 1:length(sc);
            for i5 = 1:length(Lc)
                a1 = alfa1(i1);
                a2 = alfa2(i2);
                b1 = beta1(i3);
                %scc = sc(i4);
                Lcc = Lc(i5);
                ys(k) = Y_S_secondary_flow_loss(a1,a2,b1,Lcc,0.706);
                hic(k) = Lcc;
                %figure(1);
                %plot(hic(k),ys(k),'.');
                %title(num2str(k));
                %hold on;
                %pause(0.001);
                k = k+1;
            end
            %end;
        end
    end
end
f(5) = figure(5);
plot(hic,ys,'o');
xlabel('L/c');
ylabel('Ys');

%% Validate Y_cl
clearvars

alfa1 = +(90-(20:5:40));
alfa2 = -(90-(30:5: 50));
Lc = [1.5 3 5];
dL = 0.002:0.0005:0.012;
s_c_design = 0.706;

y_nshr = zeros(length(alfa1)*length(alfa2)*length(dL)*length(Lc), 1);
y_wshr = zeros(length(alfa1)*length(alfa2)*length(dL)*length(Lc), 1);
dLk = zeros(length(alfa1)*length(alfa2)*length(dL)*length(Lc), 1);

k = 1;
for i1 = 1:length(alfa1)
    for i2 = 1:length(alfa2)
        %for i3 = 1:length(beta1),
        for i4 = 1:length(dL)
            for i5 = 1:length(Lc)
                a1 = alfa1(i1);
                a2 = alfa2(i2);
                Lcc = Lc(i5);
                dLL = dL(i4);
                y_nshr(k) = Y_cl_clearance_loss(a1,a2,s_c_design,Lcc,dLL,0);
                y_wshr(k) = Y_cl_clearance_loss(a1,a2,s_c_design,Lcc,dLL,1);
                dLk(k) = dLL;
%                 figure(1);
%                 plot(dLk(k),y_nshr(k),'.',dLk(k),y_wshr(k),'.');
%                 title(num2str(k));
%                 hold on;
%                 pause(0.001);
                k = k+1;
            end
        end
        %end;
    end
end

f(6) = figure(6);
plot(dLk,y_nshr,'.',dLk,y_wshr,'.')
grid on
xlabel('Clearance \delta/L');
ylabel('Y_{cl}');