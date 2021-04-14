clc
clear
close all

%%
%Parametri
Cp = 4181.3; %toplotni kapacitet [J/kgK]
ro = 1000; %gustina vode [kg/m^3]
Cdh =0.65; %koeficijent praznjenja ventila za hladnu vodu
Cdt = 0.65; %koeficijent praznjenja ventila za toplu vodu
Cd1 = 0.65; %koeficijent praznjenja ventila za rezervoar 1
Cd2 = 0.65; %koeficijent praznjenja ventila za rezervoar 2
wt = 0.03; %gradijent povrsine ventila [m^2/m]
wh = 0.03; 
w1 = 0.1;
w2 = 0.1;
theta_h = 10 + 273; % temperatura hladne vode [K]
theta_t = 70 + 273; % temperatura tople vode [K]
Kvh = 1/1000; % konstanta proporcionalnosti [m/v]
Kvt = 1/1000;
Kv1 = 1.5/1000;
Kv2 = 1.5/1000;
Ar1 = 0.04; %povrsina poprecnog preseka rezervoara 1 [m^2]
Ar2 = 0.04; %povrsina poprecnog preseka rezervoara 2 [m^2]
H1_max = 1.5; %visina rezervoara 1 [m]
H1_min = 0;
H2_max = 1.5; %visina rezervoara 2 [m]
H2_min = 0;
H1_N = 0.7; %nominalni nivo u rezervoaru 1 [m]
H2_N = 0.8; %nominalni nivo u rezervoaru 2 [m]
g = 9.81; %ubrzanje sile zemljine teze [m/s^2]
Uv1 = 3; %upravljacki napon za ventil hladne vode [V]
Ph = 3.25*10^5; %natpritisak hladne vode [Pa]
Pt_N = 3.25*10^5; %nominalna vrednost natpritiska tople vode [Pa]
theta_m_N = 30 + 273; %nominalna temperatura mesavine, tj. nominalna temperatura na izlazu iz rezervoara 1 [K]
Uvt_min = 0; %minimalna vrednost upravljackog napona za ventil tople vode [V]
Uvt_max = 10; %maksimalna vrednost upravljackog napona za ventil tople vode [V]
Uvh_min = 0; 
Uvh_max = 10;
Uv2_min = 0; 
Uv2_max = 10;
%%
%Nominalne vrednosti upravljanja
Uvt_N = 0;
syms Ut_N Uh_N U2_N 
disp('Nominalne vrednosti upravljanja:')
%jednacine za rezervoar 1 u nominalnom rezimu (dH1_N/dt = 0, dtheta_m_N/dt = 0)
e1 = Cdh*wh*Kvh*Uh_N*sqrt(2*Ph/ro) + Cdt*wt*Kvt*Ut_N*sqrt(2*Pt_N/ro) - Cd1*w1*Kv1*Uv1*sqrt(2*g*H1_N);
e2 = Cdh*wh*Kvh*Uh_N*sqrt(2*Ph/ro)*(theta_h - theta_m_N) + Cdt*wt*Kvt*Ut_N*sqrt(2*Pt_N/ro)*(theta_t - theta_m_N);
sol1 = solve(e1 == 0, Uh_N);
sol2 = solve(e2 == 0, Uh_N);
Uvt_N = double(solve(sol1 == sol2, Uvt_N))
f1n = subs(e1, Ut_N, Uvt_N);
Uvh_N = double(solve(f1n == 0, Uh_N))

%jednacina za rezervoar 2 u nominalnom rezimu (dH2_N/dt = 0)
e3 = Cd1*w1*Kv1*Uv1*sqrt(2*g*H1_N) - Cd2*w2*Kv2*U2_N*sqrt(2*g*H2_N);
Uv2_N = double(solve(e3 == 0, U2_N))
%%
%Linearizacija
syms Uvt Uvh Uv2 Pt H1 H2 theta_m uvt uvh uv2 pt h1 h2 t_m
%Uvt, Uvh, Uv2, Pt, H1, H2, theta_m: vrednosti promenljivih u totalnim koordinatama
%uvt, uvh, uv2, pt, h1, h2, t_m: odstupanja
%X = x + X_N => totalna koordinata je zbir nominalne vrednosti i odstupanja

%Jednacina kontinuiteta za rezervoar 1
f1 = (Cdh*wh*Kvh*Uvh*sqrt(2*Ph/ro) + Cdt*wt*Kvt*Uvt*sqrt(2*Pt/ro) - Cd1*w1*Kv1*Uv1*sqrt(2*g*H1))/Ar1;
%izvod po Uvt
df1_Uvt = diff(f1, Uvt, 1);
df1_Uvt = double(subs(df1_Uvt, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po Uvh
df1_Uvh = diff(f1, Uvh, 1);
df1_Uvh = double(subs(df1_Uvh, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po H1
df1_H1 = diff(f1, H1, 1);
df1_H1 = double(subs(df1_H1, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po Uv2
df1_Uv2 = diff(f1, Uv2, 1);
df1_Uv2 = double(subs(df1_Uv2, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po H2
df1_H2 = diff(f1, H2, 1);
df1_H2 = double(subs(df1_H2, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po theta_m
df1_theta_m = diff(f1, theta_m, 1);
df1_theta_m = double(subs(df1_theta_m, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po Pt
df1_Pt = diff(f1, Pt, 1);
df1_Pt = double(subs(df1_Pt, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%linearizovana jednacina
disp('Linearizovana jednacina kontinuiteta za rezervoar 1:')
dh1 = df1_Uvt*uvt + df1_Uvh*uvh + df1_Uv2*uv2 + df1_H1*h1 + df1_H2*h2 + df1_Pt*pt + df1_theta_m*t_m;
dh1 = vpa(dh1, 4)

%Jednacina kontinuiteta za rezervoar 2
f2 = (Cd1*w1*Kv1*Uv1*sqrt(2*g*H1) - Cd2*w2*Kv2*Uv2*sqrt(2*g*H2))/Ar2;
%izvod po Uvt
df2_Uvt = diff(f2, Uvt, 1);
df2_Uvt = double(subs(df2_Uvt, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po Uvh
df2_Uvh = diff(f2, Uvh, 1);
df2_Uvh = double(subs(df2_Uvh, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po H1
df2_H1 = diff(f2, H1, 1);
df2_H1 = double(subs(df2_H1, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po Uv2
df2_Uv2 = diff(f2, Uv2, 1);
df2_Uv2 = double(subs(df2_Uv2, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po H2
df2_H2 = diff(f2, H2, 1);
df2_H2 = double(subs(df2_H2, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po theta_m
df2_theta_m = diff(f2, theta_m, 1);
df2_theta_m = double(subs(df2_theta_m, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po Pt
df2_Pt= diff(f2, Pt, 1);
df2_Pt = double(subs(df2_Pt, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%linearizovana jednacina
disp('Linearizovana jednacina kontinuiteta za rezervoar 2:')
dh2 = df2_Uvt*uvt + df2_Uvh*uvh + df2_Uv2*uv2 + df2_H1*h1 + df2_H2*h2 + df2_Pt*pt + df2_theta_m*t_m;
dh2 = vpa(dh2, 4)

%Jednacina energije za rezervoar 1
f3 = (Cdh*wh*Kvh*Uvh*sqrt(2*Ph/ro)*(theta_h - theta_m) + Cdt*wt*Kvt*Uvt*sqrt(2*Pt/ro)*(theta_t - theta_m))/(Ar1*H1);
%izvod po Uvt
df3_Uvt = diff(f3, Uvt, 1);
df3_Uvt = double(subs(df3_Uvt, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po Uvh
df3_Uvh = diff(f3, Uvh, 1);
df3_Uvh = double(subs(df3_Uvh, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po H1
df3_H1 = diff(f3, H1, 1);
df3_H1 = double(subs(df3_H1, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po Uv2
df3_Uv2 = diff(f3, Uv2, 1);
df3_Uv2 = double(subs(df3_Uv2, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po H2
df3_H2 = diff(f3, H2, 1);
df3_H2 = double(subs(df3_H2, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po theta_m
df3_theta_m = diff(f3, theta_m, 1);
df3_theta_m = double(subs(df3_theta_m, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%izvod po Pt
df3_Pt= diff(f3, Pt, 1);
df3_Pt = double(subs(df3_Pt, [Uvt Uvh Uv2 Pt H1 H2 theta_m], [Uvt_N Uvh_N Uv2_N Pt_N H1_N H2_N theta_m_N]));
%linearizovana jednacina
disp('Linearizovana jednacina energije za rezervoar 1:')
%u ovoj jednacini se javlja izvod po H1 sa koeficijentom velicine reda
%10^-72 pa se moze zanemariti
dh3 = df3_Uvt*uvt + df3_Uvh*uvh + df3_Uv2*uv2 + df3_H1*h1 + df3_H2*h2 + df3_Pt*pt + df3_theta_m*t_m;
dh3 = vpa(dh3, 4)
%%
%Kalmanov regulator
%dx = Ax + Bu

%vektor upravljanja je: u = [uvh; uvt; uv2], jer su to velicine kojima se
%kontrolisu ventili, pt je poremecaj a ne upravljacka velicina
%dimenzije matrica su onda: A(3x3), B(3x3)

Ak = [df1_H1 df1_H2 df1_theta_m; df2_H1 df2_H2 df2_theta_m; df3_H1 df3_H2 df3_theta_m];
Bk = [df1_Uvh df1_Uvt df1_Uv2; df2_Uvh df2_Uvt df2_Uv2; df3_Uvh df3_Uvt df3_Uv2];
Ck = [1 0 0; 0 1 0; 0 0 1];
Dk = [0 0 0; 0 0 0; 0 0 0];
[K S E] = lqr(Ak, Bk, eye(3), eye(3));
disp('Pojacanje u povratnoj grani:')
K = K 
%%
%Poredjenje odziva nelinearnog i linearnog modela

load('R2_nivo.mat')
load('R1_temperatura.mat')
load('R1_nivo.mat')

%nivo u rezervoaru 1
figure(1)
plot(H1(1, :), H1(2, :), '-', H1(1, :), H1(3, :), '-.', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('H1(t) [m], h1(t) [m]', 'fontsize', 15)
legend('nelinearni', 'linearni', 'location', 'bestoutside', 'fontsize', 15)

%nivo u rezervoaru 2
figure(2)
plot(H2(1, :), H2(2, :), '-', H2(1, :), H2(3, :), '-.', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('H2(t) [m], h2(t) [m]', 'fontsize', 15)
legend('nelinearni', 'linearni', 'location', 'bestoutside', 'fontsize', 15)

%temperatura u rezervoaru 1
figure(3)
plot(theta(1, :), theta(2, :), '-', theta(1, :), theta(3, :), '-.', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('\Theta_m(t) [\circC], \theta_m(t) [\circC]', 'fontsize', 15)
legend('nelinearni', 'linearni', 'location', 'bestoutside', 'fontsize', 15)
%%
%Prikaz konvergencije stanja ka nultom ravnoteznom stanju pri primeni
%Kalmanovog regulatora

load('kalmanov_regulator.mat')
figure(4)
plot(k(1, :), k(2, :), '-', k(1, :), k(3, :), '-.', k(1, :), k(4, :), '--', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('h1(t) [m], h2(t) [m] , \theta_m(t) [m]', 'fontsize', 15)
legend('h1(t)', 'h2(t)', '\theta_m(t)', 'location', 'bestoutside', 'fontsize', 15)






