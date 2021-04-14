clc
clear
close all
%Izbacen je cevovod tople vode
%%
%Parametri
Cp = 4181.3; %toplotni kapacitet [J/kgK]
ro = 1000; %gustina vode [kg/m^3]
Cdh =0.65; %koeficijent praznjenja ventila za hladnu vodu
Cd1 = 0.65; %koeficijent praznjenja ventila za rezervoar 1
Cd2 = 0.65; %koeficijent praznjenja ventila za rezervoar 2
wh = 0.03; 
w1 = 0.1;
w2 = 0.1;
theta_h = 10 + 273; % temperatura hladne vode [K]
Kvh = 1/1000; % konstanta proporcionalnosti [m/v]
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
Uvh_min = 0; 
Uvh_max = 10;
Uv2_min = 0; 
Uv2_max = 10;
%%
%Nominalne vrednosti upravljanja
syms Uh_N U2_N 
disp('Nominalne vrednosti upravljanja:')
%jednacina za rezervoar 1 u nominalnom rezimu (dH1_N/dt = 0)
e1 = Cdh*wh*Kvh*Uh_N*sqrt(2*Ph/ro) - Cd1*w1*Kv1*Uv1*sqrt(2*g*H1_N);
Uvh_N = double(solve(e1 == 0, Uh_N))
%jednacina za rezervoar 2 u nominalnom rezimu (dH2_N/dt = 0)
e2 = Cd1*w1*Kv1*Uv1*sqrt(2*g*H1_N) - Cd2*w2*Kv2*U2_N*sqrt(2*g*H2_N);
Uv2_N = double(solve(e2 == 0, U2_N))
%%
%Linearizacija
clc

syms Uvh Uv2 H1 H2 uvh uv2 h1 h2 
%Uvh, Uv2, H1, H2: vrednosti promenljivih u totalnim koordinatama
%uvh, uv2, h1, h2: odstupanja
%X = x + X_N => totalna koordinata je zbir nominalne vrednosti i odstupanja

%Jednacina kontinuiteta za rezervoar 1
f1 = (Cdh*wh*Kvh*Uvh*sqrt(2*Ph/ro) - Cd1*w1*Kv1*Uv1*sqrt(2*g*H1))/Ar1;
%izvod po Uvh
df1_Uvh = diff(f1, Uvh, 1);
df1_Uvh = double(subs(df1_Uvh, [Uvh Uv2 H1 H2], [Uvh_N Uv2_N H1_N H2_N]));
%izvod po H1
df1_H1 = diff(f1, H1, 1);
df1_H1 = double(subs(df1_H1, [Uvh Uv2 H1 H2], [Uvh_N Uv2_N H1_N H2_N]));
%izvod po Uv2
df1_Uv2 = diff(f1, Uv2, 1);
df1_Uv2 = double(subs(df1_Uv2, [Uvh Uv2 H1 H2], [Uvh_N Uv2_N H1_N H2_N]));
%izvod po H2
df1_H2 = diff(f1, H2, 1);
df1_H2 = double(subs(df1_H2, [Uvh Uv2 H1 H2], [Uvh_N Uv2_N H1_N H2_N]));
%linearizovana jednacina
disp('Linearizovana jednacina kontinuiteta za rezervoar 1:')
dh1 = df1_Uvh*uvh + df1_Uv2*uv2 + df1_H1*h1 + df1_H2*h2;
dh1 = vpa(dh1, 4)

%Jednacina kontinuiteta za rezervoar 2
f2 = (Cd1*w1*Kv1*Uv1*sqrt(2*g*H1) - Cd2*w2*Kv2*Uv2*sqrt(2*g*H2))/Ar2;
%izvod po Uvh
df2_Uvh = diff(f2, Uvh, 1);
df2_Uvh = double(subs(df2_Uvh, [Uvh Uv2 H1 H2], [Uvh_N Uv2_N H1_N H2_N]));
%izvod po H1
df2_H1 = diff(f2, H1, 1);
df2_H1 = double(subs(df2_H1, [Uvh Uv2 H1 H2], [Uvh_N Uv2_N H1_N H2_N]));
%izvod po Uv2
df2_Uv2 = diff(f2, Uv2, 1);
df2_Uv2 = double(subs(df2_Uv2, [Uvh Uv2 H1 H2], [Uvh_N Uv2_N H1_N H2_N]));
%izvod po H2
df2_H2 = diff(f2, H2, 1);
df2_H2 = double(subs(df2_H2, [Uvh Uv2 H1 H2], [Uvh_N Uv2_N H1_N H2_N]));
%linearizovana jednacina
disp('Linearizovana jednacina kontinuiteta za rezervoar 2:')
dh2 = df2_Uvh*uvh + df2_Uv2*uv2 + df2_H1*h1 + df2_H2*h2;
dh2 = vpa(dh2, 4)
%%
%Matematicki model u prostoru stanja
clc
A = [df1_H1 df1_H2; df2_H1 df2_H2];
B = [df1_Uvh df1_Uv2; df2_Uvh df2_Uv2];
C = [1 0; 0 1];
D = [0 0; 0 0];

%Prenosna matrica
sys = ss(A, B, C, D);
[num1 den1] = ss2tf(A, B, C, D, 1);
W11 = tf(num1(1, :), den1); %h1 u odnosu na uvh
W21 = tf(num1(2, :), den1); %h2 u odnosu na uvh
[num2 den2] = ss2tf(A, B, C, D, 2);
W12 = tf(num2(1, :), den2); %h1 u odnosu na uv2
W22 = tf(num2(2, :), den2); %h2 u odnosu na uv2
W_o = [W11 W12; W21 W22]
%%
%Dekuplovanje
clc

%Trazimo W_dekuplera za visestruko prenosni sistem, to nam
%omogucava da dobijemo dijagonalnu prenosnu matricu objekta, tj. da dobijemo
%prenosnu matricu tako da samo jedan ulaz deluje samo na jedan izlaz
%Wdek = Wo^-1 * Wo_dijagonalno
%Wo_dijagonalno = [W11 0; 0 W22]

%odredjivanje koeficijenata brojioca i imenioca svake prenosne funkcije
[num11 den11] = tfdata(W11, 'v');
[num12 den12] = tfdata(W12, 'v');
[num21 den21] = tfdata(W21, 'v');
[num22 den22] = tfdata(W22, 'v');
%prebacivanje u simbolicki oblik koji sadrzi promenljivu "s"
syms s
num11 = poly2sym(num11, s);
den11 = poly2sym(den11, s);
num12 = poly2sym(num12, s);
den12 = poly2sym(den12, s);
num21 = poly2sym(num21, s);
den21 = poly2sym(den21, s);
num22 = poly2sym(num22, s);
den22 = poly2sym(den22, s);
%formiranje prenosnih funkcija sa simbolickom promenljivom
W11_s = num11/den11;
W12_s = num12/den12;
W21_s = num21/den21;
W22_s = num22/den22;
W_ou_s = [W11_s W12_s; W21_s W22_s];
%sada se "s" menja sa 0
W_ou_0 = double(subs(W_ou_s, s, 0));
%dijagonalna prenosna matrica objekta
W_ou_d = [W11_s 0; 0 W22_s];
W_ou_d_0 = double(subs(W_ou_d, s, 0));
%korekcioni faktor
%K_dek = W_dek(0)
K_dek = inv(W_ou_0)*W_ou_d_0;
a1 = K_dek(1, 1);
a2 = K_dek(1, 2);
a3 = K_dek(2, 1);
a4 = K_dek(2, 2);

%%
%Diskretizacija prenosne funkcije regulatora
clc

%pojacanja regulatora dobijena iz PID bloka opcijom "tune"
Kp1 = 6.58388823392313; 
Ki1 = 0.318107260315565;
disp('Vremenska konstanta integralnog clana regulatora za ventil 1:')
Ti1 = Kp1/Ki1
Wreg1 = tf([Kp1 Ki1], [1 0]);
Kp2 = -9.63504962105945; 
Ki2 = -0.943120785670637;
disp('Vremenska konstanta integralnog clana regulatora za ventil 1:')
Ti2 = Kp2/Ki2
Wreg2 = tf([Kp2 Ki2], [1 0]);
%period odabiranja
dT = 0.1;
%diskretne prenosne funkcije
Wreg1_dis = c2d(Wreg1, dT, 'tustin');
[num_d1 den_d1] = tfdata(Wreg1_dis, 'v');
Wreg2_dis = c2d(Wreg2, dT, 'tustin');
[num_d2 den_d2] = tfdata(Wreg2_dis, 'v');
%%
%Obrada slika
clc

%Poredjenje nelinearnog i linearnog neupravljanog objekta
%nivo u rezervoaru 1
load('R1_nivo_neupravljani.mat')
load('R2_nivo_neupravljani.mat')
figure(1)
plot(H1_neupravljani(1, :), H1_neupravljani(2, :), '-', H1_neupravljani(1, :), H1_neupravljani(3, :), '-.', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('H1(t) [m], h1(t) [m]', 'fontsize', 15)
legend('nelinearni', 'linearni', 'location', 'bestoutside', 'fontsize', 15)
%nivu u rezervoaru 2
figure(2)
plot(H2_neupravljani(1, :), H2_neupravljani(2, :), '-', H2_neupravljani(1, :), H2_neupravljani(3, :), '-.', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('H2(t) [m], h2(t) [m]', 'fontsize', 15)
legend('nelinearni', 'linearni', 'location', 'bestoutside', 'fontsize', 15)

%Poredjenje nelinearnog i linearnog upravljanog objekta
%nivo u rezervoaru 1
load('R1_nivo_PI_kontinualni.mat')
load('R2_nivo_PI_kontinualni.mat')
load('Uvh_PI_kontinualni.mat')
load('Uv2_PI_kontinualni.mat')
figure(3)
plot(H1_PI_kont(1, :), H1_PI_kont(2, :), '-', H1_PI_kont(1, :), H1_PI_kont(3, :), '-.', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('H1(t) [m], h1(t) [m]', 'fontsize', 15)
legend('nelinearni', 'linearni', 'location', 'bestoutside', 'fontsize', 15)
%nivo u rezervoaru 2
figure(4)
plot(H2_PI_kont(1, :), H2_PI_kont(2, :), '-', H2_PI_kont(1, :), H2_PI_kont(3, :), '-.', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('H2(t) [m], h2(t) [m]', 'fontsize', 15)
legend('nelinearni', 'linearni', 'location', 'bestoutside', 'fontsize', 15)
%upravljanje Uvh
figure(5)
plot(Uvh_PI_kont(1, :), Uvh_PI_kont(2, :), '-', Uvh_PI_kont(1, :), Uvh_PI_kont(3, :), '-.', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('Uvh(t) [V]', 'fontsize', 15)
legend('nelinearni', 'linearni', 'location', 'bestoutside', 'fontsize', 15)
%upravljanje Uv2
figure(6)
plot(Uv2_PI_kont(1, :), Uv2_PI_kont(2, :), '-', Uv2_PI_kont(1, :), Uv2_PI_kont(3, :), '-.', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('Uv2(t) [V]', 'fontsize', 15)
legend('nelinearni', 'linearni', 'location', 'bestoutside', 'fontsize', 15)

%SAR - kontinualni, diskretni_tf i diskretni_fcn
load('SAR_H1.mat')
load('SAR_H2.mat')
load('SAR_Uvh.mat')
load('SAR_Uv2.mat')
%nivo u rezervoaru 1
figure(7)
plot(H1_r(1, :), H1_r(2, :), '-', H1_r(1, :), H1_r(3, :), '-.', H1_r(1, :), H1_r(4, :), '--', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('H1(t) [m]', 'fontsize', 15)
%nivo u rezervoaru 2
figure(8)
plot(H2_r(1, :), H2_r(2, :), '-', H2_r(1, :), H2_r(3, :), '-.', H2_r(1, :), H2_r(4, :), '--', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('H2(t) [m]', 'fontsize', 15)
%upravljanje ventilom cevovoda hladne vode
figure(9)
plot(Uvh_r(1, :), Uvh_r(2, :), '-', Uvh_r(1, :), Uvh_r(3, :), '-.', Uvh_r(1, :), Uvh_r(4, :), '--', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('Uvh(t) [V]', 'fontsize', 15)
%upravljanje ventilom na izlazu rezervoara 2
figure(10)
plot(Uv2_r(1, :), Uv2_r(2, :), '-', Uv2_r(1, :), Uv2_r(3, :), '-.', Uv2_r(1, :), Uv2_r(4, :), '--', 'linewidth', 3)
grid
xlabel('t [s]', 'fontsize', 15)
ylabel('Uv2(t) [V]', 'fontsize', 15)

