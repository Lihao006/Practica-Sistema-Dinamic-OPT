clear all
clc

% =========================================================================
% PARÁMETROS FÍSICOS
% =========================================================================
r = 0.006;
M_c = 0.135;
I = 0.0007176;
l = 0.2;       
g = 9.81;
b = 0.00007892;
Rm = 12.5;
kb = 0.031;
kt = 0.031;
c = 0.63;
m = 0.1;
M = 0.136;

% =========================================================================
% MATRICES DE ESPACIO DE ESTADOS
% =========================================================================
AA = I*(M+m) + M*m*(l^2);
aa = (((m*l)^2)*g)/AA;
bb = ((I+m*(l^2))/AA)*(c + (kb*kt)/(Rm*(r^2)));
cc = (b*m*l)/AA;
dd = (m*g*l*(M+m))/AA;
ee = ((m*l)/AA)*(c + (kb*kt)/(Rm*(r^2)));
ff = ((M+m)*b)/AA;
mm = ((I+m*(l^2))*kt)/(AA*Rm*r);
nn = (m*l*kt)/(AA*Rm*r);   % <-- corregido: era kt no t

A = [0  0   1    0;
     0  0   0    1;
     0  aa  -bb  -cc;
     0  dd  -ee  -ff];

B = [0; 0; mm; nn];

C = [1 0 0 0;    % Medimos posición x
     0 1 0 0];   % Medimos ángulo theta
D = zeros(2,1);

% =========================================================================
% 1. CONTROLABILIDAD Y OBSERVABILIDAD
% =========================================================================
fprintf('--- Análisis del Sistema ---\n');
fprintf('Rango Controlabilidad: %d (deseado: 4)\n', rank(ctrb(A,B)));
fprintf('Rango Observabilidad:  %d (deseado: 4)\n', rank(obsv(A,C)));

% =========================================================================
% 2. CONTROLADOR LQR
% =========================================================================
Q = diag([1200, 1500, 0, 0]);
R = 0.035;
KK = lqr(A, B, Q, R);
fprintf('\n--- Ganancia LQR ---\n');
disp(KK)

% Asignación de polos (para comparar en el report)
p1 = [1i*2.8; -1i*2.8; 1i*1.5; -1i*1.5];       % Oscilatorio
p2 = [-8+1i*2; -8-1i*2; -7+1i*2; -7-1i*2];      % Subamortiguado
p3 = [-8; -10; -4.5; -5.8];                       % Estable
p4 = [-20; -15.5; -45.5; -4.8];                   % Agresivo
k_pole = place(A, B, p3);

% =========================================================================
% 3. FILTRO DE KALMAN (LQG)
% =========================================================================
% G: el ruido de proceso entra por las aceleraciones (estados 3 y 4)
G = [0 0;
     0 0;
     1 0;
     0 1];

% Covarianza ruido de proceso (vibraciones mecánicas en aceleraciones)
Q_v = diag([0.05, 0.05]);

% Covarianza ruido de medida (imprecisión encoders de x y theta)
R_w = diag([0.001, 0.002]);

% Sistema aumentado con canales de ruido
System_Noise = ss(A, [B G], C, [D zeros(2,2)]);

% Cálculo de la ganancia del estimador de Kalman
[~, L, P] = kalman(System_Noise, Q_v, R_w);

fprintf('--- Ganancia Kalman L ---\n');
disp(L)

% =========================================================================
% 4. PREALIMENTACIÓN Nbar (CASOS LIBRES)
% =========================================================================
% Caso 1: referencia de posición
Nbar_pos = -1 / ([1 0 0 0] * ((A - B*KK) \ B));

% Caso 2: referencia de ángulo
Nbar_ang = -1 / ([0 1 0 0] * ((A - B*KK) \ B));

fprintf('\n--- Prealimentación Nbar ---\n');
fprintf('Nbar posición: %f\n', Nbar_pos);
fprintf('Nbar ángulo:   %f\n', Nbar_ang);