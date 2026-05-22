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
nn = (m*l*kt)/(AA*Rm*r);

A = [0  0   1    0;
     0  0   0    1;
     0  aa  -bb  -cc;
     0  dd  -ee  -ff];
B = [0; 0; mm; nn];

% C e D para LQR (estado completo, como en el PDF de Jitendra)
C_lqr = eye(4);
D_lqr = zeros(4,1);

% C e D para LQG (solo medimos x y theta)
C_lqg = [1 0 0 0;
         0 1 0 0];
D_lqg = zeros(2,1);

% =========================================================================
% 1. CONTROLABILIDAD Y OBSERVABILIDAD
% =========================================================================
fprintf('--- Análisis del Sistema ---\n');
fprintf('Rango Controlabilidad: %d (deseado: 4)\n', rank(ctrb(A,B)));
fprintf('Rango Observabilidad: %d (deseado: 4)\n', rank(obsv(A,C_lqr)));

fprintf('Rango Observabilidad (C 2x4): %d (deseado: 4)\n', rank(obsv(A,C_lqg)));

% =========================================================================
% 2. CONTROLADOR PID
% =========================================================================
% Función de transferencia Voltaje -> theta usando C_lqg
sys_tf  = tf(ss(A, B, C_lqg, D_lqg));
G_theta = sys_tf(2, 1);  % Canal Voltaje -> theta
function y = fcn(u)
    y = mod(u * 180/pi, 360);
end

% PID externo - controla posición x
G_x = sys_tf(1, 1);
opts_x = pidtuneOptions('CrossoverFrequency', 3, 'PhaseMargin', 60);  % añadir esta línea
[PID_x, ~] = pidtune(G_x, 'PID', opts_x);
Kp_x = PID_x.Kp;
Ki_x = PID_x.Ki;
Kd_x = PID_x.Kd;

fprintf('\n--- Ganancias PID posición x ---\n');
fprintf('Kp_x = %f\n', Kp_x);
fprintf('Ki_x = %f\n', Ki_x);
fprintf('Kd_x = %f\n', Kd_x);

[PID_theta, ~] = pidtune(G_theta, 'PID');
Kp = PID_theta.Kp;
Ki = PID_theta.Ki;
Kd = PID_theta.Kd;
fprintf('\n--- Ganancias PID ---\n');
fprintf('Kp = %f\n', Kp);
fprintf('Ki = %f\n', Ki);
fprintf('Kd = %f\n', Kd);


% =========================================================================
% 3. CONTROLADOR LQR (estado completo)
% =========================================================================
Q = diag([500, 1000, 0, 0]);
R = 0.008;
KK = lqr(A, B, Q, R);
fprintf('\n--- Ganancia LQR ---\n');
disp(KK)

% Asignación de polos (para comparar con LQR)
p1 = [1i*2.8; -1i*2.8; 1i*1.5; -1i*1.5];    % Oscilatorio
p2 = [-8+1i*2; -8-1i*2; -7+1i*2; -7-1i*2];  % Subamortiguado
p3 = [-8; -10; -4.5; -5.8];                   % Estable
p4 = [-20; -15.5; -45.5; -4.8];               % Agresivo
k_pole = place(A, B, p3);
fprintf('--- Ganancia Place (p3 estable) ---\n');
disp(k_pole)

% =========================================================================
% 4. FILTRO DE KALMAN para LQG
% =========================================================================
% Ruido de proceso entra por las aceleraciones (estados 3 y 4)
G_noise = [0 0;
           0 0;
           1 0;
           0 1];

% original: Q_v = diag([0.05, 0.05])
Q_v = diag([500, 500]);   % Covarianza ruido de proceso
R_w = diag([0.001, 0.002]); % Covarianza ruido de medida

% Definir sistema para el bloque nativo
sys_kalman = ss(A, B, C_lqg, D_lqg);

% Q de 4x4 que se usa para el bloque de Kalman nativo del Simulink
Q_kalman = G_noise * Q_v * G_noise';
% Q_kalman =
% [0, 0, 0, 0]
% [0, 0, 0, 0]
% [0, 0, Q_v, 0]
% [0, 0, 0, Q_v]

sys_kalman = ss(A, B, C_lqg, D_lqg);

System_Noise = ss(A, [B G_noise], C_lqg, [D_lqg zeros(2,2)]);
[~, L, ~] = kalman(System_Noise, Q_v, R_w);
fprintf('\n--- Ganancia Kalman L ---\n');
disp(L)

% =========================================================================
% 5. PREALIMENTACIÓN Nbar
% =========================================================================
Nbar_pos = -1 / ([1 0 0 0] * ((A - B*KK) \ B));
Nbar_ang = -1 / ([0 1 0 0] * ((A - B*KK) \ B));
fprintf('\n--- Prealimentación Nbar ---\n');
fprintf('Nbar posición: %f\n', Nbar_pos);
fprintf('Nbar ángulo:   %f\n', Nbar_ang);

% =========================================================================
% 6. Valores propios
% =========================================================================
fprintf('\n--- Valores propios del sistema (lazo abierto) ---\n');
disp(eig(A))

fprintf('--- Valores propios con LQR (lazo cerrado) ---\n');
disp(eig(A - B*KK))

fprintf('--- Valores propios con Place p3 (lazo cerrado) ---\n');
disp(eig(A - B*k_pole))

fprintf('--- Valores propios del observador Kalman ---\n');
disp(eig(A - L*C_lqg))







% =========================================================================
% MENÚ PRINCIPAL
% =========================================================================
fprintf('=========================================\n');
fprintf('   CONTROL LQR — PÉNDULO INVERTIDO      \n');
fprintf('=========================================\n');
fprintf('  [1] Caso 1 — Perturbaciones externas  \n');
fprintf('  [2] Caso 2 — Robustez paramétrica (l) \n');
fprintf('=========================================\n\n');

caso = input('Selecciona el caso (1 o 2): ');

if caso == 1
    run_caso1();
elseif caso == 2
    run_caso2();
else
    error('❌ Opción no válida. Escribe 1 o 2.');
end

% =========================================================================
% CASO 1 — PERTURBACIONES EXTERNAS
% =========================================================================
function run_caso1()
    fprintf('\n=========================================\n');
    fprintf('   CONTROL LQR — CASO 1                 \n');
    fprintf('   Test de robustez ante perturbaciones  \n');
    fprintf('=========================================\n\n');

    F_valores_todos = [2, 4, 6, 8, 10, 12, 15];
    dur_empuje      = 0.2;
    t_inicio        = 3;

    fprintf('¿Qué quieres simular?\n');
    fprintf('  [1] Todos los empujes incrementales (2N → 15N)\n');
    fprintf('  [2] Un solo empuje (una única vez)\n');
    fprintf('  [3] Un solo empuje repetido cada X segundos\n');
    opcion = input('Selecciona opción (1, 2 o 3): ');

    if opcion == 1
        F_valores       = F_valores_todos;
        t_entre_empujes = 3;
        t_stop          = 25;
        fprintf('✅ Simulando todos los empujes | Stop time: %ds\n', t_stop);

    elseif opcion == 2
        F_single = input(sprintf('Introduce la fuerza en N %s: ', mat2str(F_valores_todos)));
        if ~ismember(F_single, F_valores_todos)
            warning('Fuerza no estándar, se usará igualmente.');
        end
        F_valores       = F_single;
        t_entre_empujes = 3;
        t_stop          = 10;
        fprintf('✅ Empuje único de %.0fN | Stop time: %ds\n', F_single, t_stop);

    elseif opcion == 3
        F_single        = input(sprintf('Introduce la fuerza en N %s: ', mat2str(F_valores_todos)));
        if ~ismember(F_single, F_valores_todos)
            warning('Fuerza no estándar, se usará igualmente.');
        end
        t_entre_empujes = input('Cada cuántos segundos se repite el empuje: ');
        t_stop          = input('Stop time (s): ');
        n_rep           = floor((t_stop - t_inicio) / t_entre_empujes);
        F_valores       = repmat(F_single, 1, n_rep);
        fprintf('✅ %.0fN cada %ds → %d empujes | Stop time: %ds\n', ...
            F_single, t_entre_empujes, n_rep, t_stop);
    else
        error('❌ Opción no válida.');
    end

    % Generación de la señal
    t_total = t_inicio + length(F_valores) * t_entre_empujes;
    t_sim   = 0:0.001:max(t_total, t_stop);
    F_push  = zeros(size(t_sim));

    fprintf('\n=== Secuencia de empujes ===\n');
    for i = 1:length(F_valores)
        tp    = t_inicio + (i-1)*t_entre_empujes;
        amp   = F_valores(i);
        sigma = dur_empuje/4;
        idx   = abs(t_sim - tp) < 3*sigma;
        F_push(idx) = amp * exp(-((t_sim(idx)-tp).^2)/(2*sigma^2));
        fprintf('Empuje %2d | t=%5.1fs | F=+%.1fN\n', i, tp, amp);
    end

    assignin('base', 'F_perturbacion', [t_sim', F_push']);
    set_param('IP_LQR_Design_Caso1', 'StopTime', num2str(t_stop));
    fprintf('\n✅ F_perturbacion lista | Stop time: %ds\n', t_stop);
    fprintf('   ▶ Ahora simula en Simulink\n');

    % Visualización señal
    figure('Name','Caso 1 — Señal de perturbación','Color','w');
    plot(t_sim, F_push, 'b', 'LineWidth', 2); hold on;
    yline(0,'k--');
    for i = 1:length(F_valores)
        tp = t_inicio + (i-1)*t_entre_empujes;
        xline(tp,'r--');
        text(tp+0.1, F_valores(i)*0.85, sprintf('%.0fN', F_valores(i)), ...
            'Color','r','FontSize',10,'FontWeight','bold');
    end
    xlabel('Tiempo (s)'); ylabel('Fuerza (N)');
    title('Caso 1: Señal de perturbación generada');
    grid on;
end

% =========================================================================
% CASO 2 — ROBUSTEZ PARAMÉTRICA
% =========================================================================
function run_caso2()
    fprintf('\n=========================================\n');
    fprintf('   CONTROL LQR — CASO 2                 \n');
    fprintf('   Robustez paramétrica                  \n');
    fprintf('   Variación de longitud del péndulo (l) \n');
    fprintf('=========================================\n\n');

    % Parámetros fijos
    r  = 0.006;  I  = 0.0007176;  g  = 9.81;
    b  = 0.00007892;  Rm = 12.5;
    kb = 0.031;  kt = 0.031;  c  = 0.63;
    m  = 0.1;   M  = 0.136;
    l_nominal = 0.2;

    % LQR nominal (no cambia)
    AA   = I*(M+m) + M*m*(l_nominal^2);
    aa   = (((m*l_nominal)^2)*g)/AA;
    bb   = ((I+m*(l_nominal^2))/AA)*(c+(kb*kt)/(Rm*(r^2)));
    cc_v = (b*m*l_nominal)/AA;
    dd   = (m*g*l_nominal*(M+m))/AA;
    ee   = ((m*l_nominal)/AA)*(c+(kb*kt)/(Rm*(r^2)));
    ff   = ((M+m)*b)/AA;
    mm_v = ((I+m*(l_nominal^2))*kt)/(AA*Rm*r);
    nn   = (m*l_nominal*kt)/(AA*Rm*r);

    A_nom = [0  0    1     0;
             0  0    0     1;
             0  aa  -bb   -cc_v;
             0  dd  -ee   -ff];
    B_nom = [0; 0; mm_v; nn];

    Q  = diag([500, 1000, 0, 0]);
    R  = 0.008;
    KK = lqr(A_nom, B_nom, Q, R);
    fprintf('KK diseñado con l_nominal = %.2fm (fijo)\n\n', l_nominal);

    % -------------------------------------------------------------------------
    % PREGUNTA AL USUARIO
    % -------------------------------------------------------------------------
    fprintf('Longitud nominal del péndulo: %.2fm\n', l_nominal);
    fprintf('Introduce la longitud que quieres probar (en metros).\n');
    fprintf('Ejemplos: 0.10, 0.15, 0.20, 0.30, 0.45\n\n');
    l = input('Longitud del péndulo (m): ');

    if l <= 0
        error('❌ La longitud debe ser un valor positivo.');
    end

    % -------------------------------------------------------------------------
    % CALCULAR ESTABILIDAD CON LA LONGITUD ELEGIDA
    % -------------------------------------------------------------------------
    AA   = I*(M+m) + M*m*(l^2);
    aa   = (((m*l)^2)*g)/AA;
    bb   = ((I+m*(l^2))/AA)*(c+(kb*kt)/(Rm*(r^2)));
    cc_v = (b*m*l)/AA;
    dd   = (m*g*l*(M+m))/AA;
    ee   = ((m*l)/AA)*(c+(kb*kt)/(Rm*(r^2)));
    ff   = ((M+m)*b)/AA;
    mm_v = ((I+m*(l^2))*kt)/(AA*Rm*r);
    nn_v = (m*l*kt)/(AA*Rm*r);

    A_l = [0  0    1     0;
           0  0    0     1;
           0  aa  -bb   -cc_v;
           0  dd  -ee   -ff];
    B_l = [0; 0; mm_v; nn_v];

    eigs_lc = eig(A_l - B_l*KK);
    estable = all(real(eigs_lc) < 0);

    % -------------------------------------------------------------------------
    % RESUMEN EN CONSOLA
    % -------------------------------------------------------------------------
    fprintf('\n=== Resultado para l = %.2fm ===\n', l);
    fprintf('Valores propios (lazo cerrado):\n');
    for k = 1:length(eigs_lc)
        fprintf('  λ%d = %.4f %+.4fi\n', k, real(eigs_lc(k)), imag(eigs_lc(k)));
    end
    fprintf('\n');
    if estable
        fprintf('✅ Sistema ESTABLE con l = %.2fm\n', l);
        if l < l_nominal
            fprintf('   Péndulo más corto que el nominal (%.2fm)\n', l_nominal);
        elseif l > l_nominal
            fprintf('   Péndulo más largo que el nominal (%.2fm)\n', l_nominal);
        else
            fprintf('   Longitud nominal — comportamiento esperado\n');
        end
    else
        fprintf('❌ Sistema INESTABLE con l = %.2fm\n', l);
        fprintf('   El LQR diseñado para l=%.2fm no puede estabilizar este péndulo.\n', l_nominal);
    end

    % -------------------------------------------------------------------------
    % PREGUNTA FUERZA EXTERNA
    % -------------------------------------------------------------------------
    fprintf('\n¿Quieres aplicar una fuerza externa?\n');
    fprintf('  [1] Sí\n');
    fprintf('  [2] No\n');
    opcion_F = input('Selecciona opción (1 o 2): ');
    
    if opcion_F == 1
        F_val     = input('Fuerza a aplicar (N): ');
        t_inicio  = input('En qué segundo aplicar el empuje: ');
        dur_empuje = 0.2;
        t_stop_F  = input('Stop time (s): ');
    
        t_sim_v = 0:0.001:t_stop_F;
        sigma   = dur_empuje/4;
        F_push  = F_val * exp(-((t_sim_v - t_inicio).^2)/(2*sigma^2));
        assignin('base', 'F_perturbacion', [t_sim_v', F_push']);
    
        fprintf('\n✅ Empuje de %.0fN en t=%.0fs generado\n', F_val, t_inicio);
    
        if estable
            fprintf('   Sistema teóricamente estable → debería recuperar\n');
        else
            fprintf('   ⚠️  Sistema teóricamente inestable → puede divergir\n');
        end
    
        set_param('IP_LQR_Design_Caso2', 'StopTime', num2str(t_stop_F));
    
    else
        % Sin fuerza — señal cero
        t_stop_F = 15;
        t_sim_v  = 0:0.001:t_stop_F;
        assignin('base', 'F_perturbacion', [t_sim_v', zeros(size(t_sim_v))']);
        set_param('IP_LQR_Design_Caso2', 'StopTime', num2str(t_stop_F));
        fprintf('\n✅ Sin fuerza externa\n');
    end


    % -------------------------------------------------------------------------
    % ACTUALIZAR WORKSPACE Y SIMULINK
    % -------------------------------------------------------------------------
    assignin('base', 'l',  l);
    assignin('base', 'KK', KK);
    set_param('IP_LQR_Design_Caso2', 'StopTime', '15');

    fprintf('\n✅ l = %.2fm cargado en Simulink | Stop time: 15s\n', l);
    fprintf('   ▶ Ahora simula en Simulink para ver la respuesta\n');
end
