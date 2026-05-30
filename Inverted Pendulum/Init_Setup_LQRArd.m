clear all
clc
% =========================================================================
% PARÁMETROS FÍSICOS NOMINALES
% =========================================================================
m = 0.1;
M_c = 0.135;
l = 0.2;
I = 0.0007176;
M = 0.136;
g = 9.81;
c = 0.63;
b = 0.00007892;
kb = 0.031;
kt = 0.031;
Rm = 12.5;
r = 0.006;
% =========================================================================
% CONFIGURACIONES INICIALES
% =========================================================================
% --- CONDICIONES INICIALES REALES (MUNDO FÍSICO) ---
x0_real     = 0.2;              % Posición inicial del carro (metros)
theta0_real = 160 * (pi/180);   % Ángulo inicial real (radianes)
x_dot0      = 0;              % Velocidad inicial carro
theta_dot0  = 0;              % Velocidad inicial péndulo

% --- CONDICIONES INICIALES PARA EL KALMAN (VARIABLES DE DESVIACIÓN) ---
% El Kalman necesita saber la desviación respecto a pi
X0_kalman = [x0_real; theta0_real - pi; x_dot0; theta_dot0];

% =========================================================================
% MATRICES DE ESPACIO DE ESTADOS NOMINALES
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

% C e D para LQR (estado completo)
C_lqr = eye(4);        % conocemos los cuatro estados
D_lqr = zeros(4,1);

% C e D para LQG (solo medimos x y theta)
% Solo medimos los estados de la posición y del ángulo. También sería
% correcto si se utilizara la 'y' que contiene la x, theta y las
% velocidades; no obstante, se debería de utilizar la matriz Id_4 (eye(4))
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
sys_tf  = tf(ss(A, B, C_lqg, D_lqg));
G_theta = sys_tf(2, 1);
G_x     = sys_tf(1, 1);

opts_x = pidtuneOptions('CrossoverFrequency', 3, 'PhaseMargin', 60);
[PID_x, ~] = pidtune(G_x, 'PID', opts_x);
Kp_x = PID_x.Kp; Ki_x = PID_x.Ki; Kd_x = PID_x.Kd;
fprintf('\n--- Ganancias PID posición x ---\n');
fprintf('Kp_x = %f\n', Kp_x);
fprintf('Ki_x = %f\n', Ki_x);
fprintf('Kd_x = %f\n', Kd_x);

[PID_theta, ~] = pidtune(G_theta, 'PID');
Kp = PID_theta.Kp; Ki = PID_theta.Ki; Kd = PID_theta.Kd;
fprintf('\n--- Ganancias PID ángulo theta ---\n');
fprintf('Kp = %f\n', Kp);
fprintf('Ki = %f\n', Ki);
fprintf('Kd = %f\n', Kd);

% =========================================================================
% 3. CONTROLADOR LQR (estado completo)
% =========================================================================
Q = diag([500, 1000, 0, 0]);
R  = 0.008;
K = lqr(A, B, Q, R);
fprintf('\n--- Ganancia LQR ---\n');
disp(K)

% Asignación de polos (para comparar con LQR)
p1 = [1i*2.8; -1i*2.8; 1i*1.5; -1i*1.5];
p2 = [-8+1i*2; -8-1i*2; -7+1i*2; -7-1i*2];
p3 = [-8; -10; -4.5; -5.8];
p4 = [-20; -15.5; -45.5; -4.8];
k_pole = place(A, B, p3);
fprintf('--- Ganancia Place (p3 estable) ---\n');
disp(k_pole)
% =========================================================================
% 4. FILTRO DE KALMAN para LQG
% =========================================================================
% =========================================================================
% 4. FILTRO DE KALMAN para LQG
% =========================================================================
G_noise = [0 0;
           0 0;
           1 0;
           0 1];

% Parámetros de ruido
Ts = 0.01;

% Varianzas (2% del rango operacional)
sigma2_x     = (0.02 * 0.5)^2;   % = 1e-4   m²
sigma2_theta = (0.02 * 0.87)^2;  % = 2.89e-4 rad²
sigma2_F     = (0.01 * 12)^2;    % 1% del rango

% Covarianza ruido de medida (2x2)
R_w = diag([sigma2_x, sigma2_theta]);

% Covarianza ruido de proceso (2x2) via densidad espectral
Q_v = diag([sigma2_F / Ts, sigma2_F / Ts]);  % 2x2

% --- Sistema para kalman() con G_noise explícito ---
System_Noise = ss(A, [B G_noise], C_lqg, [D_lqg zeros(2,2)]);
[~, L, ~] = kalman(System_Noise, Q_v, R_w);
fprintf('\n--- Ganancia Kalman L (4x2) ---\n');
disp(L)

% --- sys_kalman para bloque nativo Simulink ---
sys_kalman = ss(A, B, C_lqg, D_lqg);

% =========================================================================
% 5. PREALIMENTACIÓN Nbar
% =========================================================================
Nbar_pos = -1 / ([1 0 0 0] * ((A - B*K) \ B));
Nbar_ang = -1 / ([0 1 0 0] * ((A - B*K) \ B));
fprintf('\n--- Prealimentación Nbar ---\n');
fprintf('Nbar posición: %f\n', Nbar_pos);
fprintf('Nbar ángulo:   %f\n', Nbar_ang);

% =========================================================================
% 6. VALORES PROPIOS
% =========================================================================
fprintf('\n--- Valores propios del sistema (lazo abierto) ---\n');
disp(eig(A))
fprintf('--- Valores propios con LQR (lazo cerrado) ---\n');
disp(eig(A - B*K))
fprintf('--- Valores propios con Place p3 (lazo cerrado) ---\n');
disp(eig(A - B*k_pole))
fprintf('--- Valores propios del observador Kalman ---\n');
disp(eig(A - L*C_lqg))

% =========================================================================
% EXPORTAR AL WORKSPACE DE SIMULINK
% =========================================================================
assignin('base', 'K',           K);
assignin('base', 'L',            L);
assignin('base', 'A',            A);
assignin('base', 'B',            B);
assignin('base', 'sys_kalman',   sys_kalman);
%assignin('base', 'Q_kalman_new', Q_kalman_new);
assignin('base', 'R_w',          R_w);

% =========================================================================
% MENÚ PRINCIPAL — SELECCIÓN DE CONTROLADOR
% =========================================================================
fprintf('\n=========================================\n');
fprintf('   PÉNDULO INVERTIDO — SELECCIÓN CONTROL \n');
fprintf('=========================================\n');
fprintf('  [1] LQR (estado completo)              \n');
fprintf('  [2] LQG (Kalman + LQR)                 \n');
fprintf('  [3] Salir                               \n');
fprintf('=========================================\n\n');

ctrl = input('Selecciona el controlador (1, 2 o 3): ');

if ctrl == 1
    modelo_caso1 = 'IP_LQR_Design_Caso1';
    modelo_caso2 = 'IP_LQR_Design_Caso2';
    es_lqg = false;
    fprintf('✅ Controlador: LQR\n\n');
elseif ctrl == 2
    modelo_caso1 = 'IP_LQG_Design_Caso1';
    modelo_caso2 = 'IP_LQG_Design_Caso2';
    es_lqg = true;
    fprintf('✅ Controlador: LQG\n\n');
elseif ctrl == 3
    fprintf('👋 Saliendo. Parámetros calculados en el Workspace.\n');
    return;
else
    error('❌ Opción no válida. Escribe 1, 2 o 3.');
end

assignin('base', 'es_lqg',      es_lqg);
assignin('base', 'modelo_caso1', modelo_caso1);
assignin('base', 'modelo_caso2', modelo_caso2);

% =========================================================================
% SUBMENÚ — SELECCIÓN DE CASO
% =========================================================================
fprintf('=========================================\n');
fprintf('   SELECCIÓN DE CASO DE ESTUDIO          \n');
fprintf('=========================================\n');
fprintf('  [1] Caso 1 — Perturbaciones externas  \n');
fprintf('  [2] Caso 2 — Robustez paramétrica     \n');
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
    es_lqg   = evalin('base', 'es_lqg');
    str_ctrl = 'LQR'; if es_lqg, str_ctrl = 'LQG'; end
    modelo   = evalin('base', 'modelo_caso1');

    fprintf('\n=========================================\n');
    fprintf('   CONTROL %s — CASO 1                 \n', str_ctrl);
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

    % Generación de la señal gaussiana
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
    set_param(modelo, 'StopTime', num2str(t_stop));

    fprintf('\n✅ F_perturbacion lista | Stop time: %ds | Modelo: %s\n', t_stop, modelo);
    fprintf('   ▶ Ahora simula en Simulink\n');

    % Visualización de la señal generada
    figure('Name', sprintf('Caso 1 (%s) — Señal de perturbación', str_ctrl), 'Color', 'w');
    plot(t_sim, F_push, 'b', 'LineWidth', 2); hold on;
    yline(0, 'k--');
    for i = 1:length(F_valores)
        tp = t_inicio + (i-1)*t_entre_empujes;
        xline(tp, 'r--');
        text(tp+0.1, F_valores(i)*0.85, sprintf('%.0fN', F_valores(i)), ...
            'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
    end
    xlabel('Tiempo (s)'); ylabel('Fuerza (N)');
    title(sprintf('Caso 1 (%s): Señal de perturbación generada', str_ctrl));
    grid on;
end

% =========================================================================
% CASO 2 — ROBUSTEZ PARAMÉTRICA (ANÁLISIS MATRICIAL REAL)
% =========================================================================
function run_caso2()
    es_lqg   = evalin('base', 'es_lqg');
    str_ctrl = 'LQR'; if es_lqg, str_ctrl = 'LQG'; end
    modelo   = evalin('base', 'modelo_caso2');
    K       = evalin('base', 'K');
    L_nom    = evalin('base', 'L');
    A_nom    = evalin('base', 'A');
    B_nom    = evalin('base', 'B');

    fprintf('\n=========================================\n');
    fprintf('   CONTROL %s — CASO 2                 \n', str_ctrl);
    fprintf('   Análisis de Robustez Paramétrica Real \n');
    fprintf('=========================================\n\n');

    % Parámetros fijos nominales
    r  = 0.006;  I  = 0.0007176;  g  = 9.81;
    Rm = 12.5;   kb = 0.031;      kt = 0.031;  c = 0.63;
    l_nom = 0.2; m_nom = 0.1; M_nom = 0.136; b_nom = 0.00007892;
    C_nom = [1 0 0 0; 0 1 0 0];

    % Selección del parámetro a variar
    fprintf('¿Qué parámetro del sistema quieres alterar?\n');
    fprintf('  [1] Longitud del péndulo (l)  [Nominal: %.2fm]\n',   l_nom);
    fprintf('  [2] Masa del péndulo (m)      [Nominal: %.2fkg]\n',  m_nom);
    fprintf('  [3] Masa del carro (M)        [Nominal: %.3fkg]\n',  M_nom);
    param_opc = input('Selecciona una opción (1-3): ');

    switch param_opc
        case 1, str_p = 'Longitud (l)'; valores_test = unique([linspace(0.05, 0.60, 10), l_nom]); v_nom = l_nom;
        case 2, str_p = 'Masa pénd. (m)'; valores_test = linspace(0.02, 0.50, 10); v_nom = m_nom;
        case 3, str_p = 'Masa carro (M)'; valores_test = linspace(0.05, 0.60, 10); v_nom = M_nom;
        otherwise, error('❌ Opción inválida.');
    end

    % Barrido de estabilidad
    fprintf('\n📊 EXPERIMENTO: Evaluando estabilidad real de %s...\n', str_p);
    fprintf('-------------------------------------------------------------\n');
    fprintf('%-15s | %-20s | %-12s\n', str_p, 'Max Real(eigs)', 'Estado');
    fprintf('-------------------------------------------------------------\n');

    for i = 1:length(valores_test)
        l = l_nom; m = m_nom; M = M_nom; b = b_nom;
        switch param_opc
            case 1, l = valores_test(i);
            case 2, m = valores_test(i);
            case 3, M = valores_test(i);
        end

        % Recalcular planta real modificada
        AA_   = I*(M+m) + M*m*(l^2);
        aa_   = (((m*l)^2)*g)/AA_;
        bb_v  = ((I+m*(l^2))/AA_)*(c+(kb*kt)/(Rm*(r^2)));
        cc_v  = (b*m*l)/AA_;
        dd_   = (m*g*l*(M+m))/AA_;
        ee_   = ((m*l)/AA_)*(c+(kb*kt)/(Rm*(r^2)));
        ff_   = ((M+m)*b)/AA_;
        mm_v  = ((I+m*(l^2))*kt)/(AA_*Rm*r);
        nn_v  = (m*l*kt)/(AA_*Rm*r);

        A_real = [0 0 1 0; 0 0 0 1; 0 aa_ -bb_v -cc_v; 0 dd_ -ee_ -ff_];
        B_real = [0; 0; mm_v; nn_v];
        C_real = [1 0 0 0; 0 1 0 0];

        if es_lqg
            % Sistema aumentado 8x8: [planta real | estimador nominal fijo]
            F1 = [A_real,        -B_real * K              ];
            F2 = [L_nom * C_real, A_nom - B_nom*K - L_nom*C_nom];
            eigs_total = eig([F1; F2]);
        else
            eigs_total = eig(A_real - B_real * K);
        end

        max_real = max(real(eigs_total));
        status   = '✅ ESTABLE'; if max_real >= 0, status = '❌ INESTABLE'; end

        marker = '';
        if abs(valores_test(i) - v_nom) < 1e-9
            marker = ' ← Nominal';
        end
        fprintf('%-15.6f | %-20.4f | %s%s\n', valores_test(i), max_real, status, marker);
    end
    fprintf('-------------------------------------------------------------\n');

    % Valor definitivo para simular
    fprintf('\nIntroduce el valor exacto que deseas simular.\n');
    v_elegido = input(sprintf('Valor para %s: ', str_p));
    if v_elegido <= 0, error('❌ Debe ser un valor positivo.'); end

    l = l_nom; m = m_nom; M = M_nom; b = b_nom;
    switch param_opc
        case 1, l = v_elegido;
        case 2, m = v_elegido;
        case 3, M = v_elegido;
    end

    % Verificación final de estabilidad con valor elegido
    AA_  = I*(M+m) + M*m*(l^2);
    aa_  = (((m*l)^2)*g)/AA_;
    bb_v = ((I+m*(l^2))/AA_)*(c+(kb*kt)/(Rm*(r^2)));
    cc_v = (b*m*l)/AA_;
    dd_  = (m*g*l*(M+m))/AA_;
    ee_  = ((m*l)/AA_)*(c+(kb*kt)/(Rm*(r^2)));
    ff_  = ((M+m)*b)/AA_;
    mm_v = ((I+m*(l^2))*kt)/(AA_*Rm*r);
    nn_v = (m*l*kt)/(AA_*Rm*r);

    A_el = [0 0 1 0; 0 0 0 1; 0 aa_ -bb_v -cc_v; 0 dd_ -ee_ -ff_];
    B_el = [0; 0; mm_v; nn_v];
    C_el = [1 0 0 0; 0 1 0 0];

    if es_lqg
        F1 = [A_el,        -B_el * K                ];
        F2 = [L_nom * C_el, A_nom - B_nom*K - L_nom*C_nom];
        eigs_fin = eig([F1; F2]);
    else
        eigs_fin = eig(A_el - B_el * K);
    end

    fprintf('\n=== Resultado para %s = %.4f ===\n', str_p, v_elegido);
    if all(real(eigs_fin) < 0)
        fprintf('✅ Sistema ESTABLE con el valor elegido\n');
    else
        fprintf('❌ Sistema INESTABLE con el valor elegido\n');
    end

    % Fuerza externa (cero por defecto en Caso 2)
    t_stop_F = 5;
    t_sim_v  = 0:0.001:t_stop_F;
    assignin('base', 'F_perturbacion', [t_sim_v', zeros(size(t_sim_v))']);
    set_param(modelo, 'StopTime', num2str(t_stop_F));

    % Exportar parámetros modificados al Workspace
    assignin('base', 'l', l);
    assignin('base', 'm', m);
    assignin('base', 'M', M);
   
    fprintf('\n✅ Parámetros %s cargados con éxito.\n', str_ctrl);
    fprintf('   ▶ Ejecuta la simulación en Simulink: %s\n', modelo);
end