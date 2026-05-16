clear all
clc
r = 0.006;
M_c = 0.135;
I = 0.0007176;
l = 0.2;   % podemos cambiar para que el péndulo es 0.6m
% Más complicado de controlar o no?
% En cuanto la 'x' no cambia (lo de moverme me da igual). La primera
% posición.

% Pero ahora para controlar el péndulo tengo que aplicar más fuerza (para
% controlar el ángulo). La segunda posición.

% Gasto más si el carro es gordo o la masa de arriba es gordo (acción de
% control más grande)

g = 9.81;
b= 0.00007892;
L = 0.046;
Rm = 12.5;
kb = 0.031;
kt = 0.031;
c = 0.63;
m = 0.1;    % la masa que está arriba (m=0.4)
% La bola que está arriba es 4 veces más pesada


M = 0.136;
Er = 2*m*g*l;
n= 3;
AA = I*(M+m) + M*m*(l^2);
aa = (((m*l)^2)*g)/AA;
bb = ((I +m*(l^2))/AA)*(c + (kb*kt)/(Rm*(r^2)));
cc  = (b*m*l)/AA;
dd  = (m*g*l*(M+m))/AA;
ee  = ((m*l)/AA)*(c + (kb*kt)/(Rm*(r^2)));
ff  = ((M+m)*b)/AA;
mm = ((I +m*(l^2))*kt)/(AA*Rm*r);
nn = (m*l*kt)/(AA*Rm*r);
A  =  [0 0 1 0; 0 0 0 1; 0 aa -bb -cc; 0 dd -ee -ff];
B = [0;0; mm; nn]; 
Q = diag([1200 1500 0 0]);
R  = 0.035;
KK = lqr(A,B,Q,R) 
p1 = [i*2.8; -i*2.8; i*1.5; -i*1.5]; % oscillatory
p2 =[-8+i*2; -8-i*2; -7+i*2; -7-i*2 ]; % underdamped
p3 =[-8; -10; -4.5; -5.8];  % stable. 
p4 =[-20; -15.5; -45.5; -4.8];  % Fast or Aggressive.
k = place(A,B,p3);



% Peras-manzanas (posición-ángulo)
% 100 metros, me puedo equivocar en 100, ahora multiplico 100 por 1000
% la grua va colgando para llevar peso -> se mueve (muy poco) 20 grados,
% como trabajo en SI estoy trabajando con radianes -> 0'8. 
%

% Si fuera 0'8 y 100, en vez de 100 sería 1000 más 2 ceros (desde un punto
% de vista de 0/1). Tenemos que trabajar en un escalado de 0/1 o -1/1.
% TENEMOS QUE ESCALAR (TENEMOS QUE ESCALAR BIEN LOS VALORES). DETERMINAR
% BIEN LA ZONA DE TRABAJO.

% R = 0.035 (valor petit). 1. vull que sigui petit o 2. perquè la zona en
% la qual treballa l'acció de control no és 0-1. 
%
% Corriente continua, 24V coche (estamos trabajando entre 0-24). Este 0,035
% para escalarlo lo tendríamos que multiplicar por 24 para compararlo con
% lo del 1000 100??  TENEMOS QUE DEFINIRLO PORQUE ES IMPORTANTE.

% Si no tenemos en cuenta esto es porque NO HEMOS TENIDO EN CUENTA TODO
% ESTO.
%
