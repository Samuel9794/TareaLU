clc; clear;

% ==========================================
% INGRESA AQUI LOS DATOS DE TU EJERCICIO
% (REEMPLAZA CON LOS DE LA FOTO)
% ==========================================
A = [ 10  2  -1;
      -3  -6  2;
      1  1  5 ];

b = [27; -61.5; -21.5];

n = length(b);

disp("Matriz A:");
disp(A);
disp("Vector b:");
disp(b);

% ==========================================
% METODO DE GAUSS
% ==========================================
disp("=====================================");
disp("METODO DE GAUSS");

Ag = A;
bg = b;

for i = 1:n-1
  if Ag(i,i) == 0
    error("Pivote cero detectado");
  endif

  for k = i+1:n
    factor = Ag(k,i)/Ag(i,i);

    printf("F%d = F%d - (%.3f)F%d\n", k, k, factor, i);

    Ag(k,i:n) = Ag(k,i:n) - factor*Ag(i,i:n);
    bg(k) = bg(k) - factor*bg(i);
  endfor
endfor

% Sustitución hacia atrás
x_gauss = zeros(n,1);
for i = n:-1:1
  suma = bg(i);
  for j = i+1:n
    suma = suma - Ag(i,j)*x_gauss(j);
  endfor
  x_gauss(i) = suma / Ag(i,i);
endfor

disp("Solución por Gauss:");
disp(x_gauss);

% ==========================================
% METODO LU (DOOLITTLE)
% ==========================================
disp("=====================================");
disp("METODO LU - DOOLITTLE");

L = eye(n);
U = zeros(n);

for j = 1:n
  for i = 1:j
    suma = 0;
    for k = 1:i-1
      suma += L(i,k)*U(k,j);
    endfor
    U(i,j) = A(i,j) - suma;
  endfor

  for i = j+1:n
    suma = 0;
    for k = 1:j-1
      suma += L(i,k)*U(k,j);
    endfor

    if U(j,j) == 0
      error("Pivote cero en LU");
    endif

    L(i,j) = (A(i,j) - suma)/U(j,j);
  endfor
endfor

disp("Matriz L:");
disp(L);
disp("Matriz U:");
disp(U);

% Ly = b
y = zeros(n,1);
for i = 1:n
  suma = b(i);
  for j = 1:i-1
    suma -= L(i,j)*y(j);
  endfor
  y(i) = suma;
endfor

% Ux = y
x_lu = zeros(n,1);
for i = n:-1:1
  suma = y(i);
  for j = i+1:n
    suma -= U(i,j)*x_lu(j);
  endfor
  x_lu(i) = suma/U(i,i);
endfor

disp("Solución por LU (Doolittle):");
disp(x_lu);

disp("=====================================");
disp("FIN");
