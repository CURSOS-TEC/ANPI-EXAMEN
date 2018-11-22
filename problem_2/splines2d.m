#!/usr/bin/octave-cli

## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## II Semestre 2018
## Examen Final

## PROBLEMA 2

## NOMBRE: JUan Esteban Navarro Camacho
## CARNE: 201236227
2;


## Construya algunos datos 2D para el problema
## N: Número de datos.  Cada fila de la matriz tendrá un punto.
##    La primera columna tendrá la coordenada x y la segunda columna
##    la coordenada y.
function points = createData(N)
  astep = 360/N;
  
  angles = (0:astep:360-astep)';
  anoise = rand(size(angles))*astep/4;
  rnoise = rand(size(angles))*1.5;
  
  radii  = 2+rnoise;
  angles = deg2rad(angles + anoise);
  
  points = [radii.*cos(angles) radii.*sin(angles)];
endfunction

## Calcule las segundas derivadas
## x: posiciones x de las muestras
## f: valores de la función en cada x
## retorne fpp con los valores de la segunda derivada en cada posición t
function fpp=findDerivs(t,f)
  assert(size(f)==size(t));

  N=length(t)-1; # Número de subintervalos
  
  ## Arme el sistema de ecuaciones
  M=eye(N+1,N+1);
  fpp=zeros(N+1,1);
  b  =zeros(N+1,1);
  

  ## ################## 
  ## ## Problema 2.4 ##
  ## ################## 
  for i = 2:N
   #X_i es t_i+1
   M(i,i-1)   = t(i)-t(i-1);
   M(i,i) = 2*(t(i+1)-t(i-1));
   M(i,i+1) = t(i+1)-t(i);
   
   A = f(i+1) - f(i);
   B = t(i+1) - t(i);
   
   C = f(i) - f(i-1);
   D = t(i) - t(i-1);
   
   b(i) = 6*(A/B) - 6*(C/D);
   
  endfor
  #primera fila
  M(1,N+1) = t(1)-t(N+1);
  M(1,1) =2*( t(2)-t(N+1));
  M(1,2) = t(2)-t(1);
  
  b(1) = 6*(f(2)-f(1))/(t(2)-t(1)) - 6*(f(1)-f(N+1))/(t(1)-t(N+1));
  
  
  
  #última fila
  M(N+1,N) = t(N+1)-t(N);
  M(N+1,N+1) =2*( t(1)-t(N));
  M(N+1,1) = t(1)-t(N+1);
  
  b(N+1) =  6*(f(1)-f(N+1))/(t(1)-t(N+1)) - 6*(f(N+1)-f(N))/(t(N+1)-t(N));

  ## >>> Ponga su solución aquí <<<
  
  
  ## Resuelva el sistema
  fpp = M\b; 
  
endfunction

## Interpole los valores fi(xi) usando los puntos x y sus valores f(x)
## t: valores de soporte conocidos
## f: valores de la función en los x conocidos
## ts: valores en donde debe encontrarse la función interpolada
## retorna fs: valores de la función en los xs dados
function fs=interpole(t,f,ts)
  assert(size(t)==size(f));
  ts=ts(:); ## Asegúrese de que es un vector columna
  ## Encuentre las segundas derivadas
  fpp=findDerivs(t,f);

  ## ##################
  ## ## Problema 2.5 ##
  ## ##################

  ## >>> Ponga su solución aquí <<<
  
  ## Sugerencia: Puede serle muy útil el uso de la función 'lookup'
  ##             para encontrar cuál subintervalo utilizar.
  fs=zeros(size(ts));
  for i=1:numel(ts)
    xi1 = lookup(t, ts(i));
    xi = xi1+1;
    if(xi>numel(t))
      xi=1;
    endif
    fpp_xi1 = fpp(xi1);
    fpp_xi  = fpp(xi);
    fxi_1   = f(xi1);
    fxi     = f(xi1);
    
    x = ts(i);
    x_i = t(xi);
    x_i_1 = t(xi1);
    
    c1 = fpp_xi1*(x - t(xi))^3/(6*(t(xi1)-t(xi)))+fpp_xi*(ts(i) - t(xi1))^3/(6*(t(xi)-t(xi1)));
    c2 = (f(xi1)/(t(xi1)-t(xi))-fpp_xi1*(t(xi1)-t(xi))/6)*(ts(i)-t(xi));
    c3 = (f(xi)/(t(xi)-t(xi1))-fpp_xi*(t(xi)-t(xi1))/6)*(ts(i)-t(xi1));
    fs(i) = c1+c2+c3;
    
  endfor
  

  
  
endfunction

## Depuración
figure(2,"name","Interpolación simple cerrada (depuración)");
x=[0,1,2,3];
f=[1,2,1,0.5];

hold off;
plot(x,f,'rx-;original;',"linewidth",2);

step=0.1;
xs=0:step:4-step;

fs=interpole(x,f,xs);

hold on;
plot(xs,fs,'bo-;interpolado;');
grid on;
xlabel("t");
ylabel("f(t)");

## El caso completo
N=10;
D = createData(N);
figure(1,"name","Interpolación 2D cerrada");
hold off;
plot(D(:,1),D(:,2),'rx-',"linewidth",2);

step=0.1;
t=0:step:N-step;
xs=interpole([0:N-1]',D(:,1),t);
ys=interpole([0:N-1]',D(:,2),t);

## ##################
## ## Problema 2.6 ##
## ##################

## >>> Ponga su solución aquí <<<

hold on;
plot(xs,ys,'bo-');
xlabel("x");
ylabel("y");
grid;
