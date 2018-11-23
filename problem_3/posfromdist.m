#!/usr/bin/octave-cli

## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## II Semestre 2018
## Examen Final

## PROBLEMA 3

## NOMBRE:Juan Esteban Navarro Camacho
## CARNE:201236227

3;
## Cargue los datos de distancias
##
## Cada columna de D tiene las distancias a un sensor fijo en el
## espacio. 
D=load("-ascii","dists.dat");

## Las posiciones de los emisores están dadas por las _columnas_
## de la matriz S
E=load("-ascii","emisors.dat");

## Función calcula posición a partir de distancias a emisores.
##
## dists:     Distancias del punto a cada emisor como vector fila.
## emisorPos: Posiciones de los emisores, cada emisor en una columna.
## option:    Si solo hay 3 emisores, cuál de las soluciones retorna:
##            1: use + en la cuadrática
##            2: use - en la cuadrática
##            3: decida cuál signo usar dependiendo de la velocidad y
##               aceleración
##
## La posición i de dists contiene la distancia al emisor en la
## columna i.
##
function p=calcPosition(dists,emisorPos,option=1)

  ## Número de emisores usados
  dim = min(columns(emisorPos),length(dists));

  ## ##################
  ## ## Problema 3.1 ##
  ## ##################
  ## Construya la matriz M
  ## <<< Ponga su solución aquí
  M  = zeros(dim,4); 
  #iterar sobre la matriz de posición de emisores
  for i=1:size(M)(1)#filas
    for j=1:size(M)(2)#columnas
        if(j == 1)
            M(i,j) = 1;
        else
            M(i,j) = -2*emisorPos(j-1,i); 
        endif

    endfor
  endfor 

  M=[M];
  
  M; # impresión para depurar, #el breakpoint va aquí
  

  ## ##################
  ## ## Problema 3.2 ##
  ## ##################
  ## Construya el vector b
   b = zeros(dim,1); 
   
   for k = 1:size(b)(1)
     #se define el i-ésimo emisor
     emisor = [emisorPos(1,k) emisorPos(2,k) emisorPos(3,k)];
     b(k,1) = (dists(k)*dists(k)) - (norm(emisor)*norm(emisor));
   endfor
   b1 = b;
  
   b1;
  
  

  ## ##################
  ## ## Problema 3.3 ##
  ## ##################

  ## Calcule la matriz seudo-inversa utilizando SVD
  iM=pinv(M); ## <<< Ponga su solución aquí (se busca calcular esto con SVD)
  #Retorna un arreglo de tres valores
  [X,Y,Z] = svd(M);
  
  x = resize(X,rows(M),columns(M));
  
  y = resize(Y,columns(M));
  
  z = resize(Z,columns(M));
  
  pinv(y);
 
  
  iM = z * pinv(y) * x';
  

  ## Verifique que iM y pinv(M) son lo mismo
  if (norm(iM-pinv(M),"fro") > 1e-6)
    error("Matriz inversa calculada con SVD incorrecta");
  endif

  ## ##################
  ## ## Problema 3.4 ##
  ## ##################
  
  ## Calcule la solución particular
  hatp= iM * b1; ## <<< Ponga su solución aquí
    
  
  ## El caso de 3 dimensiones tiene dos posibles soluciones:
  if (dim==3)
    ## Con 3 emisores, calcule las dos posibles posiciones
    
    ## ##################
    ## ## Problema 3.5 ##
    ## ##################

    ## >>> Ponga su solución aquí <<<
    vec = z(:,4);
    
    v1 = (vec(2,1)^2) + (vec(3,1)^2) + (vec(4,1)^2);
    
    v2 = ((2*hatp(2,1)*vec(2,1)) + (2*hatp(3,1)*vec(3,1)) + (2*hatp(4,1)*vec(4,1))) - vec(1,1);
    
    v3 = ((hatp(2,1)^2) + (hatp(3,1)^2) + (hatp(4,1)^2)) - hatp(1,1);
    
    if(option == 1)
      q = 2*v1;
      dis = 4*v1*v3;
      l = (-v2 + sqrt((v2^2) - (dis))) / q; 
      
      p = (hatp + l )(1:3);
      
      
    else
      q = 2*v1;
      dis = 4*v1*v3;
      l = (-v2 - sqrt((v2^2) - (dis))) / q; 
      
      p = (hatp + l )(1:3);
    endif
  else 
    ## Caso general de más de tres dimensiones

    ## ##################
    ## ## Problema 3.6 ##
    ## ##################

    ## >>> Ponga su solución aquí <<<
    p = hatp(2:4);
  endif
  
endfunction

## Calcule la trayectoria de posiciones para una matriz de
## distancias que contiene en cada fila el vector de distancias
## a cada emisor.
## 
## Sea N el número de datos y S el número de emisores
## dists:     matriz de distancias de tamaño N x S
## emisorPos: matrix de tamaño 3xS
## option:    si S==3, entonces cuál de las dos soluciones devolver
##            option=1 implica usar + en cuadrática, y otra cosa usa -
function p=calcPositions(dists,emisorPos,option=1)
  n=min(columns(emisorPos),columns(dists));
  k=rows(dists);

  ## Reserve memoria para todos los puntos en la trayectoria
  p=zeros(k,3);

  if (option==3)
    
    predPos=[2;0;1]; ## Predicción inicial de posición

    predVel=[0;0;0];
    predAcc=[0;0;0];
    ## Para cada punto en la secuencia
    for i=1:k

      ## Calcule las dos posibilidades
      sol1=calcPosition(dists(i,1:n),emisorPos(:,1:n),1)';
      sol2=calcPosition(dists(i,1:n),emisorPos(:,1:n),2)';
      
      ## Verifique cuál solución está más cerca de la predicción
      if (norm(predPos-sol1) < norm(predPos-sol2))
	p(i,:) = sol1;
      else
	p(i,:) = sol2;
      endif

      ## ##################
      ## ## Problema 3.7 ##
      ## ##################
      
      ## >>> Complete su solución aquí <<<
      
      ## Actualice la predicción
      
    endfor
  else
    ## Para todos los puntos en la trayectoria
    for i=1:k
      ## Calcule la posición del punto, dadas las distancias a los emisores
      p(i,:)=calcPosition(dists(i,1:n),emisorPos(:,1:n),option)';
    endfor
  endif

endfunction

## Función calcula distancia de un punto a todos los emisores
## pos:       Posición del objeto (en las filas)
## emisorPos: Posiciones de los emisores, en las columnas.
function d=calcDistances(pos,emisorPos)
  if (columns(pos)!=3)
    error("Position must be a 3D vector");
  endif

  ## Calcule la distancia a cada emisor, dada la posición "pos" del
  ## objeto y las posiciones de los emisores.
  
  ## ##################
  ## ## Problema 3.8 ##
  ## ##################

  ## >>> Ponga su solución aquí <<<
  d = zeros(rows(pos),5);
  m = 1;
  n = 1;
  for m = 1:rows(pos)
    for n = 1:columns(pos)
      d(m,n) = sqrt(((pos(m,1)-emisorPos(1,n))^2) + ((pos(m,2)-emisorPos(2,n))^2) + ((pos(m,3)-emisorPos(3,n))^2));
    endfor
  endfor
  
endfunction


## #################################
## ## Pruebas con varios emisores ##
## #################################

## Pruebe todo el esquema usando el número de emisores dado
## n:    número de emisores que debe ser 3, 4 o 5
## D:    datos de distancia
## E:    posición de emisores
## fig1: número de primera figura a utilizar
function testCase(n,D,E,fig1)

  assert(n>2 && n<=columns(D));

  printf("Probando el caso de %i emisores\n",n);

  ## Prueba unitaria: calcule la posición del primer dato
  p=calcPosition(D(1,1:n),E)
  
  ## Calcule las posiciones a partir de las distancias
  p=calcPositions(D(:,1:n),E);

  ## Distancias a posiciones estimadas
  ed=calcDistances(p,E);

  ## Errores
  err=(ed-D);

  ## Promedio de errores
  error_average = sum(err)/rows(err)

  ## Muestre los errores de las distancias a los primeros tres emisores
  figure(fig1,"name",cstrcat("Errores con ",num2str(n)," emisores"));
  hold off;
  plot(err(:,1),"r");
  hold on;
  plot(err(:,2),"b");
  plot(err(:,3),"k");
  
  ## Grafique las posiciones encontradas
  figure(fig1+1,"name",...
	 cstrcat("Trayectoria estimada con ",num2str(n)," emisores"));
  hold off;
  
  ## Grafique la trayectoria en azul
  plot3(p(:,1),p(:,2),p(:,3),'b');
  xlabel("x");
  ylabel("y");
  zlabel("z");
  
  hold on;
  
  if (n==3)  
    ## Si n==3 hay otros posibles valores para la trayectoria

    ## Use - en la cuadrática
    p=calcPositions(D(:,1:n),E,2);
    plot3(p(:,1),p(:,2),p(:,3),'r');

    ## Use la predicción para auto-seleccionar cuál solución usar y
    ## píntela con cuadrados magenta
    p=calcPositions(D(:,1:n),E,3);
    plot3(p(:,1),p(:,2),p(:,3),'ms');
    
    ## A la fuerza use n=5 para ver la trayectoria verdadera
    ## y píntela usando circulos verdes
    p=calcPositions(D(:,1:5),E);
    plot3(p(:,1),p(:,2),p(:,3),'go');
  endif

  ## Grafique los emisores en postes
  plot3([E(1,1:n);E(1,1:n)],...
	[E(2,1:n);E(2,1:n)],...
	[E(3,1:n);zeros(1,n)],...
	"linewidth",3);
  
  plot3(E(1,1:n),E(2,1:n),E(3,1:n),'rx');
endfunction;

## Fuerce probar con 3, 4 y 5 emisores
testCase(3,D,E,1); ## 3 emisores, en las figuras 1 y 2
testCase(4,D,E,3); ## 4 emisores, en las figuras 3 y 4
testCase(5,D,E,5); ## 5 emisores, en las figuras 5 y 6
