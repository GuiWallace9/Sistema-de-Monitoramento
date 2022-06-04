#include <stdio.h>
#include <stdlib.h>
#define PI 3.14156
#define RAIO_AP 200
#define RAIO_ZA 2000
#define DELTA_ALARME 60
#define EPS_COS 0.000001
#define EPS 0.01



/*FUNCTION DECLARATION*/
int iguais(double x, double y);
double cosseno(double theta, double epsilon);
int localiza (double xi, double yi, double div, double xj, double yj, double djv, double xk, double yk, double dkv, double *xv, double *yv);
double raiz (double x, double epsilon);
double velocidade(double x0, double y0, double x1, double y1, double deltaT);
double distancia(double x0, double x1, double y0, double y1);
int intercepta(double x0, double y0, double x1, double y1, double *x, double *y);


/*FUNÇÃO PARA DETERMINAR SE DOIS NÚMEROS SÃO IGUAIS*/
int iguais(double x, double y) {
  if(x-y < EPS && y-x < EPS)
  return 1;
  else
  return 0;
}


/*FUNÇÃO PARA CALCULAR O COSSENO ATRAVÉS DA SÉRIE DE TAYLOS*/

double cosseno(double theta, double epsilon) {
  double cosseno = 1.0, fraction = 1.0, numerador = 1.0, denominador = 1.0, grau;
  int i,j, n=1;
  grau = theta * PI / 180;
  while (fraction >= epsilon) {
    for (i=2*n; i > 0; i--) {
      numerador = numerador * grau;
    }
    for (j=2*n;j > 0; j--) {
      denominador = denominador * j;
    }
    fraction = numerador / denominador;
    if (n%2) cosseno = cosseno - fraction;
    else cosseno = cosseno + fraction;
    n++;
    denominador = 1.0;
    numerador = 1.0;

  }
  return cosseno;

}


/*FUNÇÃO PARA LOCALIZAR A POSIÇÃO EXATA DO VEÍCULO*/

int localiza (double xi, double yi, double div, double xj, double yj, double djv, double xk, double yk, double dkv, double *xv, double *yv) {
    double pki, pkj, qki, qkj, pji, pjk, qji, qjk, pij, pik, qij, qik;
   if (iguais(xi,xj) && iguais(xi,xk)) {
     return 0;
   } else if (iguais(xi,xj)) {
     pki = (((xk*xk)-(xi*xi)+(yk*yk)-(yi*yi)-(dkv*dkv)+(div*div))/(2*(xk-xi)));
     pkj = (((xk*xk)-(xj*xj)+(yk*yk)-(yj*yj)-(dkv*dkv)+(djv*djv))/(2*(xk-xj)));
     qki = ((yi - yk)/(xk-xi));
     qkj = ((yj - yk)/(xk-xj));
     *yv = ((pki - pkj)/(qkj - qki));
     *xv = pki + qki* (*yv);
     return 1;
   } else if (iguais(xi,xk)) {
     pji = (((xj*xj)-(xi*xi)+(yj*yj)-(yi*yi)-(djv*djv)+(div*div))/(2*(xj-xi)));
     pjk = (((xj*xj)-(xk*xk)+(yj*yj)-(yk*yk)-(djv*djv)+(dkv*dkv))/(2*(xj-xk)));
     qji = ((yi - yj)/(xj-xi));
     qjk = ((yk - yj)/(xj-xk));
     *yv = ((pji - pjk)/(qjk - qji));
     *xv = pji + qji* (*yv);
     return 1;
   } else {
     pij = (((xi*xi)-(xj*xj)+(yi*yi)-(yj*yj)-(div*div)+(djv*djv))/(2*(xi-xj)));
     pik = (((xi*xi)-(xk*xk)+(yi*yi)-(yk*yk)-(div*div)+(dkv*dkv))/(2*(xi-xk)));
     qij = ((yj - yi)/(xi-xj));
     qik = ((yk - yi)/(xi-xk));
     *yv = ((pij - pik)/(qik - qij));
     *xv = pij + qij * (*yv);
     return 1;
   }
}

double raiz (double x, double epsilon) {
  double r0, ri;
  r0 = x;
  ri = (0.5)*(r0+(x/r0));
  while (abs(ri-r0) > epsilon) {
    r0 = ri;
    ri = (0.5)*(r0+(x/r0));
  }
  return ri;
}


/*FUNÇÃO PARA CALCULAR A VELOCIDADE*/

double velocidade(double x0, double y0, double x1, double y1, double deltaT) {
  double dist;
  dist = distancia(x0, x1, y0, y1);
  return (dist / deltaT);
}

/*FUNÇÃO PARA CALCULAR A DISTÂNCIA ENTRE DOIS PONTOS*/
double distancia(double x0, double x1, double y0, double y1) {
  double dist, d;
  d = (y1-y0)*(y1-y0) + (x1-x0)*(x1-x0);
  dist = raiz(d, EPS);
  if (x0 == x1 && y0 == y1) {
    return 0;
  } else {
    return dist;
  }
}

/*INTERCEPTA A AP*/

int intercepta(double x0, double y0, double x1, double y1, double *x, double *y) {
  double m, n, a, b, c, delta, d0, d1, pontoy1, pontoy2, pontox1, pontox2, d3, d4;
  if (iguais(x0,x1) && iguais(y0,y1)) {
    /*printf("O CARRO ESTÁ PARADO \n");*/
    return 0;
  }
  if (iguais(x0, x1)) {
    if (x0 > 200 || x0 < -200) {
      /*printf("A TRAJETÓRIA NÃO INTERCEPTA A ZONA AP\n");*/
      return 0;
    } else {
      d0 = distancia(0,x0,0,y0);
      d1 = distancia(0,x1,0,y1);
      if (d0 <= d1) {
        /*printf("A TRAJETÓRIA SE AFASTA DA ZONA AP\n");*/
        return 0;
      } else {
        /*printf("A TRAJETÓRIA SE APROXIMA DA ZONA AP \n");*/
        pontoy1 = raiz(RAIO_AP*RAIO_AP-x0*x0, EPS);
        pontoy2 = (-1)*raiz(RAIO_AP*RAIO_AP-x0*x0, EPS);
        d3 = distancia(x0,x1,pontoy1,y1);
        d4 = distancia(x0,x1,pontoy2,y1);
        if (d3 >= d4) {
          *x = x0;
          *y = pontoy2;

        } else {
          *x = x0;
          *y = pontoy1;
        }
        return 1;
      }
    }
  } else {
    m = (y1-y0) / (x1-x0);
    n = y0 - m * x0;
    a = m*m+1;
    b = 2*m*n;
    c = n*n-RAIO_AP*RAIO_AP;
    delta = b*b - 4*a*c;
    if (delta < 0) {
      /*printf("A TRAJETÓRIA SE AFASTA DA ZONA AP \n");*/
      return 0;
    } else if (iguais(delta,0)){
      d0 = distancia(0,x0,0,y0);
      d1 = distancia(0,x1,0,y1);
      if (d0 <= d1) {
        /*printf("A TRAJETÓRIA SE AFASTA DA ZONA AP \n");*/
        return 0;
      } else {
        *x = (-b)/(2*a);
        *y = m*(*x)+n;
        /*printf("A TRAJETÓRIA SE APROXIMA DA ZONA AP \n");*/
        return 1;
      }
    } else {
      d0 = distancia(0,x0,0,y0);
      d1 = distancia(0,x1,0,y1);
      if (d0 <= d1) {
        /*printf("A TRAJETÓRIA SE AFASTA DA ZONA AP \n");*/
        return 0;
      } else {
        pontox1 = ((-b)+raiz(delta,EPS)) / (2*a);
        pontox2 = ((-b)-raiz(delta,EPS)) / (2*a);
        pontoy1 = m*pontox1+n;
        pontoy2 = m*pontox2+n;
        d3 = distancia(pontox1,x1,pontoy1,y1);
        d4 = distancia(pontox2,x1,pontoy2,y1);
        if(d3 >= d4) {
          *x = pontox2;
          *y = pontoy2;
        } else {
          *x = pontox1;
          *y = pontoy1;
        }
        /*printf("A TRAJETÓRIA SE APROXIMA DA ZONA AP \n");*/
        return 1;
      }
    }
  }
}
int main() {
  FILE *arq;
  char filename[256];
  int cont, n, idveiculo, antenaI, antenaJ, antenaK, perigo, teste_localiza1, teste_localiza2;
  double xi, yi, hiv, thetaiv, div, xj, yj, hjv, thetajv, djv, xk, yk, hkv, thetakv, dkv, xv, yv, deltaT, dist, speed, distOrigem, x0, y0, x1, y1, x, y, dist_toque;
  printf("Digite o nome do arquivo de entrada: \n");
  printf("\n");
  printf("\n");
  scanf("%s", filename);

  arq = fopen(filename, "r");

  if(arq == NULL) {
    printf("Não consegui abrir o arquivo %s. \n", filename);
    return 0;
  }
  /*LER A QUANTIDADE DE VEÍCULOS*/
  fscanf(arq, "%d", &n);
  printf("Número de casos a serem analisados: %d \n", n);
  printf("\n");
  printf("\n");
  cont = 0;
  while (cont < n) {
    printf("\n");
    fscanf(arq, "%d", &idveiculo);
    printf("\n");
    fscanf(arq, "%d %lf %lf %lf %lf", &antenaI, &xi, &yi, &hiv, &thetaiv);
    div = hiv * cosseno(thetaiv, EPS_COS);
    fscanf(arq, "%d %lf %lf %lf %lf", &antenaJ, &xj, &yj, &hjv, &thetajv);
    djv = hjv * cosseno(thetajv, EPS_COS);
    fscanf(arq, "%d %lf %lf %lf %lf", &antenaK, &xk, &yk, &hkv, &thetakv);
    dkv = hkv * cosseno(thetakv, EPS_COS);
    fscanf(arq, "%lf", &deltaT);
    teste_localiza1 = localiza(xi,yi,div,xj,yj,djv,xk,yk,dkv,&xv,&yv);

    x0 = xv;
    y0 = yv;
    if (!teste_localiza1) {
      printf("\n");
      printf("IDENTIFICAÇÃO: Veículo %d \n", idveiculo);
      printf("\n");
      printf("id  |       posição      |   H(m)   |  theta(graus) | distância(m) | \n");
      printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaI, xi, yi, hiv, thetaiv, div);
      printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaJ, xj, yj, hjv, thetajv, djv);
      printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaK, xk, yk, hkv, thetakv, dkv);
      printf("\n");
      printf("Nao foi possivel calcular a localizacao inicial do veiculo %d. \n", idveiculo);
    } else {
      printf("\n");
      printf("IDENTIFICAÇÃO: Veículo %d \n", idveiculo);
      printf("\n");
      printf("Antenas na posição prévia: \n");
      printf("id  |       posição      |   H(m)   |  theta(graus) | distância(m) | \n");
      printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaI, xi, yi, hiv, thetaiv, div);
      printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaJ, xj, yj, hjv, thetajv, djv);
      printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaK, xk, yk, hkv, thetakv, dkv);
      printf("\n");
      printf("Localização pŕevia (%.2lf, %.2lf) \n", x0, y0);
      printf("\n");
      printf("Intervalo de tempo: %.2lf segundos \n", deltaT);
      printf("\n");
    }


    fscanf(arq, "%d %lf %lf %lf %lf", &antenaI, &xi, &yi, &hiv, &thetaiv);
    div = hiv * cosseno(thetaiv, EPS_COS);
    fscanf(arq, "%d %lf %lf %lf %lf", &antenaJ, &xj, &yj, &hjv, &thetajv);
    djv = hjv * cosseno(thetajv, EPS_COS);
    fscanf(arq, "%d %lf %lf %lf %lf", &antenaK, &xk, &yk, &hkv, &thetakv);
    dkv = hkv * cosseno(thetakv, EPS_COS);
    teste_localiza2 = localiza(xi,yi,div,xj,yj,djv,xk,yk,dkv,&xv,&yv);
    x1 = xv;
    y1 = yv;
    if (!teste_localiza2 || !teste_localiza1) {
      printf("\n");
      if (!teste_localiza1) {
        printf("\n");
      } else {
        printf("id  |       posição      |   H(m)   |  theta(graus) | distância(m) | \n");
        printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaI, xi, yi, hiv, thetaiv, div);
        printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaJ, xj, yj, hjv, thetajv, djv);
        printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaK, xk, yk, hkv, thetakv, dkv);
        printf("\n");
        printf("Nao foi possivel calcular a localizacao final do veiculo %d. \n", idveiculo);
      }

    } else {
      printf("\n");
      printf("Antenas na posição atual: \n");
      printf("\n");
      printf("id  |       posição      |   H(m)   |  theta(graus) | distância(m) | \n");
      printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaI, xi, yi, hiv, thetaiv, div);
      printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaJ, xj, yj, hjv, thetajv, djv);
      printf("%d | ( %.2lf , %.2lf )| %.2lf | %.2lf | %.2lf \n", antenaK, xk, yk, hkv, thetakv, dkv);
      printf("\n");
      printf("Localização atual (%.2lf, %.2lf) \n", x1, y1);
      printf("\n");
    }


    /*CHAMADA DE FUNÇÕES*/
    perigo = intercepta(x0, y0, x1, y1, &x, &y);
    dist = distancia(x0,x1,y0,y1);
    speed = velocidade(x0,y0,x1,y1,deltaT);
    distOrigem = distancia(0,x1,0,y1);

    if (!teste_localiza2 || !teste_localiza1) {
      printf("\n");
      printf("Nao foi possivel determinar a situacao do veiculo %d", idveiculo);
      printf("\n");
    } else {
      printf("\n");
      printf("Distancia percorrida: %.2lf m \n", dist);
      printf("Velocidade: %.2lf m/s \n", speed);

      printf("\n");
      printf("Distancia da origem: %.2lf \n", distOrigem);
      if(distOrigem <= RAIO_AP) {
        if (dist == 0) {
          printf("\n");
          printf("Veículo está estacionado na zona AP \n");
          printf("\n");
          printf("********************************************************************************* \n");
          printf("ALERTA| ALERTA| ALERTA| ALERTA| ALERTA| ALERTA| ALERTA: \n");
          printf("Veículo está na AP \n");
          printf("********************************************************************************* \n");
          printf("\n");

        } else {
          printf("\n");
          printf("Veículo em movimento na zona AP \n");
          printf("\n");
          printf("********************************************************************************* \n");
          printf("ALERTA| ALERTA| ALERTA| ALERTA| ALERTA| ALERTA| ALERTA: \n");
          printf("Veículo está em movimento na AP \n");
          printf("********************************************************************************* \n");
          printf("\n");
        }

      } else if (distOrigem > RAIO_AP && distOrigem <= RAIO_ZA) {
        if (dist == 0) {
          printf("\n");
          printf("Veículo está estacionado na zona de ALERTA \n");
          printf("\n");
        } else {
          printf("\n");
          printf("Veículo está em movimento na zona de ALERTA \n");
          printf("\n");
        }

      } else {
        if (dist == 0) {
          printf("\n");
          printf("Veículo está estacionado FORA da zona de alerta \n");
          printf("\n");
        } else {
          printf("\n");
          printf("Veículo está em movimento FORA Da zona de alerta \n");
          printf("\n");
        }
      }


      if (perigo && distOrigem <= RAIO_ZA / 4) {
        printf("\n");
        dist_toque = distancia(x,x1,y,y1);
        printf("TRAJETORIA INTERCEPTARA A ZONA AP \n");
        printf("Distância atual à AP: %.2lf \n", dist_toque);
        printf("Interseção ocorrerá em %.2lf segundos \n", dist_toque/speed);
        printf("O toque ocorrerá na coordenada (%.2lf , %.2lf) \n", x, y);
        printf("\n");
        printf("********************************************************************************* \n");
        printf("ALERTA| ALERTA| ALERTA| ALERTA| ALERTA| ALERTA| ALERTA: \n");
        printf("           RISCO IMINENTE \n");
        printf("Carro está em direção a zona AP \n");
        printf("********************************************************************************* \n");
        printf("\n");


      } else if (perigo && distOrigem <= RAIO_ZA / 2) {
        printf("\n");
        dist_toque = distancia(x,x1,y,y1);
        printf("TRAJETORIA INTERCEPTARA A ZONA AP \n");
        printf("Distância atual à AP: %.2lf \n", dist_toque);
        printf("Interseção ocorrerá em %.2lf segundos \n", dist_toque/speed);
        printf("O toque ocorrerá na coordenada (%.2lf , %.2lf) \n", x, y);
        printf("\n");
        printf("********************************************************************************* \n");
        printf("ALERTA| ALERTA| \n");
        printf("Carro está em direção a zona AP \n");
        printf("********************************************************************************* \n");
        printf("\n");
      } else if (perigo && distOrigem <= RAIO_ZA ) {
        printf("\n");
        dist_toque = distancia(x,x1,y,y1);
        printf("TRAJETORIA INTERCEPTARA A ZONA AP \n");
        printf("Distância atual à AP: %.2lf \n", dist_toque);
        printf("Interseção ocorrerá em %.2lf segundos \n", dist_toque/speed);
        printf("O toque ocorrerá na coordenada (%.2lf , %.2lf) \n", x, y);
        printf("\n");
        printf("************************************* \n");
        printf("ALERTA| ALERTA| \n");
        printf("Risco Leve: Carro está em direção a zona AP \n");
        printf("************************************* \n");
        printf("\n");

      }
    }

    cont++;
  }
  fclose(arq);
  return 1;
}
