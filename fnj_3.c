#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <iostream.h>

struct p_node {
	long int d;
	int i,j;
};

struct otuName {
   char *name;
   int ordem;
};

typedef struct p_node p_node_type;


int father(int x){
    return x/2; }

void heap_fica(p_node_type vet[], int i, int qtde){
   int g,j=i;
   p_node_type aux;
   while( (2*j)<=qtde ){
      g=2*j;
      if ((g<qtde) && (vet[g].d > vet[g+1].d))
         g=g+1;
      if (vet[j].d <= vet[g].d)
         j=qtde;
      else {
         aux = vet[j];
         vet[j] = vet[g];
         vet[g] = aux;
         j=g;
      }
   }
   return;
}
        
void matrizQ(long int **nm, int sz, long int **d, long int D[]){
   int i,j,x=0;
   for(i=0;i<sz;i++){
      for(j=0;j<i;j++){
         nm[i][j]= ((long int)(sz - 2) * d[i][j]) - (D[i] + D[j]);
//         if(x<45) printf("(size-2) * m[i=%d] - didj = %d * %ld - %ld = %ld\n",x,(sz-2),d[i][j],(D[i]+D[j]),nm[i][j]);
         x++;
      }
   }
   //printf("MATRIZ Q\n");
   //for(i=0;i<sz;i++){
   //   for(j=0;j<i;j++)
   //     printf("%ld ",nm[i][j]);
   //   printf("\n");}

   return;
}//matrizQ()

void minimoModificado(p_node_type R[], int size, int nelem, long int **mq, double p){
   // esse metodo obtem uma lista R de n/p  pares de OTUs que possuem menor valor Q
   // carregar a matriz Q em uma lista L = (i, j, Q[i][j]) onde i, j sao indices de Q
   //printf("minimoModificado size=%d n.elem=%d\n",size,nelem);
   p_node_type aux;
   p_node_type *L;
   L=(p_node_type *) malloc((nelem+1)*sizeof(p_node_type));
   //printf("1 ");
   //p_node_type L[nelem+1];
   int indices[size];
   int k, lim, j=0, ind=0,x;
   lim = size*p;
   for (int i=0; i<size; i++) indices[i]=0;
   for (int i=0; i<=nelem; i++){
      L[i].i=-1;
      L[i].j=-1;
      L[i].d=0; }
   //printf ("\nModificado size = %d p=%f nelem=%d lim=%d size*p=%f\n",size,p,nelem,lim,size*p);
   k=1;
   for(int i=0; i<size; i++)						//monta uma lista L com as dist. em m
      for(int j=0; j<i; j++){ 
         L[k].i=i; 
         L[k].j=j; 
         L[k++].d=mq[i][j];
      }
   k--;
   // transformar L em um heap
   x=k/2;//floor((k-1)/2);
   clock_t t1,t2;
   t1=clock();	
   for (int h=x; h>0; h--){
      heap_fica(L,h,k);
   }
   t2=clock();
   printf("Ordenação %f s\n",(t2-t1)/(double)CLOCKS_PER_SEC);
   //printf("heap montado \n");
   //for (int h=1; h<k; h++) printf("L[%d]=(%d,%d)=%ld \n",h,L[h].i,L[h].j,L[h].d);
               
   // obter a lista de nos l[u] = True se o no u foi escolhido; False caso contrario.
   if (lim < 1) lim=1 ;
   //printf("\np=%f lim=%d   size=%d  nelem=%d  k=%d\n",p,lim,size,nelem,k);
   t1=clock();
   int i=0;	
   if(lim==1){
      R[0]=L[1];}
//   printf("\nR[0]=(%d,%d)=%ld\n",R[0].i,R[0].j,R[0].d);}
   else{
      i=0;
      //seleciona os menores o topo do heap 
      while((i<lim)&&(nelem>0)){
         if(indices[L[1].i]==0 && indices[L[1].j]==0){
//            printf("Selecionou %d (%d,%d)=%ld \n",i,L[1].i,L[1].j, L[1].d);
            R[i++]=L[1];
            indices[L[1].i]=1; 
            indices[L[1].j]=1; 
         }
         /*remover (i+1) do heap*/
         //printf("0,");
         L[1]=L[nelem];
         nelem--;
         //printf("1,");
         heap_fica(L, 1, nelem);
         //printf("2\n");
      }//while()
   }
   t2=clock();
   printf("Sel.Pares %f s\n",(t2-t1)/(double)CLOCKS_PER_SEC);
   //printf("\nMinimos selecionados= ");
   //for (int y=0;y<lim;y++) printf("R[%d]=(%d,%d)=%ld , ",y,R[y].i, R[y].j,R[y].d);
   //printf("\n");
   if(L) free(L);
   return;
}//minimo modificado


// Driver Code
int main(int argc, char *argv[])
{
   // *************************************ABRINDO ARQUIVO COM MATRIZ DE DISTÂNCIA
   FILE* ptr = fopen(argv[1], "r");
   double f = strtod(argv[2],NULL);

   if (ptr == NULL) {
      printf("no such file.");
      return 0;
   }
	
   int  size, size_O,i=0,j=0,k=1, ind=0;
   fscanf(ptr,"%d", &size);		      	// ****LENDO No. OTUs
   size_O=size;
   float s = 0.0;
   s = (size/2);
   s = (size+1)*s;  			//******DETERMINANDO TAMANHO MATRIZ TRIANGULAR
   int t=s-size;

   int limite = size*f;
   if(limite==0 && (size*f)>0.0) limite=1;

   printf("FNJ_SEQ SIZE=%d f=%1.4f ",size,f);
   printf("Tam Matriz t=%d Tam. Heap = %d\n",t, limite);
     

   // ****  ESTRUTURAS DE DADOS GLOBAIS
   struct otuName otus[size]; // NOMES DAS OTUS
   long int **m;				// MATRIZ DE DISTANCIAS
   long int *Q[size];      // MATRIZ Q
   long int D[size];			// SOMA DAS DISTANCIAS DE UMA OTU(linha)
   p_node_type h[limite];  // HEAP
   
   //*****  ALOCAÇÃO DE MEMÓRIA
   m=(long int**)malloc(size*sizeof(long int));   
   for(int i=1;i<size;i++)
      Q[i]= (long int *) malloc(i*sizeof(long int));

   int menor,maior;
   long int mij, Di, Dj;

   int  z=0;    
   for(i=0;i<size;i++){			// LAÇO ALOCAÇÃO LINHAS DA MATRIZ E INICIALIZAÇÃO DE VALORES
      otus[i].name = (char *) malloc(12*sizeof(char));
      D[i]=0;
      m[i]=(long int*)malloc((i+1)*sizeof(long int));

      for(j=0;j<i;j++) {m[i][j]=-1; z++;}
      sprintf(otus[i].name,"%d",i);
      otus[i].ordem=i;
   }
   
   //printf("Matriz Triangular:\n");
   i=0;k=0;
   while (i<size){				// LENDO DISTANCIAS DO ARQUIVO...
      j=0;
      while(j<=i){
         fscanf(ptr, "%ld ", &m[i][j]);		// ...inserindo na matriz
         //if(m[i][j]==1238661257573) printf("i=%d j=%d m[i][j]=%ld\n", i , j, m[i][j]);
         j++;
     }
     i++;
   }
//printf("m[0][0]=%ld\n",m[0][0]);
   printf("\nFIM Leirura da Matriz Triangular: i=%d j=%d\n", i, j);

   for(i=0;i<size;i++){					// calculando soma das distancias das OTUS - D[i]
      for(j=0;  j<i;   j++) D[i]=D[i]+m[i][j];
      for(j=i+1;j<size;j++) D[i]=D[i]+m[j][i];
   }
   //printf("Soma das Distacias na liha i:\n");
   //for(i=0;i<size;i++) printf("D[%d] = %ld \n",i,D[i]);

   //**************************************** LAÇO PRINCIPAL
   clock_t tempo=clock();
   clock_t t1,t2, t_i, t_f;
   int inter=0;
   
   while(size>=2){
      t_i=clock();
      if (size==2){
	//char aux[(strlen(otus[0].name)+strlen(otus[1].name)+3)];
   // 
   // if(otus[0].ordem<otus[1].ordem)
	//       sprintf(aux,"(%s,%s)",otus[0].name,otus[1].name);
	//    else
	//        sprintf(aux,"(%s,%s)",otus[1].name,otus[0].name);
	//    //sprintf(otus[0].name,"%s",aux);
	//    printf("ÁRVORE = %s\n",aux); //otus[0].name);
         size--;//exit(1);
	//free(aux);
     }
     else
     {
        printf("\n________________________________________________Size=%d\n",size);

//        int v1=0;
//        if(inter<2)
//        for(int k1=1;k1<size;k1++)
//           for(int j1=0; j1<k1;j1++)
//            printf("i=%d j=%d m[%d]=%ld \n",k1,j1,v1++,m[k1][j1]);

//        if(inter<2)
//        for(int k1=0;k1<size;k1++)
//            printf("D[%d]=%ld \n",k1,D[k1]);

        s = 0.0;
        s = size/2.0;
        s = s*(size+1);					//******DETERMINANDO TAMANHO MATRIZ TRIANGULAR
        limite = size*f;					// número de pares a serem selecionados
        if (limite==0) limite=1;
        t = s-size;
        //printf("Calculando Matriz Q %d\n",size);
        t1=clock();
        matrizQ(Q, size, m, D);								// CALCULANDO MATRIZ Q
        t2=clock();
        printf("Matriz Q: %f s\n",(t2-t1)/(double)CLOCKS_PER_SEC);
        int z=0;
//        if (inter==0) {
//           printf("\n\nCálculo MatrizQ        t=%f \n",(t2-t1)/(double)CLOCKS_PER_SEC);

//           for(int x=0;x<10;x++) 
//              for(int y=0;y<x;y++)
//                 printf("i=%d (%d,%d)=%ld \n",z++,x,y,Q[x][y]);
//        }
//        printf("\n\n");*/

        t1=clock();
        minimoModificado(h,size, t, Q, f);      // Seleciona Menores
        t2=clock();
        //if (inter==0) printf("Selecionar menores     t=%f \n",(t2-t1)/(double)CLOCKS_PER_SEC);
        
//        for(int z=0;z<10;z++) 
//           printf("h[%d]=(%d,%d)=%ld\n",z,h[z].i,h[z].j,h[z].d);
           
        printf("Juntando %d Pares de OTU's\n",limite);
        t1=clock();
        for(int c=0; c<limite; c++) {
           clock_t t_11,t_21;
           t_11=clock();
           if (size > 2) {
              double d_i_novo,d_j_novo=0.0;
              double t1, t2;
              if (h[c].i > h[c].j) {
                 //maior i
                 mij=m[h[c].i][h[c].j];
                 Di=D[h[c].j]; 
                 Dj=D[h[c].i];
                 menor=h[c].j;
                 maior=h[c].i; 
              } 
              else {
                 //maior j
                 mij = m[h[c].j][h[c].i];
                 Di = D[h[c].i]; 
                 Dj = D[h[c].j];
                 menor=h[c].i;
                 maior=h[c].j; 
              }
//              printf("mij=%ld didj=%ld menor=%d maior=%d\n",mij,Di+Dj, menor, maior);
              //char *aux=(char *) malloc((strlen(otus[menor].name)+strlen(otus[maior].name))*sizeof(char));
              //int x = (strlen(otus[menor].name)+ strlen(otus[maior].name)+3)%8;
              //if (x>0) x = 8-x; 
              char aux[(strlen(otus[menor].name)+strlen(otus[maior].name)+4)];
              if(otus[menor].ordem  < otus[maior].ordem)
                 sprintf(aux,"(%s,%s)",otus[menor].name,otus[maior].name);
              else{
                 sprintf(aux,"(%s,%s)",otus[maior].name,otus[menor].name);
                 otus[menor].ordem=otus[maior].ordem;
              }
                  
              //printf("\n%ld ",(strlen(otus[menor].name)+strlen(otus[maior].name)+3));
              //printf("name = %s tamenho=%ld menor=%d\n",otus[menor].name,strlen(otus[menor].name),menor);
              //printf("name = %s tamenho=%ld\n",otus[maior].name,strlen(otus[maior].name));
              //printf("aux = %s tamanho=%ld\n",aux,strlen(aux));
              free(otus[menor].name);
                  
             //else printf("Não Liberou*******\n\n");
             //if(strlen(aux)>strlen(otus[menor].name)){  printf("R ");
             int y;
             y = 8 - (strlen(aux)%8) + 1;
             //printf("%d ",y);
             // otus[menor].name=(char*)realloc(otus[menor].name,(strlen(aux)+y));
             if(!(otus[menor].name=(char*) malloc((strlen(aux)+y)*sizeof(char))))
                printf("Não alocou*************\n");
             strcpy(otus[menor].name, aux);
             printf("\njuntando c=%d otus[menor=%d].name=%s mij=%ld\n",c,menor,otus[menor].name, mij);
               
             //printf("menor=%d maior=%d size=%d\n",menor,maior,size);
             //printf("\n Nomes\n");
             //for(int k=0;k<size;k++) printf("%d- %s\n",otus[k].ordem,otus[k].name);
             //printf("Calculando distancias da nova OTU para as demais.\n"); 
             t1=(size - 2) * 2;
             t2= (Di - Dj);
             t2 = t2/t1;
             //printf("t1 = %f  t2=%f \n",t1,t2);
             d_i_novo = (0.5*mij);
             d_i_novo = d_i_novo + t2 ;
             d_j_novo = mij - d_i_novo;
             //printf("calculando distancia do novo aos demais ---- calculando a linha-coluna\n");
             //for(int k=0; k<size; k++) printf("D[%d]=%ld\n",k,D[k]);
             //printf("menor=%d maior=%d\n",menor,maior);
             
             //printf("Distâncias (menor=%d,maior=%d\n",menor,maior);
             t_11=clock();
             for(int k=0;k<menor;k++){
               if ( ((menor==2)&&(k==3)) || ((menor==3)&&(k==349)) )
                printf("A1-m[%d][%d] = ((m[%d][%d]=%ld + m[%d][%d])=%ld - mij=%ld)*0.5=%ld\n",menor,k, menor,k, m[menor][k], maior, k, m[maior][k], mij, (long int)(((m[menor][k]+m[maior][k])-mij)*0.5));
                D[menor] = D[menor] - (m[menor][k]);//+m[maior][k]);
                D[k] = D[k] - (m[menor][k]+m[maior][k]);//+m[maior][k]);
                m[menor][k] = (long int)((m[menor][k]+m[maior][k])-mij)*0.5;
                D[k] = D[k] + m[menor][k];
                D[menor] = D[menor] + m[menor][k];
             }
             for(int k=menor+1;k<maior;k++){
               if ( ((menor==2)&&(k==3)) || ((menor==3)&&(k==349)) )
                printf("A2-m[%d][%d]= (m[%d][%d]=%ld + m[%d][%d]=%ld - mij=%ld)*0.5=%ld\n",k,menor,k,menor,m[k][menor],maior,k,m[maior][k],mij, (long int)(((m[k][menor]+m[maior][k])-mij)*0.5));
                D[menor] = D[menor]-m[k][menor]; //-----
                D[k] = D[k] - (m[k][menor]+m[maior][k]);//+m[maior][k]);
                m[k][menor]= (long int)((m[k][menor]+m[maior][k])-mij)*0.5;
                D[k] = D[k] + m[k][menor];
                D[menor]=D[menor]+m[k][menor];//---
             }
             for(int k=maior+1;k<size;k++){
               if ( ((menor==2)&&(k==3)) || ((menor==3)&&(k==349)) )
                printf("A3-m[%d][%d] = ((m[%d][%d]=%ld + m[%d][%d])=%ld - mij=%ld)*0.5=%ld\n",k,menor,k,menor,m[k][menor],k,maior,m[k][maior],mij, (long int)(((m[k][menor]+m[k][maior])-mij)*0.5));
                D[menor] = D[menor]-m[k][menor]; //-----
                D[k] = D[k] - (m[k][menor]+m[k][maior]);
                m[k][menor] = (long int)((m[k][menor]+m[k][maior])-mij)*0.5;
                D[k] = D[k] + m[k][menor];
                D[menor]=D[menor]+m[k][menor];//---
             }
             t_21=clock();
//             if (c==0 && inter==0) printf("Distância para demais   t=%f \n",(t_21-t_11)/(double)CLOCKS_PER_SEC);
           if ((menor==2))
              printf("D-m[3][2]= %ld\n", m[3][2]);
           if ((menor==3))
              printf("D-m[349][3]= %ld\n", m[349][3]);

//           printf("D-m[v_somai[159]+9].i=159 m[v_somai[159]+9].j=9 \n");

             //printf("Eliminando linha-coluna com cópia.\n");
//           if((menor==2) || (menor==3)) printf("Copia linha (maior=%d,última=%d\n",maior,size-1);
//   for(int j1=0;j1<maior;j1++)			//percorre a linha maior até coluna maior copiando da última
//      printf("m[maior=%d][j=%d]=m[ultima=%d][j=%d]=%ld\n",maior, j1, size-1,j1, m[size-1][j1]);


   
//   for(int i1=maior+1;i1<size-1;i1++){              //percorre linhas de maior até a penúltima copiando da uĺtima
//      printf("m[i=%d][maior=%d]=m[ultima=%d][i=%d]=%ld\n",i1, maior, size-1, i1, m[size-1][i1]);
//   }

             t_11=clock();
             for(int k=0;k<size;k++) 
                if (k>maior) m[k][maior]=m[size-1][k];
                else 
                   if(k==maior) m[maior][k]=0;
                   else m[maior][k]=m[size-1][k];
             t_21=clock();
//             if (c==0 && inter==0) printf("Eliminando linha-coluna t=%f \n",(t_21-t_11)/(double)CLOCKS_PER_SEC);
             free(otus[maior].name);
             otus[maior].name = otus[size-1].name;
             otus[maior].ordem=otus[size-1].ordem;
             D[maior]=D[size-1];
//           if ((menor==2))
//              printf("L-m[3][2] = %ld\nL-m[349][3] = %ld\n", m[3][2], m[349][3]);
//           if ((menor==3))
//              printf("L-m[3][2] = %ld\nL-m[349][3] = %ld\n", m[3][2], m[349][3]);
//              printf("L-m[v_somai[159]+9].i=159 m[v_somai[159]+9].j=9\n");
//             for(int k1=0;k1<10;k1++)
//                printf("m[maior,%d]=%ld = m[size-1,%d]=%ld\n",k1,m[maior][k1],k1,m[size-1][k1]);

             //printf("Corrigindo coordenadas.\n");
             for(int e=c+1; e<limite; e++){
                if(h[e].i==size-1)
                  if (h[e].j > maior) {
                     h[e].i=h[e].j;
                     h[e].j=maior;}
                  else h[e].i=maior;
               
                  if (h[e].j==size-1)
                     h[e].j =maior;
             }   
             t_21=clock();
//             if (c==0 && inter==0) printf("Corrigindo Coordenadas t=%f \n",(t_21-t_11)/(double)CLOCKS_PER_SEC);
              
             //printf("Recalculando somas das linhas.\n\n");
//             printf("ANTES D[MAIOR}=%ld\n",D[maior]);
//             D[maior]=0;
//             for(int k=0;k<size-1;k++)
//                D[maior]=D[maior]+m[size-1][k];
//             printf("DEPOIS D[MAIOR}=%ld\n",D[maior]);
//             if (c==0 && inter==0) printf("Demais tarefas        %8.4f\n",((double)(t2-t1))/(double)CLOCKS_PER_SEC);
          }//if (size > 2)	
          size--;
          //printf("\n++++++>>>>>>>>>>>SIZE=%d\n",size);
       }//for
       t2=clock();
       printf("Jun.Pares %f s\n",((double)(t2-t1))/(double)CLOCKS_PER_SEC);
       //free(h);
     }
     t_f=clock();
     printf("Iteração  %f s\n",((double)(t_f-t_i))/(double)CLOCKS_PER_SEC);
     inter++;
   }//while(size>1)
   float tempo2 = (clock() - tempo) / (double)CLOCKS_PER_SEC;
       
   if (otus[0].ordem<otus[1].ordem)
      printf("(%s,%s)\n",otus[0].name,otus[1].name);
   else
      printf("(%s,%s)\n",otus[1].name,otus[0].name);
      
     
   printf("Tempo(s) = %f \n",tempo2);
   
 
   free(m);
 //  free(Q);
 

   return 0;
}
