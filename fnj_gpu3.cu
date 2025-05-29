#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENACC
    #include <openacc.h>
#endif
//#include "fnjseq.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX_THREADS_PER_BLOCK 1024
#define THRESHOLD 512


struct p_node {
	long int d;
	int i,j,pv;
};

struct otuName {
   char *name;
   int ordem;
};

typedef struct p_node p_node_type;

struct m_distance{
   long int d;
   int i,j;
};
typedef struct m_distance m_dist;

int power2(int n){
   int i=1;
   while( i*2 <= n )
      i*=2;
   return(i);
}

void ordenacao(p_node_type *A, int inicio, int fim){
    int i,j;
    p_node_type key;
    //printf("\n\nORDENACAO %d - %d \n",inicio, fim);
    i=inicio;
    while(i<fim-1){

       j=i+1;
       key = A[j];
       //printf("%d, ",i);       
       while((j>0) && (A[j-1].d > key.d)) {
          //if(i>=523712) printf("j=%d ,",j);
          A[j]=A[j-1];
          j--;
       }
       A[j]=key;
       i++;
    }
    //printf("\n\n");
    return;
}

void matrizQ(p_node_type *q, int sz, m_dist *m, long int *D){
   int i,j, v=0;
   for(i=0;i<sz;i++)
      for(j=0;j<i;j++){
        q[v].i=i;
        q[v].j=j;
        q[v].pv=v;
        q[v].d =((long int)(sz - 2) * m[v].d) - (D[i] + D[j]);
        if(i>0) v++;
      }
   return;
}//matrizQ()

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
}//heap_fica


void minimoModificado(p_node_type R[], int size, int nelem, p_node_type *q, double f, int iter){
   // esse metodo obtem uma lista R de size*f  pares de OTUs que possuem menor valor q
   // carregar a matriz Q em uma lista L = (i, j, Q[i][j]) onde i, j sao indices de Q
   //printf("minimoModificado size=%d n.elem=%d\n",size,nelem);
   p_node_type L[nelem+1];
   int indices[size];
   int k, lim,x,v=0;
   lim = size*f;
   for (int i=0; i<size; i++) indices[i]=0;
   for (int i=0; i<=nelem; i++){
      L[i].i=-1;
      L[i].j=-1;
      L[i].pv=0;
      L[i].d=0; 
   }
   //printf ("\nModificado size = %d p=%f nelem=%d lim=%d size*p=%f\n",size,p,nelem,lim,size*p);
   k=1;
   for(int i=0; i<size; i++)						//monta uma lista L com as dist. em m
      for(int j=0; j<i; j++){ 
         L[k].i=i; 
         L[k].j=j; 
         L[k].pv=k;
         L[k++].d=q[v++].d;
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
   //printf("\nR[0]=(%d,%d)=%ld\n",R[0].i,R[0].j,R[0].d);}
   else{
      i=0;
      while((i<lim)&&(nelem>0)){		//seleciona os menores o topo do heap 
         //printf("i=%d nelem=%d L[1].i=%d L[1].j=%d \n",i,nelem,L[1].i,L[1].j);
         if(indices[L[1].i]==0 && indices[L[1].j]==0){
            printf("Selecionou i=%d (%d,%d)=%ld\n",i,L[1].i,L[1].j, L[1].d);
            R[i++]=L[1];
            indices[L[1].i]=1; 
            indices[L[1].j]=1; }
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
   //if(L) free(L);
   return;
}//minimo modificado


inline cudaError_t checkCuda(cudaError_t result)
{
  #if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
  #endif
  return result;
}

//GPU Kernel Implementation of Matrix Q computation
__global__ void matrizQ( int size, int t, m_dist *m, p_node_type *q, p_node_type *d, int pw ) {
    int v = blockIdx.x*blockDim.x + threadIdx.x;    // handles the data at its thread id
    while (v < pw) {
       if (v<t){
          q[v].d = (size-2)*m[v].d - d[v].d;
          q[v].i = d[v].i;
          q[v].j = d[v].j;
          q[v].pv = d[v].pv;
       }
       else{
          q[v].d = 99999999999999999;//INT_MIN;
          q[v].i = -1;
          q[v].j = -1;
          q[v].pv = -1;
       }
       v+= (blockDim.x * gridDim.x);
    }
    __syncthreads();
}//matrizQ


//GPU Kernel Implementation of Bitonic Sort
__global__ void bitonicSortGPU(p_node_type* arr, int j, int k)
{
    unsigned int i, ij;
    p_node_type temp;
    
    i = threadIdx.x + blockDim.x * blockIdx.x;

    ij = i ^ j;

    if (ij > i)
    {
        if ((i & k) == 0)
        {
            if (arr[i].d > arr[ij].d)
            {
                temp = arr[i];
                arr[i] = arr[ij];
                arr[ij] = temp;
            }
        }
        else
        {
            if (arr[i].d < arr[ij].d)
            {
                temp = arr[i];
                arr[i] = arr[ij];
                arr[ij] = temp;
            }
        }
    }
}


__global__ void init_Q(int n, p_node_type *vet) {
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride) {
    vet[i].i=-1; 
    vet[i].j=-1;
    vet[i].pv=-1;
    vet[i].d=-1;
  }
//  __syncthreads();
}//init


__global__ void init_D(int n, long int *vet) {
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i<n; i+=stride) {
    vet[i]=0;
  }
  //__syncthreads();
}//init_D

__global__ void init_m(int n, m_dist *vet) {
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i<n; i+=stride) {
    vet[i].d=0;
    vet[i].i=0;
    vet[i].j=0;
  }
  //__syncthreads();
}//init




__global__ void copyLine(int maior, int ultima, int size, m_dist *m, int *si){
   int stride = blockDim.x * gridDim.x;
   int index = threadIdx.x + blockIdx.x * blockDim.x;
   int si_i = si[maior];
   int si_ult = si[ultima];

   for(int j=index;j<maior;j+=stride)			//percorre a linha maior até coluna maior copiando da última
      m[si_i+j].d=m[si_ult+j].d;
	    
   __syncthreads();
   
   for(int i=maior+1;i<size-2;i+=stride){              //percorre linhas de maior até a penúltima copiando da uĺtima
      si_i=si[i];
      m[si_i+maior].d=m[si_ult + i].d;
   }
   __syncthreads();
}//copyLine


__global__ void copyLine_1(int maior, int ultima, int t_vet, int size, m_dist *m, int *si){
   int stride = blockDim.x * gridDim.x;
   int index = threadIdx.x + blockIdx.x * blockDim.x;

   int si_maior = si[maior];
   int si_ult = si[ultima];
   int l_maior = si[maior];
   int l_ult = si[ultima];
   int l_k=0, j_k=0;
   
   for(int k=index;k<t_vet;k+=stride){
      l_k=m[k].i;
      j_k=m[k].j;
      
      if ((l_k==l_maior)&&(j_k<maior))			//se linha igual a maior e coluna vai até maior-1
         m[k].d=m[si_ult+j_k].d;			//   copia da ultima linha até coluna maior para linha maior na coluna
      
      if(((l_k>l_maior)&&(l_k<size))&&(j_k==maior))	//se linha maior que maior e coluna igual maior  
         m[k].d=m[si_ult+j_k].d;			//   copia da ultima linha coluna maior para linha k coluna maior
     
   }
   __syncthreads();
}//copyLine



void imprimeCopiaLinhas(int maior, int ultima, int size, m_dist *m, int *si){   
   int si_i = si[maior];
   int si_ult = si[ultima];

   for(int j=0;j<maior;j++)			//percorre a linha maior até coluna maior copiando da última
      printf("j=%d m[si_i+j=%d].d=m[si_ult+j=%d].d=%ld\n",j, si_i+j,si_ult+j,m[si_ult+j].d);
	    
   
   for(int i=maior+1;i<size-2;i++){              //percorre linhas de maior até a penúltima copiando da uĺtima
      si_i=si[i];
      printf("i=%d m[si_i+maior=%d].d=m[si_ult + i=%d].d=%ld\n",i,si_i+maior, si_ult+maior,m[si_ult+maior].d);
   }
}


/*//Calcula distância da nova OTU (menor,maior) para as demais e atualiza vetor D
__global__ void dPDOTUs(int menor, int maior, int size, long int mij, m_dist *m, long int *D, int *SI){
   int stride = blockDim.x * gridDim.x;
   int index = threadIdx.x + blockIdx.x * blockDim.x;

   int l_menor=0, l_maior=0, l_k=0;
   l_menor=SI[menor];
   l_maior=SI[maior];

   for(int k=index;k<size;k+=stride){
      if(k<menor){
         D[k] = D[k] - m[l_menor+k].d;
         D[menor] = D[menor] - m[l_menor+k];
         m[l_menor+k] = (long int)(m[l_menor+k]+m[l_maior+k]-mij)*0.5;
         D[menor] = D[menor] + m[l_menor+k];
         D[k] = D[k] + m[l_menor+k];
      }// col= 0 até menor
//      __syncthreads();       
      

      if((k>menor+1) && (k<maior)){
         l_k = SI[k];
         D[k] = D[k] - m[l_k+menor];
         D[menor] = D[menor] - m[l_k+menor];
         m[l_k+menor] = (long int)(m[l_k+menor]+m[l_maior+k]-mij)*0.5;
         D[menor] = D[menor] + m[l_k+menor];
         D[k] = D[k] + m[l_k+menor];
      }// otu k = menor+1 até maior
//      __syncthreads(); 
      
      if(k>maior){
         l_k = SI[k];
         D[k] = D[k] - m[l_k+menor];
         D[menor] = D[menor] - m[l_k+menor];
         m[l_k+menor] = (long int)(m[l_k+menor]+m[l_k+maior]-mij)*0.5;
         D[menor] = D[menor] + m[l_k+menor];
         D[k] = D[k] + m[l_k+menor];
      }// lin = maior+1 até size
//      __syncthreads();
   }
   __syncthreads();
}//dPDOTUS
*/


//Calcula distância da nova OTU (menor,maior) para as demais e atualiza vetor D
__global__ void dPDOTUs_1(int menor, int maior, int t_vet, int size, long int mij, m_dist *m, long int *D, int l_menor, int l_maior){
   int stride = blockDim.x * gridDim.x;
   int index = threadIdx.x + blockIdx.x * blockDim.x;

   int l_k=0, j_k=0;

   for(int k=index;k<t_vet;k+=stride){
      l_k=m[k].i;
      j_k=m[k].j;
      if((m[k].i==menor)&&(m[k].j<menor)){
         D[j_k] = D[j_k] - m[l_menor+j_k].d;
         D[menor] = D[menor] - m[l_menor+j_k].d;
         m[l_menor+j_k].d = (long int)(m[l_menor+j_k].d+m[l_maior+j_k].d-mij)*0.5;
         D[menor] = D[menor] + m[l_menor+j_k].d;
         D[j_k] = D[j_k] + m[l_menor+j_k].d;
      }// colunas(OTU) = 0 até menor

      if ( ((m[k].i>menor+1) && (m[k].i<maior)) && (m[k].j==menor) ){
         D[l_k] = D[l_k] - m[l_k+menor].d;
         D[menor] = D[menor] - m[l_k+menor].d;
         m[l_k+menor].d = (long int)(m[l_k+menor].d+m[l_maior+l_k].d-mij)*0.5;
         D[menor] = D[menor] + m[l_k+menor].d;
         D[l_k] = D[l_k] + m[l_k+menor].d;
      }// otu k = menor+1 até maior
      
      if ( ((m[k].i>maior)&&(m[k].i<size)) && (m[k].j==menor) ){
         D[l_k] = D[l_k] - m[l_k+menor].d;
         D[menor] = D[menor] - m[l_k+menor].d;
         m[l_k+menor].d = (long int)(m[l_k+menor].d+m[l_k+maior].d-mij)*0.5;
         D[menor] = D[menor] + m[l_k+menor].d;
         D[k] = D[k] + m[l_k+menor].d;
      }// otu = maior+1 até size
   }
   __syncthreads();
}//distanciaParaDemaisOTUS-dPDOTUs_1
//**************************************
// Driver Code
//**************************************
int main(int argc, char *argv[])
{
   FILE* ptr = fopen(argv[1], "r");	//arquivo contendo matriz de distancia triangular
   double f = strtod(argv[2],NULL);	//f=p/q percentual de pares a serem selecionados de uma vez
   int BlockSize = atoi(argv[3]);
   int device=-1;

   clock_t t_start, t_end, t_1,t_2, t_i,t_f;
   double dt;

   cudaGetDevice(&device);

   if(ptr == NULL) {
     printf("no such file.\n");
     return 0;
   }

   int  size, size_O;		// No. OTU's
   int t0=fscanf(ptr,"%d", &size);
   size_O=size;
   float s = 0.0;
   s = (size/2);
   s = (size+1)*s;
   int t_vet = s-size;	//tamanho da matriz triangular
   
   int i,j,k;
   int menor,maior;
   long int mij, Di, Dj;

   int limite=(size*f);

   struct otuName otus[size];   // armazena nome dos nós-ramos da árvore
   p_node_type *pares = (p_node_type*)malloc(limite*sizeof(p_node_type));   // Lista COM PARES DE menor distância
   int *indices=(int *) malloc(size*sizeof(int));
   int v_somai[size];
   int soma_i=0;
   m_dist *m;    // MATRIZ DE DISTANCIAS
   long int *D;    // SOMA DAS DISTANCIAS DE UMA OTU(linha)
   p_node_type *didj; // vetor D[i]+D[j] para cálculo da matriz Q
   p_node_type *vetQ; //MatrizQ-vetor de ordenação
   p_node_type *d_vet; //Auxiliar para ordenação


   int pw=power2(t_vet);
   pw=pw*2;


   printf("fnj_gpu_%d SIZE=%d P/Q=%f t_vet=%d pw=%d\n",BlockSize, size,f,t_vet,pw);
   
   char str[24];
   int GridSize = (t_vet+BlockSize-1)/BlockSize;

//   printf("ALOCAÇÃO LINHAS DAS MATRIZES E INICIALIZAÇÃO\n");
   cudaMallocManaged(&D, size*sizeof(long int));          //matriz de distância - unified memory
   GridSize = (size+BlockSize-1)/BlockSize;
   init_D<<<GridSize,BlockSize>>>(size,D);     //INICIALIZA D
   cudaDeviceSynchronize();
   
   for(int i=0;i<size;i++){
      sprintf(str,"%d",i);
      otus[i].name = (char *) malloc(12*sizeof(char));
      sprintf(otus[i].name,"%d",i);
      otus[i].ordem=i;
      indices[i]=0;
//      D[i]=0;
   }

   cudaMallocManaged(&m, t_vet*sizeof(m_dist));        //matriz de distância - unified memory
   GridSize = (t_vet+BlockSize-1)/BlockSize;
   init_m<<<GridSize,BlockSize>>>(t_vet,m);     //INICIALIZA M
   cudaDeviceSynchronize();
   
   //for(int ix=0;ix<10;ix++){
   //   printf("m[%d]=%ld ",ix, m[ix]);
   //   printf("D[%d]=%ld\n",ix, D[ix]);
   //}

   //printf("Leitura da Matriz Triangular:\n");

   i=0;k=0;
   int v=0,somai=0;
   long int temp;
   while (i<size){		         // LENDO DISTANCIAS DO ARQUIVO..
     j=0;
     while((j<=i)&&(v<t_vet)){
       //printf("i=%d j=%d\n, ", i, j);
       if( (fscanf(ptr, "%ld ", &temp))==0) printf("Read Error\n");   // ...inserindo na matriz
       if((j<i)&&(i<size)) {
          m[v].d=temp;
          m[v].j=j;
          m[v].i=i;
          D[i]=D[i]+m[v].d;
          D[j]=D[j]+m[v].d;
//          if(v==9) printf("m[v=%d].i =%d .j=%d\n",v,m[v].i, m[v].j);
//          if(i==5) printf("m[v=%d].i=%d .j=%d\n",v,m[v].i, m[v].j);
          v++;}           //Soma das Linhas
       j++;
     }
     v_somai[k++]=soma_i;
     soma_i+=i;
     i++;
   }
   fclose(ptr);
   printf("\nFIM Leitura da Matriz Triangular: i=%d j=%d\n", i, j);
//   printf("m[v_somai[159]]+9].i=%d m[v_somai[159]]+9].j=%d \n",m[v_somai[159]+9].i,m[v_somai[159]+9].j);
   //Criando vetores na GPU
   long int *d_D;    //vetor D
   p_node_type *d_d; //vetor didj
   int *d_v; //auxiliar v_somai[size]
   
   cudaMallocManaged(&vetQ, pw*sizeof(p_node_type));  //matriz Q - unified memory
   GridSize = (pw+BlockSize-1)/BlockSize;
   init_Q<<<GridSize,BlockSize>>>(pw,vetQ);    //INICIALIZA VETQ
   cudaDeviceSynchronize();   
   
   checkCuda( cudaMalloc( (void**)&d_d,  t_vet * sizeof(p_node_type) ) ); 
   checkCuda( cudaMalloc( (void**)&d_vet,pw * sizeof(p_node_type) ) );
   checkCuda( cudaMalloc( (void**)&d_v,  size * sizeof(int) ) );  
   checkCuda( cudaMalloc( (void**)&d_D,  size * sizeof(long int) ));
   
   k=0;
   int soma_maior, soma_menor;
   int iteracao=1;

   didj = (p_node_type*) malloc( t_vet * sizeof(p_node_type));  // soma D[i]+D[j]

   limite = (size*f);
   if (limite==0) limite=1;
   k=0;
   t_1 = clock();
   t_start = clock();
   
   while(size > THRESHOLD){
      t_i=clock();
      printf("\n________________________________________________Size=%d\n",size);
         
      limite = (size*f);
      if (limite==0) limite=1;
      t_1=clock();
      soma_i=0;k=0;
      #pragma acc data copy(D)
      #pragma acc kernel
      {
         #pragma acc loop
         for(i=0;i<size;i++){
	    indices[i]=0;
            #pragma acc loop
            for(j=0;j<i;j++){
               didj[soma_i+j].d=D[i]+D[j];
               didj[soma_i+j].i=i;
               didj[soma_i+j].j=j;
               didj[soma_i+j].pv=soma_i+j;
            }
            soma_i+=i;
         }
      }
      #pragma acc update(didj,v_somai)
      //sincronizar
      
      t_2=clock();
//      for(int k1=0;k1<size;k1++)
//         printf("v_somai[%d]=%d\n",k1,v_somai[k1]);

//      int ci=1, cj=0;
//      if(iteracao<=2)
//      for(int k1=0;k1<t_vet;k1++){
//         printf("i=%d j=%d m[%d]=%ld \n",ci,cj,k1,m[k1].d);
//         cj++;
//         if(cj==ci){ ci++;cj=0;}
//      }

//      if(iteracao<=2)
//      for(int k1=0;k1<size;k1++)
//         printf("D[%d]=%ld \n",k1,D[k1]);

      pw=power2(t_vet);
      pw=pw*2;

/*      k=0;j=0;
      int x=0;
      for(i=0;i<45;i++){
         somai = didj[k].pv;
         printf("(size-2) * m[i=%d] - didj = %d * %ld - %ld = %ld\n",i, (size-2), m[i].d, didj[somai].d,( ( (size-2)*m[i].d) - didj[i].d) );
         if(j==k) {
            j=0;
         }
         else j++;
         k++;
      }
*/
      //printf("// copiando os vetores 'm' e 'D' para a GPU\n");
      checkCuda( cudaMemcpy( d_d,   didj, t_vet * sizeof(p_node_type), cudaMemcpyHostToDevice ) );
      cudaMemPrefetchAsync(m, t_vet*sizeof(long int), device, NULL);	
      cudaMemPrefetchAsync(vetQ, pw*sizeof(p_node_type), device, NULL);	
      
      cudaEvent_t startGPU, stopGPU;
      cudaEventCreate(&startGPU);
      cudaEventCreate(&stopGPU);
      float millisecondsGPU = 0;

      //printf("//Cálculo Matriz Q G=%d B=%d t_vet=%d size=%d pw=%d\n",GridSize, BlockSize,t_vet,size,pw);
      cudaEventRecord(startGPU);
      GridSize = (pw+BlockSize-1)/BlockSize;
      matrizQ<<< GridSize, BlockSize >>>( size, t_vet, m, vetQ, d_d, pw );
      cudaDeviceSynchronize();
      cudaEventRecord(stopGPU);
      cudaEventSynchronize(stopGPU);
      cudaEventElapsedTime(&millisecondsGPU, startGPU, stopGPU); 
      t_2=clock();
      printf("Matriz Q: %f s   ",((double)(t_2-t_1))/(double)CLOCKS_PER_SEC);
      printf("GPU Time: %f s\n", millisecondsGPU/1000 );
     
      cudaMemPrefetchAsync(vetQ, pw*sizeof(p_node_type), device, NULL);

//      if (iteracao==1) {
//        printf("\n\n t_vet=%d pw=%d\n",t_vet,pw);
//        for(i=0;i<45;i++)           
//              printf("i=%d (%d,%d)=%ld \n",i,vetQ[i].i,vetQ[i].j,vetQ[i].d);
//          printf("[%d]=%ld \n",i,vetQ[i].d); 
//      }
//      printf("\n\n");

      //printf("Ordenação na GPU com bitonicsort %d \n",t_vet);
      cudaEventCreate(&startGPU);
      cudaEventCreate(&stopGPU);
      millisecondsGPU = 0;
      t_1=clock();
      // Copiando o vetor para GPU
//      checkCuda( cudaMemcpy(d_vet, vetQ, t_vet * sizeof(p_node_type), cudaMemcpyHostToDevice) );
      //cudaMemPrefetchAsync(vetQ, t_vet*sizeof(p_node_type), device, NULL);

      cudaEventRecord(startGPU); 
      GridSize = (pw+BlockSize-1)/BlockSize;
      for (int k1=2; k1<=pw; k1<<=1)
      {
          for (int j1 = k1>>1 ; j1 > 0; j1 = j1 >> 1)
          {
              bitonicSortGPU << <GridSize, BlockSize >> > (vetQ, j1, k1);
          }
      }
      cudaDeviceSynchronize();
      //Copiando vetor ordenado
//      checkCuda( cudaMemcpy(vetQ, d_vet, t_vet * sizeof(p_node_type), cudaMemcpyDeviceToHost) );
      cudaEventRecord(stopGPU);
      cudaEventSynchronize(stopGPU);
      cudaEventElapsedTime(&millisecondsGPU, startGPU, stopGPU);
      t_2=clock();

      printf("Ordenação %f s   ",((double)(t_2-t_1))/(double)CLOCKS_PER_SEC);
      printf("GPU Time: %f s \n", millisecondsGPU/1000 );

//        printf("\n\n t_vet=%d pw=%d\n",t_vet,pw);
//        for(i=0;i<10;i++)
//              printf("i=%d (%d,%d)=%ld\n",i,vetQ[i].i,vetQ[i].j,vetQ[i].d);


//      if (iteracao==1)
//         for(i=0;i<5;i++)
//            printf("vet[%d].d=%ld vet[].i=%d vet[].j=%d vet[].pv=%d\n",i,vetQ[i].d,vetQ[i].i,vetQ[i].j,vetQ[i].pv);

      //printf("//selecionando os menores %d do heap z\n",limite);
      t_1=clock();
      if(limite == 1) 
         pares[0]=vetQ[1];
      else{
         int l;
         j=0;l=0;
         while((l<limite)&&(j<t_vet)){                //seleciona os menores o inicio do vetor
            if(indices[vetQ[j].i]==0 && indices[vetQ[j].j]==0){
//              printf("Selecionou %d (%d,%d)=%ld \n",j,vetQ[j].i,vetQ[j].j,vetQ[j].d);//pv=%d v_somia[j]=%d  ,vetQ[j].pv,v_somai[vetQ[j].j]
              pares[l++]=vetQ[j];
//              printf("Selecionou i=%d (%d,%d)=%ld \n",l-1,pares[l-1].i,pares[l-1].j,pares[l-1].d);
              //printf("Selecionou i=%d (%d,%d)=%ld pv=%d v_somia[i]=%d\n",l-1,pares[l-1].i,pares[l-1].j,pares[l-1].d,pares[l-1].pv,v_somai[pares[l-1].i]);
              indices[vetQ[j].i]=1;
              indices[vetQ[j].j]=1; }
            j++;
         }//while()
      }//if(limite == 1)
      t_2=clock();
      printf("Sel.Pares %f s\n",((double)(t_2-t_1))/(double)CLOCKS_PER_SEC);
      //for(int x1=0;x1<10;x1++)
      //   printf("%d pares[].i=%d .j=%d .pv=%d .d=%ld\n",x1, pares[x1].i,pares[x1].j,pares[x1].pv,pares[x1].d); 

      t_1=clock();
//      checkCuda(cudaMemcpy( d_v, v_somai,   size*sizeof(int), cudaMemcpyHostToDevice ));
      printf("Juntando %d Pares de OTU's\n",limite);
      for(int c=0; c<limite; c++) {
         //if(size >2 ){
            double d_i_novo,d_j_novo=0.0;
            double t1, t2;
            int pv;
            //printf("// determina menor coordenada (%d,%d)- nomenclaturas %d \n",pares[c].i,pares[c].j,c);
            if (pares[c].i > pares[c].j) { // determina menor coordenada - nomenclaturas
               mij = m[pares[c].pv].d;
               pv=pares[c].pv;
               Di=D[pares[c].j]; 
               Dj=D[pares[c].i];
               menor=pares[c].j;
               maior=pares[c].i; 
            }
            else{
               mij = m[pares[c].pv].d;
               pv=pares[c].pv;
               Di = D[pares[c].i]; 
               Dj = D[pares[c].j];
               menor=pares[c].i;
               maior=pares[c].j;
            }//if (pares[c].i > pares[c].j)
//            printf("mij=%ld didj=%ld menor=%d maior=%d\n",mij,Di+Dj, menor, maior);
            
//            printf("ENTRE c=%d pares[c].pv = %d mij=%ld m[pares[c].pv]=%ld\n",c, pares[c].pv, mij, m[pares[c].pv].d);

            //printf("//nomeando e determinando a ordem da nova OTU\n");
            char aux[(strlen(otus[menor].name)+strlen(otus[maior].name)+4)];
            if (otus[menor].ordem < otus[maior].ordem)
               sprintf(aux, "(%s,%s)",otus[menor].name,otus[maior].name);
            else{
               sprintf(aux,"(%s,%s)",otus[maior].name,otus[menor].name);
               otus[menor].ordem=otus[maior].ordem;
            }
            free(otus[menor].name);
            int y;
            y = 8 - (strlen(aux)%8) + 1;
            if(!(otus[menor].name=(char*) malloc((strlen(aux)+y)*sizeof(char))))
               printf("Não alocou*************\n");
            strcpy(otus[menor].name, aux);
            
            //printf("Calculando distancias entre a nova OTU\n");
            t1=(size-c - 2) * 2;
            t2= (Di - Dj);
            t2 = t2/t1;
            d_i_novo = (0.5*mij);
            d_i_novo = d_i_novo + t2 ;
            d_j_novo = mij - d_i_novo;
            // printf("calculando distancia do novo (%d,%d) aos demais %d - calculando a linha-coluna\n",menor,maior,size);
            //printf("D[menor=%d]=%ld\n D[maior=%d]=%ld\n",menor,D[menor],maior,D[maior]);
            soma_menor=v_somai[menor];
            soma_maior=v_somai[maior];
            somai=soma_menor;
            i=menor;j=0;
//            D[menor]=0;    // Zera D[menor] para recalcular

//   checkCuda(cudaMemcpy( d_D, D,  size * sizeof(long int), cudaMemcpyHostToDevice ));
//   printf("menor=%d maior=%d somai=%d size=%d soma_maior=%d\n",menor,maior,somai, size, soma_maior);
   // Decrementa m[menor] da soma de todos D                                        
//   distancesPTO_1<<< GridSize, BlockSize >>>(menor, soma_menor, size, m, d_D, d_v);
           printf("\njuntando c=%d otus[menor=%d].name=%s mij=%ld\n",c,menor,otus[menor].name,mij);
           int l_menor=v_somai[menor],l_maior=v_somai[maior],l_k, j_k;

           for(int k1=0;k1<t_vet;k1++){
              l_k=m[k1].i;
              j_k=m[k1].j;
              int teste=0;
//              teste = (int)(l_k*(l_k-1))/2;
              //if(k1!=teste) printf("ERRO ERRO !!!! ");
              //((c==13)||(c==28))&&
//              if((((teste+j_k)!=k1)))printf("A0-m[l_k=%d k1=%d].i=%d m[l_k=%d k1=%d].j=%d j_k=%d menor=%d teste=%d\n",l_k, k1, m[k1].i, l_k, k1,m[k1].j, j_k, menor,teste);
      
//              if((m[k1].i==menor)&&(m[k1].j<menor)){
                 if( ((m[k1].i==menor==3)&&(m[k1].j==349)) || ((m[k1].i==menor==3)&&(m[k1].j==2)) )
                    printf("A1-m[%d][%d]= (m[%d][%d]=%ld + m[%d][%d]=%ld - mij=%ld)*0.5=%ld\n",m[l_menor+j_k].i,m[l_menor+j_k].j, m[l_menor+j_k].i, m[l_menor+j_k].j, m[l_menor+j_k].d, m[l_maior+j_k].i, m[l_maior+j_k].j, m[l_maior+j_k].d, mij, (long int)((m[l_menor+j_k].d+m[l_maior+j_k].d-mij)*0.5));
//              }// colunas(OTU) = 0 até menor

//              if ( ((m[k1].i>menor+1) && (m[k1].i<maior)) && (m[k1].j==menor) ){
                 if ( ((menor==2)&&(k1==3)) || ((menor==3)&&(k1==349)) )
                    printf("A2-m[%d][%d]= (m[%d][%d]=%ld + m[%d][%d]=%ld - mij=%ld)*0.5=%ld\n",m[l_k+menor].i,m[l_k+menor].j,m[l_k+menor].i, m[l_k+menor].j, m[l_k+menor].d, m[l_maior+l_k].i, m[l_maior+l_k].j, m[l_maior+l_k].d, mij, (long int)((m[l_k+menor].d+m[l_maior+l_k].d-mij)*0.5));
//              }// otu k = menor+1 até maior
      
//              if ( ((m[k1].i>maior)&&(m[k].i<size)) && (m[k1].j==menor) ){
                 if ( ((menor==2)&&(k==3)) || ((menor==3)&&(k==349)) )
                    printf("A3-m[%d][%d]= (m[%d][%d]=%ld + m[%d][%d]=%ld - mij=%ld)*0.5=%ld/n" , m[l_k+menor].i, m[l_k+menor].j, m[l_k+menor].i, m[l_k+menor].j, m[l_k+menor].d, m[l_k+maior].i, m[l_k+maior].j, m[l_k+maior].d, mij, (long int)((m[l_k+menor].d+m[l_k+maior].d-mij)*0.5));
//              }// otu = maior+1 até size
           }

           //printf("Distâncias (menor=%d,maior=%d\n",menor,maior);
           GridSize = (size+BlockSize-1)/BlockSize;
           dPDOTUs_1<<< GridSize, BlockSize >>>(menor, maior, t_vet, size, mij, m, D, v_somai[menor], v_somai[maior]);
           cudaDeviceSynchronize();
//           checkCuda(cudaMemcpy( m, m,   t_vet*sizeof(long int), cudaMemcpyDeviceToHost ));
           cudaDeviceSynchronize();
           if ((menor==2))
              printf("D-m[3][2]= %ld\n", m[v_somai[3]+2].d);
           if ((menor==3))
              printf("D-m[349][3]= %ld\n", m[v_somai[349]+3].d);

//           printf("D-m[v_somai[159]+9].i=%d m[v_somai[159]+9].j=%d \n",m[v_somai[159]+9].i,m[v_somai[159]+9].j);

//           printf("depois c=%d pares[c].pv = %d mij=%ld m[60313]=%ld\n",c, pares[c].pv, mij, m[60313]);
//   checkCuda(cudaMemcpy( D, d_D, size * sizeof(long int), cudaMemcpyDeviceToHost ));
//   checkCuda(cudaMemcpy( m , m, t_vet * sizeof(long int), cudaMemcpyDeviceToHost ));
   


//            if((iteracao==1)&&(c==0))
//               for(int g=0;g<10;g++)
//                  printf("D-%d - D[g]=%ld \n",g,D[g]);
                  //printf("%d - m[soma_menor+g=%d]=%ld \n",g,soma_menor+g,m[soma_menor+g]);

//            printf("D - D[0]=%ld m[soma_menor=%d+0]=%ld D[menor=%d]=%ld\n",D[0],soma_menor,m[soma_menor+0],menor,D[menor]);

//            printf("Eliminando linha-coluna com cópia (%d,%d) c=%d.\n",menor,maior,c);
           if (sizeof(otus[maior].name) < sizeof(otus[size-1])){
              // free(otus[maior].name);
              otus[maior].name = (char *)malloc(sizeof(otus[size-1].name)*sizeof(char));
           }
           otus[maior].name = otus[size-1].name;
           otus[maior].ordem= otus[size-1].ordem;

           int soma_ult=v_somai[size-1];
           i=maior;somai=soma_maior;                 //vai percorrer a matriz a partir da linha maior
           D[maior]=D[size-1];                       //copia soma da ultima linha para linha maior
           k=maior+1;                                // k indica a coluna a frente de maior na ultima linha
           j=0;                                      // j percorre todas as colunas quando linha for maior..
            
//            cudaMemPrefetchAsync(m, t_vet*sizeof(long int), device, NULL);	
            //if((menor==2) || (menor==3)) 

//            printf("Copia linha (maior=%d[%d],última=%d[%d]\n",maior,v_somai[maior],size-1,v_somai[size-1]);

//   for(int j1=0;j1<maior;j1++)			//percorre a linha maior até coluna maior copiando da última
//      printf("m[maior=%d][j=%d]=m[ultima=%d][j=%d]=%ld\n",maior, j1, size-1,j1, m[v_somai[size-1]+j1]);
	    
   
//   for(int i1=maior+1;i1<size-1;i1++){              //percorre linhas de maior até a penúltima copiando da uĺtima
//      somai=v_somai[i1];
//      printf("m[i=%d][maior=%d]=m[ultima=%d][i=%d]=%ld\n",i1, maior, size-1, i1, m[v_somai[size-1] + i1]);
//   }
            
/*            int si_ult=m[size-1].i;
            for(int k1=0;k1<t_vet;k1++){
             if(k1==9){
               int l_k=m[k1].i;
               int j_k=m[k1].j;
       
               if ((l_k==l_maior)&&(j_k<maior))			//se linha igual a maior e coluna vai até maior-1
                  printf("C-m[k=%d].d=m[si_ult=%d+j_k=%d].d=%ld\n",k1,si_ult,j_k,m[si_ult+j_k].d); // copia da ultima linha até coluna maior para linha maior na coluna
      
               if(((l_k>l_maior)&&(l_k<size))&&(j_k==maior))	//se linha maior que maior e coluna igual maior  
                  printf("C-m[k=%d].d=m[si_ult=%d+j_k=%d].d=%ld\n",k1,si_ult,j_k,m[si_ult+j_k].d); // copia da ultima linha coluna maior para linha k coluna maior
             }
            }
*/
            int teste1=m[10].i;
            GridSize = (size+BlockSize-1)/BlockSize;
            copyLine_1<<< GridSize, BlockSize >>>(maior, size-1, t_vet, size, m, d_v);
            cudaDeviceSynchronize();
            if(m[10].i!=teste1) printf("TESTE 1 - ERRADO\n");
           checkCuda(cudaMemcpy( m, m,   t_vet*sizeof(long int), cudaMemcpyDeviceToHost ));
           cudaDeviceSynchronize();
           v_somai[size-1]=v_somai[maior];
//           if ((menor==2))
//              printf("L-m[3][2] = %ld\nL-m[349][3] = %ld\n", m[v_somai[3]+2].d, m[v_somai[349]+3].d);
//           if ((menor==3))
//              printf("L-m[3][2] = %ld\nL-m[349][3] = %ld\n", m[v_somai[3]+2].d, m[v_somai[349]+3].d);
//              printf("L-m[v_somai[159]].i=%d m[v_somai[159]].j=%d = %ld\n", m[159].i, m[159].j, m[m[159].i+9].d);
//            printf("v_somai[ultima=%d]=%d v_somai[maior=%d]=%d\n",size-1, v_somai[size-1],maior,v_somai[maior]);
            
//printf("v_somai[ultima=%d]=%d v_somai[maior=%d]=%d\n",size-1, v_somai[size-1],maior,v_somai[maior]);

//             for(int k1=0;k1<10;k1++)
//                printf("m[maior,%d]=%ld = m[size-1,%d]=%ld\n",k1,m[v_somai[maior]+k1],k1,m[v_somai[size-1]+k1]);


//            v_somai[maior]=v_somai[size-1];
//            printf("Corrigindo coordenadas. c=%d limite=%d %d \n",c,limite,pares[c+1].i);
//            #pragma acc parallel
            {
//               #pragma acc loop
               for(int e=c+1; e<limite; e++){            //percorre pares selecionados ainda no vetor   
                  if(pares[e].i==size-1){                // e aqueles que estavam na última linha
                     printf("CORRIGINDO PARES[E=%d].pv=%d  ",e,pares[e].pv);
                     if(pares[e].j > maior) {            // se estavam na coluna a frente de maior
                        pares[e].i=pares[e].j;           // foram copiados para a linha j e ...
                        pares[e].j=maior;                // ..coluna maior
                     }
                     else pares[e].i=maior;              // em colunas até maior apenas trocam de linha para maior

                     pares[e].pv=v_somai[pares[e].i] + pares[e].j;                     
//                     printf(" => pv= %d\n",pares[e].pv);

                     printf("(%d,%d).pv=%d\n",pares[e].i, pares[e].j, pares[e].pv);
                  }//if(pares[e].i==size-1)
               }//for(int e=c+1; e<limite; e++)
            }
            size--; 
         //}//if(size >2 ){
      }//for(int c=0; c<limite; c++)
      t_2=clock();
      printf("Jun.Pares %f s\n",((double)(t_2-t_1))/(double)CLOCKS_PER_SEC);

//        printf("\n\n t_vet=%d pw=%d\n",t_vet,pw);
//        for(i=0;i<pw;i++)
//           if(((vetQ[i].i==-1) || vetQ[i].j==-1) && (i<=t_vet))
//              printf("i=%d (%d,%d)\n",i,vetQ[i].i,vetQ[i].j);

      s = (size/2);
      s = (size+1)*s;
      t_vet = s-size;
      t_f=clock();
//      printf("size=%d t_vet=%d\n",size,t_vet);
      printf("Iteração  %f s\n",((double)(t_f-t_i))/(double)CLOCKS_PER_SEC);
      iteracao++;
   }//while(size > THRESHOLD)
   //________________________________________________________________________________________
   
   printf("\n\nPARALELO %f s\n",((double)(clock()-t_start))/(double)CLOCKS_PER_SEC);
   clock_t t1=clock();

   free(pares);
   free(indices);
   free(didj);

   int t=s-size;
   printf("SEQUENCIAL SIZE=%d f=%1.4f ",size,f);
   printf("Tam Matriz t=%d Tam. Heap = %d\n",t, limite);

//****************************************************************************
//****************************************************************************
   int inter=0;
   p_node_type h[limite];  // HEAP
   clock_t t_11, t_21, t2;
   //double temp1, temp2;
   inter = 0;
   while(size>2){
      t_i=clock();
      printf("\n________________________________________________Size=%d\n",size);
      s = 0.0;
      s = size/2.0;
      s = s*(size+1);					//******DETERMINANDO TAMANHO MATRIZ TRIANGULAR
      limite = size*f;					// número de pares a serem selecionados
      if (limite==0) limite=1;
      int t = s-size;
      i=0;v=0,somai=0;

      //printf("Iteração=%d, SIZE=%d t_vet=%d ===============\n",inter,size,t);
      while (i<size){
         j=0;
         while((j<i)&&(v<t)){
            D[i]=D[i]+m[v].d;
            D[j]=D[j]+m[v].d;
            v++;
            j++;
         }
         somai+=i;
         v_somai[i+1]=somai;
         //printf("v_somai[%d]=%d\n",i, v_somai[i]);         
         i++;
      }

      //printf("Calculando Matriz Q %d\n",size);
      t_1=clock();
      matrizQ(vetQ, size, m, D);			// CALCULANDO MATRIZ Q
      t_2=clock();
      printf("Matriz Q: %f s\n",(t2-t1)/(double)CLOCKS_PER_SEC);
      
      t_1=clock();
      minimoModificado(h,size, t, vetQ, f, inter);      // Seleciona Menores
      t_2=clock();
      //if (inter==0) printf("Selecionar menores     t=%f \n",(t_2-t_1)/(double)CLOCKS_PER_SEC);

      //printf("Juntando Pares de OTU's limite=%d \n",limite);
      printf("Juntando %d Pares de OTU's\n",limite);
      t_1=clock();
      for(int c=0; c<limite; c++) {
         t_11=clock();
         if (size > 2) {
            //double d_i_novo,d_j_novo=0.0;
            if (h[c].i > h[c].j) {
               mij=m[h[c].pv].d; 
               Di=D[h[c].j]; 
               Dj=D[h[c].i];
               menor=h[c].j;
               maior=h[c].i;
            } 
            else {
               mij = m[h[c].pv].d;
               Di = D[h[c].i]; 
               Dj = D[h[c].j];
               menor=h[c].i;
               maior=h[c].j; 
            }
            char aux[(strlen(otus[menor].name)+strlen(otus[maior].name)+4)];
            if(otus[menor].ordem  < otus[maior].ordem)
               sprintf(aux,"(%s,%s)",otus[menor].name,otus[maior].name);
            else{
               sprintf(aux,"(%s,%s)",otus[maior].name,otus[menor].name);
               otus[menor].ordem=otus[maior].ordem;
            }
            free(otus[menor].name);
            int y;
            y = 8 - (strlen(aux)%8) + 1;
            if(!(otus[menor].name=(char*) malloc((strlen(aux)+y)*sizeof(char))))
                printf("Não alocou*************\n");
            strcpy(otus[menor].name, aux);
            printf("juntando c=%d otus[menor=%d].name=%s\n",c,menor,otus[menor].name);
//            printf("otus[menor=%d].name=%s\n",menor,otus[menor].name);

            //CÁLCULO DA DISTÂNCIA DOS RAMOS DO PAR AGRUPADO--------------
            //temp1=(size - 2) * 2;
            //temp2= (Di - Dj);
            //temp2 = t2/t1;
            //d_i_novo = (0.5*mij);
            //d_i_novo = d_i_novo + temp2 ;
            //d_j_novo = mij - d_i_novo;
            
            //printf("calculando distancia do novo aos demais ---- calculando a linha-coluna menor=%d maior=%d sizze=%d\n",menor,maior,size);
            int somai_menor=v_somai[menor];
            int somai_maior=v_somai[maior];
                                    
            t_11=clock();
            //printf("0-menor=%d\n",menor);
            for(int k=0;k<menor;k++){
               D[menor] = D[menor] - (m[somai_menor+k].d);
               D[k] = D[k] - (m[somai_menor+k].d+m[somai_maior+k].d);
               m[somai_menor+k].d = (int)((m[somai_menor+k].d+m[somai_maior+k].d)-mij)*0.5;
               D[k] = D[k] + m[somai_menor+k].d;
               D[menor] = D[menor] + m[somai_menor+k].d;
            }
            
            //printf("menor+1 - maior=%d\n",maior);
            for(int k=menor+1;k<maior;k++){
               D[menor] = D[menor]-m[v_somai[k]+menor].d; //-----
               D[k] = D[k] - (m[v_somai[k]+menor].d+m[somai_maior+k].d);
               m[v_somai[k]+menor].d= (int)((m[v_somai[k]+menor].d+m[somai_maior+k].d)-mij)*0.5;
               D[k] = D[k] + m[v_somai[k]+menor].d;
               D[menor]=D[menor]+m[v_somai[k]+menor].d;//---
            }
            
            //printf("maior=%d-size=%d menor=%d v_somai[maior]=%d \n",maior,size,menor,v_somai[maior]);
            for(int k=maior;k<size;k++){
               //printf("maior=%d k=%d v_somai[k]=%d   \n",maior, k, v_somai[k]);
               D[menor] = D[menor]-m[v_somai[k]+menor].d; //printf("a,");//-----
               D[k] = D[k] - (m[v_somai[k]+menor].d+m[v_somai[k]+maior].d); //printf("b,");
               m[v_somai[k]+menor].d = (int)((m[v_somai[k]+menor].d+m[v_somai[k]+maior].d)-mij)*0.5; //printf("c,");
               D[k] = D[k] + m[v_somai[k]+menor].d;//printf("d,");
               D[menor]=D[menor]+m[v_somai[k]+menor].d;//printf("e.    ");//---
            }
            t_21=clock();
            //if (c==0 && inter==0) printf("\nDistância para demais   t=%f \n",(t_21-t_11)/(double)CLOCKS_PER_SEC);
            
            
            //printf("Eliminando linha-coluna com cópia. \n");
            t_11=clock();
            for(int k=0;k<size;k++) 
               if (k>maior) m[v_somai[k]+maior].d=m[v_somai[size-1]+k].d;
               else if(k==maior) m[somai_maior+k].d=0;
                    else m[somai_maior+k].d=m[v_somai[size-1]+k].d;
            t_21=clock();
            //if (c==0 && inter==0) printf("Eliminando linha-coluna t=%f \n",(t_21-t_11)/(double)CLOCKS_PER_SEC);
            
            free(otus[maior].name);
            otus[maior].name = otus[size-1].name;
            otus[maior].ordem=otus[size-1].ordem;
         

            //printf("Corrigindo coordenadas.\n");
            for(int e=c+1; e<limite; e++){
               if(h[e].i==size-1){
                  if (h[e].j > maior) {
                     h[e].i=h[e].j;
                     h[e].j=maior;
                  }
                  else h[e].i=maior;
               
                  //if (h[e].j==size-1)
                  //   h[e].j =maior;
                  h[e].pv=v_somai[h[e].i] + h[e].j;   
               }
            }   
            t_21=clock();
            //if (c==0 && inter==0) printf("Corrigindo Coordenadas t=%f \n",(t_21-t_11)/(double)CLOCKS_PER_SEC);
             
         }//if (size > 2) {
         size--;
      }//for(int c=0; c<limite; c++)
      t_2=clock();
      printf("Jun.Pares %f s\n",((double)(t_2-t_1))/(double)CLOCKS_PER_SEC);
      t_f=clock();
      printf("Iteração  %f s\n",((double)(t_f-t_i))/(double)CLOCKS_PER_SEC);
      inter++;
   }//while(size>2)

//****************************************************************************
//****************************************************************************

   printf("SEQUENCIAL %f s\n",((double)(clock()-t1))/(double)CLOCKS_PER_SEC);
   t_end = clock();
   dt = ((double)(t_end-t_start)) / (double)CLOCKS_PER_SEC;
   printf("\nTime =%f s\n",dt);      


//   printf("[0]\n %s\n",otus[0].name);
//   printf("[1]\n %s\n",otus[1].name);
//   if (otus[0].ordem<otus[1].ordem)
//      printf("(%s,%s)\n",otus[0].name,otus[1].name);
//   else
//      printf("(%s,%s)\n",otus[1].name,otus[0].name);   

   char aux[(strlen(otus[0].name)+strlen(otus[1].name)+4+1)];
   if(otus[0].ordem  < otus[1].ordem)
     sprintf(aux,"(%s,%s)",otus[0].name,otus[1].name);
   else
     sprintf(aux,"(%s,%s)",otus[1].name,otus[0].name);
   printf("%s\n",aux);
   
   cudaFree(d_d);
   cudaFree(d_v);
   cudaFree(d_vet);
   cudaFree(d_D);
   cudaFree(vetQ);
   cudaFree(m);
   cudaFree(D);
//   for(int i=0;i<size_O;i++)
//     free(otus[i].name);

   //cudaFree(d_m);
   //cudaFree(d_vet);
   return 0;
}//main()







