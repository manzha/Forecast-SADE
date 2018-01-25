#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <limits.h>

// Problem Properties
typedef struct {
	unsigned int dim;           // Solution size (alpha, beta, phi1, ..., phin)
	unsigned int populationSize;// Population size
	double t0;                  // Initial temperature.
	double tf;                  // Final temperature.
	double alpha;               // Alpha coefficient to update the temperature.
} sa_t;

//Timeserie struct
typedef struct {
	double *t; //Timeserie data
	int size; //Size of timeserie
} timeserie;

//Solution struct
typedef struct {
	double **solution; //Sol. population
	double *energy; //MAPE for each solution
	double *smape; //SMAPE for each solution
	double *error; //% error for each solution
	double *DA;  // DA for each solution (Directional Accuracy) as in Santamaria-bonfil article...
} sa_state;

//New solution. DE struct
typedef struct {
	double *solution; //Just one solution
	double energy; //MAPE
	double smape; //SMAPE
	double error; //Error
	double DA; //Hit Rate (DA)
} ts_sol;


void sa_init(sa_t *t, sa_state *current, sa_state *next, int degree, int pobSize);
void sa_genInitialSolution(sa_t *t, sa_state *state, timeserie *serie, int start, int end);
void sa_calculateEnergySol (sa_t *t, ts_sol *sol, timeserie *serie, int start, int range, double errorD);
void sa_calculateEnergySolTest (sa_t *t, ts_sol *sol, timeserie *serie, int N, double error);
void sa_copySolution(sa_t *t, sa_state *from, sa_state *to);
void sa_copyOneSolution(sa_t *t,sa_state *state, double* S,int idx);
void sa_printSolution(sa_t *t, sa_state *state);
void sa_printBestSolution(sa_t *t, sa_state *state, int *bestIdx,float seconds,double bfsenergy, double bfssmape, double bfsDA);
void sa_findBestSolution(sa_t *t, sa_state *state, int *bestIdx);
double sa_updateTemperature(sa_t *t, double temperature);
void de_newSample(sa_t *t, sa_state *state, int p0, double F, double CR, ts_sol* S);
void ts_read(timeserie *t,timeserie *t2);
void ts_print(sa_t *t, sa_state *state,timeserie *serie,int bestIdx, int start, float seconds);
double expSmooth(double *errorForecast, int n);
void de_newPob(sa_t *sa, ts_sol *S,timeserie *train,sa_state *next, sa_state *current,int iteracion, int range,double errorForecast,double F,double CR, double t, double bestEnergyinPob, int *iBetter,ts_sol *bestFoundSolution);
int dynEq(sa_state *next, int iBetter, int l,double Lk,int Laux, int iter,double Summary,double Esummary,double k1, double k2,double m, double epsilon, int* contadorG);

void sa_init(sa_t *t, sa_state *current, sa_state *next, int degree, int pobSize) {
	int i;
	t->t0 = 57118.5772913; //This was tuned previously
	//t->t0 = 1000;
	t->tf = 0.00001101567; //Tuned previously
	//t->tf=0.01;
	t->alpha = 0.95; // 95%
	t->dim = degree+2; //As in the article.
	t->populationSize = pobSize;
	printf("Population size: %lf\n",t->populationSize);
	printf("Initial temperature: %lf\n", t->t0);
	printf("Final temperature: %lf\n", t->tf);
	printf("Alpha: %lf\n", t->alpha);

	current->solution = malloc(t->populationSize * sizeof(double*));
	next->solution = malloc(t->populationSize * sizeof(double*));

	current->energy = malloc(t->populationSize * sizeof(double));
	next->energy = malloc(t->populationSize * sizeof(double));

	current->smape = malloc(t->populationSize * sizeof(double));
	next->smape = malloc(t->populationSize * sizeof(double));

	current->error = malloc(t->populationSize * sizeof(double));
	next->error = malloc(t->populationSize * sizeof(double));

	current->DA = malloc(t->populationSize * sizeof(double));
	next->DA = malloc(t->populationSize * sizeof(double));

	for(i = 0; i < t->populationSize ; i ++){
		current->solution[i] = malloc(t->dim * sizeof(double));
		next->solution[i] = malloc(t->dim * sizeof(double));
	}
}

void sa_genInitialSolution(sa_t *t, sa_state *state, timeserie *serie, int start, int end) {
	int N=serie->size;
	//start= N*(start/100);
	//end= N*(end/100);
	end=N;
	double MAD=0;
  double epsilon=10.8;
  double sumX=0, sumY=0,sumX2=0, sumXY=0;
  double A,B;
  int i,j;
	double* forecasted = malloc(sizeof(double)*N);
	//DESCARTA DATOS ANTIGUOS ?????
	double aux;
  for(i=start;i<end;i++){
		//printf("%lf ",serie->t[i]);
		if(serie->t[i]<=0){
			if(serie->t[i]==0){
				aux=0;
			}else{
				aux=log(serie->t[i]*-1);
			}
		}else{
			aux=log(serie->t[i]);
		}
    sumX+=(i+1);
    sumY+=aux;
    sumX2+=pow((i+1),2);
    sumXY+=((i+1)*aux);
		//printf("%lf %lf %lf %lf\n",sumX,sumY,sumX2,sumXY);
  }
	//printf("%lf %lf %lf %lf\n",sumX,sumY,sumX2,sumXY);
  //B = (sumXY - (sumX*sumY)/N)/(sumX2 + (sumX*sumX)/N);
  //A = (sumY - sumX*B)/N;
	B = (sumX*sumY-N*sumXY)/(sumX*sumX-N*sumX2);
//	B = 0.0011;
	A = (sumY-sumX*B)/N;
  A = exp(A);

	//Suavizamiento de datos originales --CORREGIR TODO ESTO
	//REVISAR TODO ESTO!!!
	FILE *solFile = fopen("test.serie","w");
  for(i=0 ; i < serie->size ; i++){
		//printf("%d ",i);
    forecasted[i] = A*exp((B*(i+1)));
    //forecasted[i] = exp(A+(B*(i+1)));
    MAD = fabs((serie->t[i]-forecasted[i])/serie->t[i]);
    //printf("%lf %lf\n",state->forecast[0].t[i],MAD);
		//printf("MAD: %lf\n",MAD);
		//TODO
		//CORREGIR!!!
    if(MAD>epsilon){
      if(serie->t[i] > forecasted[i]){
        serie->t[i]=forecasted[i]+epsilon;
      }else{
        serie->t[i]=forecasted[i]-epsilon;
      }
    }
		//fprintf(solFile,"%lf\n", serie->t[i]);
		fprintf(solFile,"%lf\n", forecasted[i]);
  }
	printf("A = %lf ; B= %lf\n",A,B);
	fclose(solFile);

	//Population
	//int K = 20;
	//int degree=4;
	//double population[K][t->dim];
	state->solution[0][0]=A;
	state->solution[0][1]=B;
	//population[0][0]=A;
	//population[0][1]=B;
	//float xRand = (float)(rand())/(float)(RAND_MAX/a);
	for(i = 2; i < t->dim ;i++){
		state->solution[0][i]=-1+(float)(rand())/(float)(RAND_MAX/2);
		printf("Phi: %lf",state->solution[0][i]);
	}
  // Mutate Sinit, K times to create Population 0 (P0)
	double range=0,delta=0;
	int l;
	double ai;
	for(i = 1 ; i < t->populationSize ; i++){
		for(j = 0; j < t->dim ; j++){
			//population[i][j]= population[0][j];
			state->solution[i][j]=state->solution[0][j];
			ai=rand()%2;
			if(ai==1){
				ai=((float)rand())/RAND_MAX;
				state->solution[i][j]*=ai;
			}else{
				ai=-((float)rand())/RAND_MAX;
				state->solution[i][j]*=ai;
			}
		}
	}
	printf("A = %lf\n B= %lf\n",A,B);
	free(forecasted);

}

void sa_calculateEnergySol (sa_t *t, ts_sol *S, timeserie *serie, int iteracion, int windowSize, double errorD) {
	int flag=0;
	int j,k,ti;
	int N=serie->size;
	int end;
	double error=0;
	double alpha,beta,fi=0;
	double* forecast = malloc(sizeof(double) * serie->size);
	iteracion = iteracion*(windowSize);
	end = iteracion+windowSize;
	if(end>serie->size){
		end=serie->size;
		windowSize=end-iteracion;
	}
	alpha=S->solution[0];
	beta=S->solution[1];
	fi=0;
	ti=iteracion;
	for(j=iteracion ; j< end ;j++){
		fi=0;
		for(k=2 ; k < t->dim; k++){
			fi+=(S->solution[k]*pow(ti,k-2));
		}
		forecast[j] = alpha*exp(beta*ti) + fi;
		ti=ti+1;
		//forecast[j]=forecast[j]+(errorD);
	}
	//Calculate for forecast i
	S->energy=0;
	S->smape=0;
	S->DA=0;
	for(j = iteracion ; j < end ; j++){
		S->energy += (fabs(serie->t[j] - forecast[j])/fabs(serie->t[j]));
		S->smape += (fabs(forecast[j] - serie->t[j])/((fabs(serie->t[j]) + fabs(forecast[j]))/2));
		error += (serie->t[j] - forecast[j]);
		if(j!=iteracion){
			if((serie->t[j] > serie->t[j-1] && forecast[j]>forecast[j-1])
				|| (serie->t[j] < serie->t[j-1] && forecast[j]<forecast[j-1])){
				S->DA +=1;
			}
		}
	}
	S->energy /= (windowSize*1.0);
	S->energy *=100;
	S->smape /= (windowSize*1.0);
	S->smape *=100;
	S->DA /= (windowSize*1.0);
	S->DA *=100;
	S->error = error/(windowSize*1.0);
	free(forecast);
}
void sa_calculateEnergySolTest (sa_t *t, ts_sol *S, timeserie *serie, int N, double error) {
	int j,k;
	double alpha,beta,fi=0,ti;
	double errorx=0;
	double* forecast = malloc(sizeof(double) * serie->size);
	double errorPtg=1;
	alpha=S->solution[0];
	beta=S->solution[1];
	fi=0;
	ti=N;
	//printf("Alpha %lf Beta %lf %lf %lf\n",alpha,beta,S->solution[2],S->solution[3]);
	for(j=0 ; j< serie->size;j++){
		fi=0;
		for(k=2 ; k < t->dim; k++){
			fi+=(S->solution[k]*pow(ti,k-2));
		}
		forecast[j] = alpha*exp(beta*ti) + fi;
		ti=ti+1;
		//forecast[j]=forecast[j]+(errorPtg*error); //AJUSTE De ERROR
		//errorPtg-=0.03;
	}
	//Calculate for forecast i
	S->energy=0;
	S->smape=0;
	S->DA=0;
	for(j = 0 ; j < serie->size ; j++){
		//printf("%lf %lf\n",serie->t[j],forecast[j]);
		S->energy += (fabs(serie->t[j] - forecast[j])/fabs(serie->t[j]));
		S->smape += (fabs(forecast[j] - serie->t[j])/((fabs(serie->t[j]) + fabs(forecast[j]))/2));
		errorx += (serie->t[j] - forecast[j]);
		if(j!=0){
			if((serie->t[j] > serie->t[j-1] && forecast[j]>forecast[j-1])
				|| (serie->t[j] < serie->t[j-1] && forecast[j]<forecast[j-1])){
				S->DA +=1.0;
			}
		}
	}
	S->energy /= (serie->size*1.0);
	S->energy *=100.0;
	S->smape /=(serie->size*1.0);
	S->smape *=100.0;
	S->DA *=100.0;
	S->DA /= (serie->size*1.0);
	S->error = errorx/serie->size;
	free(forecast);
}
void sa_copySolution(sa_t *t, sa_state *from, sa_state *to) {
	int i,j;
	for(i=0;i<t->populationSize;i++){
		for(j=0;j<t->dim;j++){
			to->solution[i][j] = from->solution[i][j];
		}
		to->smape[i] = from->smape[i];
		to->energy[i] = from->energy[i];
		to->DA[i] = from->DA[i];
		to->error[i] = from->error[i];
	}
}

void sa_copyOneSolution(sa_t *t,sa_state *state, double* S,int idx){
	int i;
	for(i=0;i<t->dim;i++){
		state->solution[idx][i]=S[i];
	}
}

void sa_printSolution(sa_t *t, sa_state *state) {
	int i,j;
	printf("Solution: \n");
	for(i=0;i<t->populationSize;i++){
		for(j=0;j<t->dim;j++){
			printf("%lf\t",state->solution[i][j]);
		}
		printf(" | MAPE: %lf | SMAPE: %lf | DA: %lf | Error: %lf\n", state->energy[i], state->smape[i],state->DA[i],state->error[i]);
	}
}
void sa_printBestSolution(sa_t *t, sa_state *state, int* bestIdx, float seconds,double bfsenergy, double bfssmape, double bfsDA) {
	FILE *solFile2 = fopen("gcarso.reg","a");
	int i,j;
	printf("Solution: \n");
	double best=state->energy[0];
	*bestIdx=0;
	for(i =1; i < t->populationSize ;i++){
		if(state->energy[i]<best){
			best=state->energy[i];
			*bestIdx=i;
		}
	}
	printf("BESTIDX: %d\n",*bestIdx);
	printf("Alpha %lf Beta %lf\n",state->solution[*bestIdx][0],state->solution[*bestIdx][1]);
	for(j=0;j<t->dim;j++){
		printf("%lf\t",state->solution[*bestIdx][j]);
	}
	printf(" | MAPE: %lf | SMAPE: %lf | DA: %lf | Error: %lf | Time: %lfsecs\n", state->energy[*bestIdx], state->smape[*bestIdx], state->DA[*bestIdx], state->error[*bestIdx],seconds);
	fprintf(solFile2, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", state->energy[*bestIdx],state->smape[*bestIdx],state->DA[*bestIdx],seconds,bfsenergy, bfssmape, bfsDA);
	fclose(solFile2);
}
void sa_findBestSolution(sa_t *t, sa_state *state, int *bestIdx) {
	int i,j;
	double best=state->energy[0];
	*bestIdx=0;
	for(i =1; i < t->populationSize ;i++){
		if(state->energy[i]<best){
			best=state->energy[i];
			*bestIdx=i;
		}
	}
}

double sa_updateTemperature(sa_t *t, double temperature) {
	return t->alpha * temperature;
}

void de_newSample(sa_t *t, sa_state *state, int p0, double F, double CR, ts_sol* S){
	int p1,p2,p3,cutpoint,i;
	int NP = t->dim;
	double *lowerBound =malloc(t->dim * sizeof(double));
	double *upperBound =malloc(t->dim * sizeof(double));
	for(i=0;i<2; i++){
		lowerBound[i]=0.0;
		upperBound[i]=100.0;
	}
	for(i=2;i<t->dim; i++){
		lowerBound[i]=-1.0;
		upperBound[i]=1.0;
	}
	double value;
	do{
		p1=rand()%t->populationSize;
	}while(p1==p0);
	do{
		p2=rand()%t->populationSize;
	}while(p2==p1 || p2==p0);
	do{
		p3=rand()%t->populationSize;
	}while(p3==p0 || p3 == p1 || p3 == p2);
	cutpoint=rand()%NP;
	for(i=0;i<NP;i++){
		if(i==cutpoint || (rand()%100)<CR){
			value=state->solution[p3][i]+ F * (state->solution[p1][i] - state->solution[p2][i]);
			if(value<lowerBound[i]){
				value = lowerBound[i];
			}
			if(value>upperBound[i]){
				value = upperBound[i];
			}
			S->solution[i]=value;
		}else{
			S->solution[i] = state->solution[p0][i];
		}
	}
	free(lowerBound);
	free(upperBound);
}
void ts_read(timeserie *train, timeserie *test){
	FILE *file;
	file = fopen("gcarso.train", "r");
	if(file == NULL) {
		printf("Could not open the train file!\n");
		exit(1);
	}
	int N,i;
	fscanf(file, "%d\n", &N);
	train->size=N;
	train->t = malloc(sizeof(double)*N);
	for(i=0;i<N;i++){
		fscanf(file, "%lf\n", &train->t[i]);
	}
	fclose(file);
	file = fopen("gcarso.test", "r");
	if(file == NULL) {
		printf("Could not open the testing file!\n");
		exit(1);
	}
	fscanf(file, "%d\n", &N);
	test->size=N;
	test->t = malloc(sizeof(double)*N);
	for(i=0;i<N;i++){
		fscanf(file, "%lf\n", &test->t[i]);
	}
	fclose(file);
}

void ts_print(sa_t *t, sa_state *state,timeserie *serie, int bestIdx, int start, float seconds){
	int i,j,k;
	double alpha,beta,fi,ti;
	FILE *solFile = fopen("gcarso.out","w");
	//FILE *solFile2 = fopen(".reg","a");
	double* forecast = malloc(sizeof(double) * serie->size);
	alpha=state->solution[bestIdx][0];
	beta=state->solution[bestIdx][1];
	fi=0;
	ti=start;
	double errorPtg=1;
	for(j=0 ; j< serie->size ;j++){
		fi=0;
		for(k=2 ; k < t->dim; k++){
			fi+=(state->solution[bestIdx][k]*pow(ti,k-2));
		}
		//forecast[j] = alpha*exp(beta*ti) + fi;
		forecast[j] = alpha*exp(beta*ti)+fi;
		//forecast[j]=forecast[j]+(errorPtg*error); //AJUSTE De ERROR
		fprintf(solFile,"%lf, %lf\n", serie->t[j],forecast[j]);
		ti=ti+1;
	}
	double mape=0,smape=0;
	double DA=0;
	for(j = 0 ; j < serie->size ; j++){
		//mape += fabs(serie->t[j] - forecast[j])/fabs(serie->t[j]);
		mape += (fabs(serie->t[j] - forecast[j])/fabs(serie->t[j]));
		//smape += fabs(forecast[j] - serie->t[j])/((fabs(serie->t[j]) + fabs(forecast[j]))/2);
		smape += (fabs(forecast[j] - serie->t[j])/((fabs(serie->t[j]) + fabs(forecast[j]))/2));
		if(j!=0){
			if((serie->t[j] > serie->t[j-1] && forecast[j]>forecast[j-1])
				|| (serie->t[j] < serie->t[j-1] && forecast[j]<forecast[j-1])){
				DA +=1;
			}
		}
	}
	//fprintf(solFile, "SUM = %lf",mape);
	mape /= serie->size;
	smape /= serie->size;
	mape *=100;
	smape *=100;
	DA /= serie->size;
	DA *=100;
	fprintf(solFile, "MAPE: %lf | SMAPE: %lf | DA: %lf | Time: %lf", mape,smape,DA,seconds);
	//fprintf(solFile2, "%lf %lf %lf\n", mape,smape,DA);
	free(forecast);
	fclose(solFile);
	//fclose(solFile2);
}

double expSmooth(double *errorForecast, int n){
	double forecast= errorForecast[0];
	double alpha=0.2;
	int i;
	//printf("Real(fcast): %lf | forecast: %lf\n",errorForecast[0],forecast);
	for(i=1;i<n+1;i++){
		forecast=(alpha*errorForecast[i-1])+((1-alpha)*forecast);
		//printf("Real(fcast): %lf | forecast: %lf\n",errorForecast[i],forecast);
	}
	return forecast;
}

void de_newPob(sa_t *sa, ts_sol *S,timeserie *train,sa_state *next, sa_state *current,
	int iteracion,int windowSize,double errorForecast,double F,double CR, double t,
	double bestEnergyinPob, int *iBetter, ts_sol *bestFoundSolution){
	int i,j;
	double randB;
	double boltzman;
	for(i=0 ; i < sa->populationSize ; i++){
		for(j=0;j<sa->dim;j++){
			S->solution[j]=current->solution[i][j];
		}
		sa_calculateEnergySol(sa,S,train,iteracion, windowSize,errorForecast);
		current->energy[i]=S->energy;
		current->smape[i]=S->smape;
		current->error[i]=S->error;
		current->DA[i]=S->DA;
		de_newSample(sa,current,i,F,CR,S);
		sa_calculateEnergySol(sa,S,train,iteracion,windowSize,errorForecast);

		randB=((float) rand() / (float) RAND_MAX);
		boltzman=exp((current->energy[i]-S->energy) /t);
		if(S->energy < current->energy[i]) { //Comparo Snew <= Sold
			if(S->energy<bestFoundSolution->energy){  //Comparo Snew <= SBest
				bestFoundSolution->energy=S->energy;
				bestFoundSolution->smape=S->smape;
				bestFoundSolution->error=S->error;
				bestFoundSolution->DA=S->DA;
			}
			sa_copyOneSolution(sa,next,S->solution,i);
			next->smape[i]=S->smape;
			next->energy[i]=S->energy;
			next->error[i]=S->error;
			next->DA[i]=S->DA;
		}else{
			if(randB < boltzman){ //Le doy otra oportunidad a Snew con boltzmann
				sa_copyOneSolution(sa,next,S->solution,i);
				next->smape[i]=S->smape;
				next->energy[i]=S->energy;
				next->error[i]=S->error;
				next->DA[i]=S->DA;
			}else{
				sa_copyOneSolution(sa,next,current->solution[i],i);
				next->smape[i]=current->smape[i];
				next->energy[i]=current->energy[i];
				next->error[i]=current->error[i];
				next->DA[i]=current->DA[i];
			}
		}
		if(next->energy[i]<bestEnergyinPob){
			*iBetter=i;
		}
	}
}

int dynEq(sa_state *next, int iBetter,int l,double Lk,int Laux, int iter,
	double Summary,double Esummary,double k1, double k2,double m, double epsilon,int* contadorG){
	int contador;

	if(l > Lk-10){
		iter=iter+1;
		Summary=Summary+iter*next->energy[iBetter];
		Esummary=Esummary+next->energy[iBetter];
		k1=12.0/(double)(iter*iter*iter-iter);
		k2 = 6.0/(double)(iter*iter+iter);
		m=k1*Summary-k2*Esummary;
		if(m<epsilon){ //Equilibrio
			if(l>=Lk){ //Si hay equilibrio y pase Lk acabo el ciclo de metropolis
				l=Lk+Laux;
				Laux=Lk*-1;
			}
		}else{
			if(l+1>=Lk+Laux){ //Aun no llegué al equilibrio
				Laux+=10;
				contador = (*contadorG)%3;;
				*contadorG= (*contadorG)+1;
				//printf("cg %d contador: %d\n",*contadorG,contador);
				if(*contadorG!= 0 && contador==0){
					//l=Lk+Laux;
					Laux=Lk*-1; //BREAK!
					//printf("Contg: %d\n",*contadorG);
				}
			}
		}
	}
	return Laux;
}

int main(){
  //srand(127);
	//gcarso
	clock_t start = clock();
	int contador=0, contadorG=0;
	srand(time(NULL));
	int i,j,k,l; //Loop variables
	int z=0;
	double t; // Temperature.
	double Lk,L1,LMax,beta; // LMax for Metropolis
	int degree=2; //Solution size = Alpha + Beta + degree
	double F=0.6, CR=40; //DE variables
	int bestSolIdx=0;
	double *errorForecast;
	double nextError;
	double suma=0;
	double randB, boltzman;

	int windowSize=50;//Tamaño de la ventana
	int iteracion=0;
	int totalIteraciones;

	L1=1;
	LMax=1400;
	Lk=L1;
	double k1,k2,m; //Stochastic equilibrium
	double Summary=0, Esummary=0;
	int iBetter;
	double bestEnergyinPob;
	int iter=1;
	int Laux=0;
	double epsilon=0.1;
	int pobSize=100;
	sa_t sa; //SA Parameters
	sa_state current, next; //Solution
	timeserie train, test;
	ts_sol S, bestFoundSolution;
	// Initialization.
	ts_read(&train, &test);

	totalIteraciones=ceil(train.size/(windowSize*1.0));
	printf("Total iteraciones: %d\n",totalIteraciones);
	errorForecast = malloc(sizeof(double)*(totalIteraciones));
	for(i=0;i<totalIteraciones;i++){
		errorForecast[i]=0;
	}
	sa_init(&sa, &current, &next, degree, pobSize); //100 es poblacion
	sa_genInitialSolution(&sa, &current, &train, 0,windowSize); //Corregir minimos cuadrados
	S.solution = malloc(sa.dim * sizeof(double));
	bestFoundSolution.solution = malloc(sa.dim *sizeof(double));
	bestFoundSolution.energy = INT_MAX;
	t = sa.t0;
	for(i=0; i < sa.populationSize ; i++){
		for(j=0;j<sa.dim;j++){
			S.solution[j]=current.solution[i][j];
		}
		sa_calculateEnergySol(&sa,&S,&train,0, windowSize, 0);
		current.energy[i]=S.energy;
		current.smape[i]=S.smape;
		current.error[i]=S.error;
		current.DA[i]=S.DA;
		if(S.energy<bestFoundSolution.energy){
			bestFoundSolution.energy=S.energy;
			bestFoundSolution.smape=S.smape;
			bestFoundSolution.error=S.error;
			bestFoundSolution.DA=S.DA;
			errorForecast[0]=S.error;
		}
	}
	beta = exp((log(LMax)- log(L1))/((log(sa.tf) - log(sa.t0))/log(sa.alpha)));
	//Recocido Simulado
	while(t > sa.tf || iteracion != 0) { //CICLO DE TEMPERATURA
		for(l=0;l<Lk+Laux;l++){ //CICLO DE METROPOLIS
			iBetter=0;
			bestEnergyinPob= INT_MAX;
			//POBLACION EVOLUCIONADA CON DE
			de_newPob(&sa,&S,&train,&next, &current,iteracion, windowSize,errorForecast[iteracion],F,CR,t, bestEnergyinPob, &iBetter,&bestFoundSolution);
			//Revision de EQUILIBRIO DINAMICO
			Laux= dynEq(&next, iBetter,l,Lk,Laux,iter,Summary,Esummary, k1,k2,m,epsilon, &contadorG);
		} //CICLO DE METROPOLIS
		Laux=0;
		Summary=0;
		Esummary=0;
		iter=1;
		//Actualizar Lk
		Lk = beta*Lk;
		//Actualizar nueva poblaciÃ³n
		//printf("Current MAPE[0]:\n");
		//printf("MAPE: %lf  SMAPE %lf DA %lf\n",current.energy[10], current.smape[10], current.DA[10]);
		sa_copySolution(&sa, &next, &current);
		//Actualizar temperatura.
		t = sa_updateTemperature(&sa, t);
		//Encontrar mejor solucion y medir el error (para un posterior suavizamiento exponencial)
		sa_findBestSolution(&sa, &current, &bestSolIdx);
		errorForecast[iteracion]=current.error[bestSolIdx];
		iteracion++;
		iteracion=iteracion%totalIteraciones;
		if(iteracion==0){
			printf("%lf > %lf\n ",t,sa.tf);
		}
	} //FIN CICLO TEMPERATURA
	printf("\n\n");
	printf("Best found solution after temperature loop:\n");
	printf("MAPE: %lf  SMAPE %lf \n",bestFoundSolution.energy, bestFoundSolution.smape);
	nextError = expSmooth(errorForecast, totalIteraciones);
	printf("NextError: %lf\n",nextError);
	for(i=0; i < sa.populationSize; i++){
		for(j=0;j<sa.dim;j++){
			S.solution[j]=current.solution[i][j];
		}
		//sa_calculateEnergySolTest(&sa,&S,&test,train.size,nextError);
		sa_calculateEnergySolTest(&sa,&S,&test,train.size,0);
		//printf("Alpha %lf Beta %lf ||",S.solution[0],S.solution[1]);
		current.energy[i]=S.energy;
		current.smape[i]=S.smape;
		current.error[i]=S.error;
		current.DA[i]=S.DA;
		if(S.energy<bestFoundSolution.energy){
			printf("ENTRE!!\n");
			bestFoundSolution.energy=S.energy;
			bestFoundSolution.smape=S.smape;
			bestFoundSolution.error=S.error;
			bestFoundSolution.DA=S.DA;
		}
	}
	// Output.
	printf("\nFinal Best ");
	bestSolIdx=0;
	printf("contadorG: %d",contadorG);
	printf("\n");
	printf("Best found solution during algorithm:\n");
	printf("MAPE: %lf  SMAPE %lf DA %lf\n",bestFoundSolution.energy, bestFoundSolution.smape, bestFoundSolution.DA);
	//printf("best idx: %d\n",bestSolIdx);

	clock_t end = clock();
	float seconds =(float)(end-start)/CLOCKS_PER_SEC;
	sa_printBestSolution(&sa, &current, &bestSolIdx,seconds,bestFoundSolution.energy, bestFoundSolution.smape, bestFoundSolution.DA);
	ts_print(&sa,&current,&test,bestSolIdx, train.size, seconds);
	printf("Time: %lf\n",seconds);
	free(S.solution);
	free(errorForecast);

  return 0;
}
