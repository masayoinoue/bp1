#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include /* Please include "Mersenne Twister" here: a randum number generator  */
/* If not, please use your random number generator */

#define yt 0.75  /* expression threshold */
#define bt pow(10.0,1.5)  /* expression sensitivity */
#define jn 0.0   /* noise amplitude */

#define Prc 0.05 /* network connection probability */
#define Imt 1.0   /* average mutation number in a network */

#define SYS 100 /* number of networks in a system */
#define Np 25   /* number of parental networks */
#define Lni 5   /* number of evolutionary inputs (< Ni) */
#define Ni 5    /* number of input-gene  */
#define Nt 5    /* number of target-genes */
#define Nm 90   /* number of middle-genes */
#define NUM (Ni+Nm+Nt)  /* total number of genes in a network*/

#define LpG 100   /* total number of generations*/
#define TIME 100  /* time for a generation */
#define Sfin 5.0  /* amplitude of external signals */

#define dt 0.5 
#define KK  3  
/*****************************************************************/
int main(int argc, char **argv)
{
  int i,k,l,n,m,j,lp,iv,li,tmpint;
  int Npar[SYS],Pcpl[(int)NUM][(int)NUM];

  double x[(int)NUM],px[(int)NUM],y[(int)NUM],cpl[(int)NUM][(int)NUM],Sext;
  double Prp,Prm,Pmt,tmp;
  double fit[Lni],tfit[SYS],Fpar[SYS],fmax[(int)NUM],var[(int)NUM],uini[(int)NUM];
  double hh[KK],ku1[(int)NUM],ku2[(int)NUM],Noi[(int)NUM];
  
  FILE *file00,*file01,*fop;
  char str00[200],str01[200],strop[200];

  unsigned long init[4]={0x152, 0x364, 0x526, 0x413}, length=4;
  init_by_array(init, length);
  

  Pmt=Imt/((((double)Nm+(double)Ni)*((double)Nm+(double)Nt)-((double)Ni*(double)Nt))*Prc); /* mutation rate */

  /* initialization */
  hh[0]=0.0;
  hh[1]=0.0;
  hh[2]=dt;

  i=0;
  k=0;
  n=0;
  m=0;
  j=0;
  lp=0;
  iv=0;
  
  Sext=0.0;
  for(j=0; j<Lni; j++) fit[j]=0.0; 
  for(j=0; j<SYS; j++){
    Npar[j]=0;
    Fpar[j]=0.0;
    tfit[j]=0.0;
  }
  for(n=0; n<(int)NUM; n++){
    x[n]=0.0;
    px[n]=x[n];
    y[n]=0.0;
    fmax[n]=0.0;
    var[n]=0.0;
    uini[n]=x[n];
    Noi[n]=0.0;
    ku1[n]=0.0;
    ku2[n]=0.0;
  }
  for(n=0; n<(int)NUM; n++){
    for(m=0; m<(int)NUM; m++){
      Pcpl[n][m]=0;
      cpl[n][m]=0.0;
    }
  }


  /* initial networks preparation */
  Prp=0.5*Prc; /* excitatory connection probability */
  Prm=0.5*Prc; /* inhibitory connection probability */

  sprintf(str01, "dattmp");
  file01 = fopen (str01,"w");
  for(j=0; j<SYS; j++){
    fprintf(file01, "%d\n", j);
    for(n=0; n<(int)NUM; n++){
      for(m=0; m<(int)NUM; m++) Pcpl[n][m]=0;
    }
    for(m=0; m<(int)Ni; m++){ /*connection from input to middle genes */
      for(n=Ni; n<(int)(Ni+Nm); n++){
	Rdm=genrand_real2(); /* a random number on [0,1)-real-interval */
	if(Rdm<Prp) Pcpl[n][m]=1;
	else if(Rdm>1.0-Prm) Pcpl[n][m]=-1;
      }
    }
    for(m=Ni; m<(int)(Ni+Nm); m++){ /* connection from middle to target genes */
      for(n=(int)(Ni+Nm); n<(int)NUM; n++){
	Rdm=genrand_real2();
	if(Rdm<Prp) Pcpl[n][m]=1;
	else if(Rdm>1.0-Prm) Pcpl[n][m]=-1;
      }
    }
    for(n=Ni; n<(int)(Ni+Nm); n++){ /* connection from middle to middle genes */
      for(m=Ni; m<(int)(Ni+Nm); m++){
	Rdm=genrand_real2();
	if(Rdm<Prp) Pcpl[n][m]=1;
	else if(Rdm>1.0-Prm) Pcpl[n][m]=-1;
      }
    }
    for(n=0; n<(int)NUM; n++){
      for(m=0; m<(int)NUM; m++) fprintf(file01, "%02d ",Pcpl[n][m]);
      fprintf(file01, "\n");
    }
  }
  fprintf(file01, "\n");
  fclose(file01);
  /* end of initial network preparation */


  sprintf(str00, "datfitness");
  file00 = fopen (str00,"w");

  for(lp=0; lp<LpG; lp++){ /* loop for counting generation number */
    Sext=0.0;
    iv=0;
      
    sprintf(strop, "dattmp");
    if( (fop=fopen(strop, "r"))==NULL){
      printf("no file dattmp\n");
      exit(1);
    }

    for(j=0; j<SYS; j++){  /* loop for counting a network number in a system */
	
      fscanf(fop, "%d\n", &tmpint);
      for(n=0; n<(int)NUM && !feof(fop); n++){
	for(m=0; m<(int)NUM && !feof(fop); m++){
	  fscanf(fop, "%d ", &tmpint);
	  cpl[n][m]=(double)tmpint;
	}
	fscanf(fop, "\n");
      }

      for(li=0; li<Lni; li++){ /*loop for input signals */
	
	Sext=0.0;
	iv=0;
	for(n=0; n<(int)NUM; n++){
	  fmax[n]=0.0;
	  var[n]=0.0;
	  x[n]=genrand_real1(); /* a random number on [0,1]-real-interval */
	  px[n]=x[n];
	  uini[n]=x[n];
	  y[n]=0.0;
	  ku1[n]=0.0;
	  ku2[n]=0.0;
	}
      
	for(i=(int)(-1.1*(double)TIME/dt); i<(int)(1.1*(double)TIME/dt); i++){ /* loop counting a time for a generation */
	  if(i==0){
	    for(n=0; n<(int)NUM; n++){
	      uini[n]=var[n]/(double)iv; /* averaged initial expression levels */
	      fmax[n]=0.0;
	      var[n]=0.0;
	    }
	    Sext=Sfin;
	  }

	  /* time-evolution of expression levels */
	  for(n=0; n<(int)NUM; n++) px[n]=x[n];
	  for(n=0; n<(int)(NUM-Nt); n++){
	    Noi[n]=0.0;
	    if(jn>0.0) Noi[n]=jn*(cos(6.283185307*genrand_real2())*sqrt(-2.0*log(1.0-genrand_real2())));
	  }
	  for(n=0; n<(int)NUM; n++){
	    y[n]=-1.0*yt;
	    if(n==li) y[n]+=Sext;
	    for(m=0; m<(int)NUM; m++) y[n]+=cpl[n][m]*px[m];
	  }
	  for(n=0; n<(int)NUM; n++){
	    ku1[n]=(1.0/(1.0+exp(-1.0*bt*y[n])))-(px[n]);
	  }
	  for(n=0; n<(int)NUM; n++){
	    y[n]=-1.0*yt;
	    if(n==li) y[n]+=Sext;
	    for(m=0; m<(int)NUM; m++) y[n]+=cpl[n][m]*(px[m]+dt*ku1[m]+sqrt(dt)*Noi[m]);
	  }
	  for(n=0; n<(int)NUM; n++){
	    ku2[n]=(1.0/(1.0+exp(-1.0*bt*y[n])))-(px[n]+dt*ku1[n]+sqrt(dt)*Noi[n]);
	  }
	  for(n=0; n<(int)NUM; n++) x[n]=px[n]+(dt/2.0)*(ku1[n]+ku2[n])+sqrt(dt)*Noi[n];
	  /* end: time-evolution of expression levels */

	  iv+=1;
	  for(n=0; n<(int)NUM; n++){
	    if(dt*(double)i<(double)TIME && fmax[n]<fabs(x[n]-uini[n])) fmax[n]=fabs(x[n]-uini[n]); /* maximal expression level */
	    var[n]+=x[n];
	  }
	
	  if(i==(int)(-0.1*(double)TIME/dt) || i==(int)((double)TIME/dt)){
	    for(n=0; n<(int)NUM; n++) var[n]=0.0;
	    iv=0;
	  }

	} /* end loop for i */

	
	/* fitness calculation */
	for(n=0; n<(int)NUM; n++){
	  var[n]=var[n]/(double)iv; /* averaged final expression levels */
	}
	fit[li]=0.0;
	tmp=0.0;
	for(n=0; n<Lni; n++){
	  if(n!=li) tmp+=var[(int)(Ni+Nm+n)];
	}
	fit[li]=(var[(int)(Ni+Nm+li)]-uini[(int)(Ni+Nm+li)])-tmp/(double)(Lni-1);

      }/* end loop for li */
	
      tfit[j]=0.0;
      for(li=0; li<Lni; li++) tfit[j]+=fit[li]/((double)Lni); /* fitness value of each network */

    }/* end loop for j */
      
    fclose(fop);

    
    /* arrange networks in descending order of the fitness */
    for(j=0; j<SYS; j++){
      Npar[j]=0;
      Fpar[j]=-1.0;
    }
    for(j=0; j<SYS; j++){
      if(Fpar[SYS-1]<=tfit[j]){
	Npar[SYS-1]=j;
	Fpar[SYS-1]=tfit[j];
	for(n=1; n<SYS; n++){
	  if(Fpar[SYS-n-1]<=Fpar[SYS-n]){
	    iv=Npar[SYS-n-1];
	    Npar[SYS-n-1]=Npar[SYS-n];
	    Npar[SYS-n]=iv;
	    tmp=Fpar[SYS-n-1];
	    Fpar[SYS-n-1]=Fpar[SYS-n];
	    Fpar[SYS-n]=tmp;
	  }
	}
      }
    }
   
    sprintf(str01, "datnetwork"); /* record networks in the descending order */
    file01 = fopen (str01,"w");
    fprintf(file01, "%d\n", lp+1);

    for(i=0; i<SYS; i++){
      sprintf(strop, "dattmp"); /* open current network */
      if( (fop=fopen(strop, "r"))==NULL){
	printf("no fil dattmp\n");
	exit(1);
      }
      for(j=0; j<1+Npar[i] && !feof(fop); j++){
	fscanf(fop, "%d\n", &iv);
	for(n=0; n<(int)NUM; n++){
	  for(m=0; m<(int)NUM; m++){
	    fscanf(fop, "%d ", &tmpint);
	    Pcpl[n][m]=tmpint;
	  }
	  fscanf(fop, "\n");
	}
	if(iv==Npar[i]){
	  if(i>0) fprintf(file01, "%d\n", i);
	  for(n=0; n<(int)NUM; n++){
	    for(m=0; m<(int)NUM; m++) fprintf(file01, "%02d ",Pcpl[n][m]);
	    fprintf(file01, "\n");
	  }
	}
      }
      fclose(fop);
      
    }
    fprintf(file01, "\n");
    fclose(file01);
    /* end: arrange networks in descending order of the fitness */
      

    fprintf(file00, "%d %lf\n", lp+1, tfit[Npar[0]] ); /* record fitness */
    
 
    /* mutation process */
    sprintf(str01, "dattmp");
    file01 = fopen (str01,"w");
  
    sprintf(strop, "datnetwork");
    if( (fop=fopen(strop, "r"))==NULL){
      printf("no file datnetwork\n");
      exit(1);
    }
      
    for(j=0; j<Np && !feof(fop); j++){ /* take out parantal networks */
      fscanf(fop, "%d\n", &tmpint);
      for(n=0; n<(int)NUM; n++){
	for(m=0; m<(int)NUM; m++){
	  fscanf(fop, "%d ", &tmpint);
	  Pcpl[n][m]=tmpint;
	}
	fscanf(fop, "\n");
      }
      for(i=0; i<(int)(SYS/Np); i++){ /* genrated mutated networks from each parent */
	for(n=0; n<(int)NUM; n++){
	  for(m=0; m<(int)NUM; m++) cpl[n][m]=(double)Pcpl[n][m];
	}
	  
	for(n=0; n<(int)NUM; n++){
	  for(m=0; m<(int)NUM; m++){
	    if(fabs((double)Pcpl[n][m])>0.0){
	      if(genrand_real2()<Pmt){ /* swap a connection */
		iv=0;
		while(iv==0){
		  k=(int)((double)NUM*genrand_real2());
		  l=(int)((double)NUM*genrand_real2());
		  if(Pcpl[k][l]==0 && (int)cpl[k][l]==0){
		    /* from input to middle */
		    if(l<(int)Ni){
		      if(k>=(int)Ni && k<(int)(Ni+Nm)) iv=1;
		    }
		    /* from middle to target */
		    if(l>=(int)Ni && l<(int)(Ni+Nm)){
		      if(k>=(int)(Ni+Nm)) iv=1;
		    }
		    /* from middle to middle */
		    if(l>=(int)Ni && l<(int)(Ni+Nm)){
		      if(k>=(int)Ni && k<(int)(Ni+Nm)) iv=1;
		    }
		  }
		}
		cpl[n][m]=0.0;
		if(genrand_real2()<0.5) cpl[k][l]=1.0;
		else cpl[k][l]=-1.0;
	      }
	    }
	  }
	}/* end loop for n */
	 
	fprintf(file01, "%d\n", (int)((SYS/Np)*j+i) );
	for(n=0; n<(int)NUM; n++){
	  for(m=0; m<(int)NUM; m++) fprintf(file01, "%02d ",(int)cpl[n][m]);
	  fprintf(file01, "\n");
	}
	
      }/* end loop for i */
    }/* end loop for j */
    
    fclose(fop);
    fprintf(file01, "%d\n", lp+1);
    fprintf(file01, "\n");
    fclose(file01);
    
  }/* end loop for lp */
  fclose(file00);
 
}/*main*/
