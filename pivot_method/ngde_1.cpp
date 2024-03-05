// NGDE- Aerosol Dynamics solver.
// Authors: Anand Prakash, Ameya P. Bapat and Michael R. Zachariah, May 01, 2003.

// This C code is used to solve for nucleation, coagulation and surface growth problems. This example problem solves for characteristics of an Aluminum aerosol. However, it can be used for any other material whose properties are known. There are four different sections in this code: 
//1) Coagulation, 
//2) Nucleation (classical) + coagulation, 
//3) Surface growth 
//4) Unified GDE with all the three phenomena combined. 

// The above phenomena constitute four different sections of the code. Each of the four sections are independent of each other. They have especially been written in such manner so that, one can identify the contribution of each phenomena to the GDE. Also, if one requires to use nucleation, coagulation or surface growth alone, it can be easily done.   

// The solution algorithm involves a new approach to solve the Aerosol GDE, where the particle volume space is divided into nodes of zero width. Using nodes allows us to cover the entire size range of aerosol (1 nm to 10 micrometer) using just 40 nodes.A more detailed description of the theory and algorithm can be found in the reference:
//A. Prakash, A. P. Bapat and M.R. Zachariah,  Aerosol Science and Technology

//The main theoretical constraints of this implementation are:
//1. Use of the Self Consistent Classical nucleation Model
//2. Free molecule collision kernel
//3. Free molecule surface growth


// The main variables in the code are v[MAX] and N [MAX]. These arrays are respectively the size (volume), and the number concentration of particles corresponding to that size. v1 is the monomer volume, which we assume, is equal to the volume of a molecule. We have tried to comment the code sufficiently to enable the user to make changes as needed. The main variables used are listed in the table below.

//**************************Note on compiling and running the code**************************
//To compile the code, use command "gcc ngde.c -lm". To run the code the user should edit the input file as required and then the command "a.out > outputfile" should be used to direct the output to a file.
//******************************************************************************************

// *************************************NOTE ON UNITS***************************************
//                            SI Units are used throught the code.
//******************************************************************************************
/*                                   ------------                                            */
/* 				     NOMENCLATURE                                            */
/* 				     ------------                                            */

/* __________________________________________________________________________________________*/
/* Variables            	            Meaning (Units)                                  */
/* __________________________________________________________________________________________*/
/* v[i]         	Volume of particles at ith node (m3/m3 of aerosol)                   */
/* N[i]         	Number of particles at ith node (#/m3 of aerosol)                    */
/* X[i][j][k]   	Splitting operator for coagulation (dimensionless)                   */
/* K[i][j]      	Collision frequency function(m3/s)                                   */
/* N1s[I]       	Saturation monomer concentration over particle of size I(#/m3)       */
/* sum1,sum2    	Addition and subtraction terms in Nk due to collision(#/(m3s))       */
/* Vav          	Average volume of particles (m3/m3 of aerosol)                       */
/* Vtotal       	Total volume of particles and monomer (m3/m3 of aerosol)             */
/* T                    Temperature	(K)                                                  */
/* coolrate     	Cooling rate	(K/s)                                                */
/* t            	Time	(s)                                                          */
/* step	                Timestep for integration  (s)                                        */
/* Ninf	                Total number of particles as got by solving GDE	(#/m3 of aerosol)    */
/* Ninf_eq       	Total number of particles by eq 7.77, Friedlander (#/m3 of aerosol)  */
/* fi           	Total volume of particles (m3/m3 of aerosol)                         */
/* Vtot	                Total volume of particles (excluding monomers)	(m3/m3 of aerosol)   */
/* addterm,subterm	Addition and subtraction terms in Nk due to surface growth(#/(m3s))  */
/* Ntot         	Total number of particles (#/m3 of aerosol)                          */
/* dpav         	Mean Diameter	(m)                                                  */
/* v1           	Monomer volume	(m3/m3 of aerosol)                                   */
/* sigma        	Surface Tension	(N/m)                                                */
/* Ps           	Saturation Pressure (Pa)                                             */
/* ns           	Monomer concentration at saturation(#/m3 of aerosol)                 */
/* kstar        	Critical node size (dimensionless)                                   */
/* dpstar       	Critical particle size (m)                                           */
/* vstar        	Volume of particles in the critical bin	(m3/m3 of aerosol)           */
/* s1           	Surface area of monomer	(m2)                                         */
/* S            	Saturation ratio (dimensionless)                                     */
/* m1           	Monomer mass(Kg)                                                     */
/* zeta[i]      	Splitting operator for nucleation at the ith node (dimensionless)    */
/* t_SPD                Time requred to reach a self-preserving distribution (s)             */

//*****************************EXAMPLE PROBLEM*********************************
// Aluminum vapor at a high temperature (1773 K) is being cooled at a rate ~1000 K/s. As the system cools down, supersaturation of vapor causes nucleation of particles, which then coagulate to form larger particles. Also, due to high concentration of monomer units (molecules), supersaturation causes surface growth of existing particles. The monomers either evaporate from the surface of existing particles, or condense on them.
//******************************************************************************


// Authors: Anand Prakash, Ameya P. Bapat and Michael R. Zachariah, Nov 27, 2002.

// Beginning of the code

#include <iostream>
#include <fstream>

#include<stdio.h>	
#include<math.h>
#define	pi	3.1415926	//pi
#define kb  	1.38e-23	//Boltzmann constant
#define R       8.314           //Universal gas constant Joule/(K.mol)
#define Na      6.0225e23       //Avogadro's Number
#define MAX1    100

 int main()
{
	int i, j, k, choice, counter, counter2, opt, printopt, flag, MAX, beta_option;

	//added dynamic allocation of memory for size splitting operator X
	double ***X;
	X = new double **[MAX1]();
	for (i = 0; i < MAX1; i++)
	{
		X[i] = new double *[MAX1]();
		for (j = 0; j < MAX1; j++)
			X[i][j] = new double[MAX1]();
	}

  // Declaration of variables
  double v[MAX1+1],N[MAX1],K[MAX1][MAX1],N1s[MAX1],sum1,sum2,temp,temp1,temp2,Vav,Vtotal;
  double t,step,temp3,Ninf,n,Vtot,addterm,subterm,Ntot,dpav,v1,fi,t_SPD,step_coag,q,P;
  double sigma,Ps,ns,theta,a,b,Jk,kstar,dpstar,vstar,s1,S,m1,c,zeta[MAX1],dp[MAX1],t1,t2,t3,t4,Kmin;
  double A,B,C,D,E,F,MW,rho,Ninf_eq,T,coolrate,d;
  double lambda,kn1,kn2,mu,D1,D2,m[MAX1],c1,c2,g1,g2,l1,l2,A_mu,B_mu;
  FILE *fptr,*fptr1,*fptr2;
  if((fptr=fopen("main.inp","r+"))==NULL)//Opening the input file that contains property data for the aerosol material.
  {
    printf("Warning: Could not open file");
    return(-1);
  }

  double time_end = 0.5;
  //Scanning the property data from the input file "main.inp"
  fscanf_s(fptr,"Number of nodes: %d Starting Temperature(Kelvin): %lf Cooling Rate(Kelvin/s): %lf Pressure(Atm): %lf Molecular Weight(Kg/mol): %lf Density(Kg/m3): %lf Surface Tension in the form A-BT(dyne/cm): A=%lf B=%lf Saturation Vapor Pressure of material in form Ps=exp(C-D/T)P: C=%lf D=%lf Enter your choice - 1)Coagulation,2)nucleation+coagulation,3)nucleation + coagulation + surface growth,4)surface growth:%d Choice of collision frequency function:1)Free Molecular Regime,2)Fuchs form for free molecular and transition regime:%d Sutherland's constants for carrier gas in the form A*T^1.5/(B+T): A=%lf B=%lf",&MAX,&T,&coolrate,&P,&MW,&rho,&A,&B,&C,&D,&choice,&beta_option,&A_mu,&B_mu);
  fclose(fptr);//Closing the main property data file.
  q=pow(10.0,(12.0/(MAX-2)));//Calculation of geometric spacing factor which depends on the number of nodes.

  v1 = MW/(rho*Na);// Volume of a monomer unit (a molecule); Na is the Avogadro's Number.
  
  //Ergebnis
  std::ofstream resultPSD;
  resultPSD.open("SizeDistribution.txt");


  //**********************************************
  //           PURE COAGULATION
  //**********************************************
  // This part of the code solves for coagulation problems. The code verifies that we obtain a self preserving distribution. It verifies the problem stated in Friedlander, 2000 - page 218, line 4 - and shows that SPD is obtained. And note that in this case node 1 is not the monomer. Its particle of size 1 nm. We do not have to worry about monomers, since we are concerned with coagulation only.
  if (choice==1)
  {
    if((fptr1=fopen("coag.inp","r+"))==NULL)//opening the input file for coagulation probelm.
    {
      printf("Warning: Could not open file");
      return(-1);
    }
      for(i=1;i<=MAX-1;i++)
	N[i]=0.0;
      //Scanning the initial monomer size and number concentration of the monodisperse distribution for coagulation.
	  fscanf_s(fptr1,"\n Initial diameter of monodisperse particles(nm): %lf  Initial number concentration(#/m3): %lf Temperature(Kelvin): %lf",&d,&N[1],&T);
      v[1] = 1e-27*pi*pow(d,3)/6;// Setting diameter of particle that the user enters as the smallest node( coagulation would never lead to decrease in particle size).

	  std::cout << "Node 0:" << v[0] << std::endl;
	  std::cout << "Node 1: " << v[1] << std::endl;

      dp[1]=pow((6*v[1]/pi),0.3333333);
       m[1]=v[1]*rho; 
      for(i=2;i<=MAX;i++) // Calculating volumes of larger nodes using a geometric spacing factor of q.
	{
	  v[i]=v[1]*pow(q,(i-1));	

	  std::cout << "Node: " << i << " Vol: " << v[i] << std::endl;

	  dp[i]=pow((6*v[i]/pi),0.3333333);
	  m[i]=v[i]*rho;
	}
      //Calculation of size splitting operators.
      for(k=1;k<=MAX-1;k++)
	{
	  for(i=1;i<=MAX-1;i++)
	    {
	      for(j=1;j<=MAX-1;j++)
		{
		  //Conditions in parentheses check if the combined volume of colliding particles is between k and k+1.
		  if(v[k]<= (v[i]+v[j]) && (v[i]+v[j])<v[k+1])
		    {
		      X[i][j][k]=(v[k+1]-v[i]-v[j])/(v[k+1]-v[k]);
		    }	     
		  else
		    {
		      if(v[k-1]<= (v[i]+v[j]) && (v[i]+v[j])<v[k])
			{
			  X[i][j][k]=(v[i]+v[j]-v[k-1])/(v[k]-v[k-1]);		      
			}
		      else 
			X[i][j][k]=0;
		    }	      
		}
	    }
	}
      t=0.0; // Initializing time.
      step=1e-12;//Initializing timestep for integration. 
      t_SPD=5/(pow(3.0/(4.0*pi),0.1666666)*pow((6.0*kb*T/rho),0.5)*pow(v[1],0.1666666)*N[1]);//Estimation of time required to reach SPD.
      printf("\n\nEstimated time required to reach SPD=%es\n",t_SPD);
      Ninf_eq=1e24;// A temporary variable that calculates N_total according to Friedlander, 2000, equation 7.77.
     // printf("\ntime\t\t Kmin\t\t Ninf\t\tN[1]\t\tN[5]\t\tN[10]\t\tN[15]");


      while(t<time_end)//Running the coagulation mode of the code until SPD is reached. The dimensionless number distribution remains same after SPD has reached.
	{
	  //Calculation of collision frequency function ( same as in " nucleation + coagulation" section of the code).
	  Kmin=1e-9;//setting an arbitrary value for Kmin to start with. 
	  if(beta_option==1)
	    {
	      for(i=1;i<=MAX-1;i++)
		{
		  for(j=1;j<=MAX-1;j++)
		    {
		      
		      temp1=1/v[i] + 1/v[j];
		      temp2=pow(v[i],0.333333) + pow(v[j],0.333333);
		      K[i][j]= pow(3.0/(4.0*pi),0.1666666)*pow((6.0*kb*T/rho),0.5)*pow(temp1,0.5)*pow(temp2,2.0);	     
		      if (K[i][j]<Kmin)
			Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time. 
		      
		    }
		}
	    }
	  if(beta_option==2)
	    {
	      mu=A_mu*pow(T,1.5)/(B_mu+T);
	      lambda = (mu/(P*101325.0))* sqrt(pi*R*T/(2.0*0.04));
	      for(i=1;i<=MAX-1;i++)
		{
		  for(j=1;j<=MAX-1;j++)
		    {
		      
		      kn1=(2.0*lambda)/dp[i];
		      kn2=(2.0*lambda)/dp[j];
		      D1=(kb*T)/(3.0*pi*mu*dp[i])*((5.0+4.0*kn1+6.0*kn1*kn1+18.0*kn1*kn1*kn1)/(5.0-kn1+(8.0+pi)*kn1*kn1));
		      D2=(kb*T)/(3.0*pi*mu*dp[j])*((5.0+4.0*kn2+6.0*kn2*kn2+18.0*kn2*kn2*kn2)/(5.0-kn2+(8.0+pi)*kn2*kn2));
		      c1=sqrt((8.0*kb*T)/(pi*m[i]));
		      c2=sqrt((8.0*kb*T)/(pi*m[j]));
		      l1=(8.0*D1)/(pi*c1);
		      l2=(8.0*D2)/(pi*c2);
		      g1=(pow((dp[i]+l1),3)-pow((dp[i]*dp[i]+l1*l1),1.5))/(3.0*dp[i]*l1)-dp[i];
		      g2=(pow((dp[j]+l2),3)-pow((dp[j]*dp[j]+l2*l2),1.5))/(3.0*dp[j]*l2)-dp[j];
		      K[i][j]=2.0*pi*(D1+D2)*(dp[i]+dp[j])/((dp[i]+dp[j])/(dp[i]+dp[j]+2.0*sqrt(g1*g1+g2*g2))+(8.0*(D1+D2))/(sqrt(c1*c1+c2*c2)*(dp[i]+dp[j])));
		      //  printf("kn %e\n",kn2);
		      if (K[i][j]<Kmin)
			Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time. 
		      
		    }
		}
	    }
	  
	  //Calculating the gain and loss terms due to coagulation.
	  for(k=1;k<=MAX-1;k++)
	    {
	      sum1=0;//Addition term when i and j collide to form a k sized particle.
	      sum2=0;//Subtraction term when k collides with any other particle.
	      for(i=1;i<=MAX-1;i++)
		{
		  sum2=sum2+K[k][i]*N[i];
		  for(j=1;j<=k;j++)
		    {
		      sum1=sum1+X[i][j][k]*K[i][j]*N[i]*N[j];
		    }
		}
	      N[k]= N[k] + step*(0.5*sum1 - N[k]*sum2);
	    }
	  Ninf=0;// Total number of particles, which is zero initially.
	  for(i=1;i<=MAX-1;i++)
	    Ninf = Ninf + N[i];//N_infinity as we get by solving the coagulation part of GDE.
	  step = 1e-3/(Kmin*Ninf);//Adaptive timestep for integration, based on characteristic coagulation time.
	  fi=0;
	  for(i=1;i<=MAX-1;i++)
	    fi = fi + N[i]*v[i];
	  //N_infinity according to equation 7.77 in Friedlander.
	  Ninf_eq = Ninf_eq - step*(3.33*pow((3/(4*pi)),0.1666666)*pow(6*kb*T/rho,0.5)*pow(fi,0.1666666)*pow(Ninf_eq,(11.0/6.0)));
	  Vtot=0;
	//  printf("\n %e \t %e \t %e \t %e \t %e\t %e\t %e",t,Kmin,Ninf,N[1],N[5],N[10],N[15]);
	  t=t+step;
	}	

      // Printing the final number distribution after SPD is reached.
      printf("\n\n***Size Distribution after %2.2e****\n", t);
      printf("\nvolume\t\tnumber");
	
	  for (i = 1; i <= MAX - 1; i++)	  {

		  printf("\n%e\t%e", v[i], N[i]);
		  resultPSD << i <<"\t" << v[i] << "\t" << N[i] << std::endl;
	  }

      fclose(fptr1);
    }
  //************************************************************
  //             Nucleation + Coagulation
  //************************************************************
  
 else if (choice==2) //If user chooses 2, code runs for nucleation + coagulation.
    {
      v[1]=v1;//Assigning first node to monomer. Note that these are not considered as particles.
      for(i=1;i<=MAX-1;i++)
	N[i]=0.0;//Setting initial number of particles = 0 in all the nodes.		
      for(i=2;i<=MAX;i++)
	//Defining volume spacing of bins. Using a geometric factor of 2 to cover 1nm to 10 micrometer size particles.
	{
	  v[i]=v[1]*pow(q,(i-1));
	  m[i]=rho*v[i];
	  dp[i]=pow(6.0*v[i]/pi,0.3333333);
	}
      dp[1]=pow(6.0*v[1]/pi,0.3333333);
      m[1]=rho*v[1];
      //Calculation of splitting operator for coagulation, X[i][j][k].
      for(k=1;k<=MAX-1;k++)
	{
	  for(i=1;i<=MAX-1;i++)
	    {
	      for(j=1;j<=MAX-1;j++)
		{
		  //Conditions in parentheses check if the combined volume of colliding particles is larger than k.
		  if(v[k]<= (v[i]+v[j]) && (v[i]+v[j])<v[k+1]) 
		    {
		      X[i][j][k]=(v[k+1]-v[i]-v[j])/(v[k+1]-v[k]);
		    }	     
		  else
		    {
		      if(v[k-1]<= (v[i]+v[j]) && (v[i]+v[j])<v[k])
			{
			  X[i][j][k]=(v[i]+v[j]-v[k-1])/(v[k]-v[k-1]);		      
			}
		      else 
			X[i][j][k]=0;
		    }	      
		}
	    }
	}
      t=0.0;//Initializing time.
      step=5e-4;//Timestep for integration.
      S=1.001;// Setting initial saturation ratio, a little larger than 1. Hereafter the saturation ratio is determined by cooling or heating rate of the aerosol. A heating rate would require a negative value.

      // Calculating surface area of the monomer, given the volume of monomer v1.
      temp1 = 6*v1/pi;
      s1 = pi*pow(temp1,0.66666666667);
      //Calculating mass of the monomer unit.
      m1=rho*v1;
      //Calculating saturation vapor pressure for the Aerosol.
      Ps=exp(13.07-36373.0/T)*101325.0*P;
      //Calculating number concentration of monomers at saturation.
      ns=Ps/(kb*T);
      // Calculating monomer concentration ( multiplying saturation ratio by saturation concentration of monomers).
      N[1]=S*ns;
      counter=0; // Temporary variable for printing after a large number of iterations.
      printf("\ntime\t\tMonomer\t\tJk\t\tS\t\tDp_mean\t\tk*\t\tNtot");
      while(T>300) // Calculating aerosol properties until the system cools down to 27 C.
	{
	  sigma=(A-B*T)*1e-3; // Surface tension. 
	  Ps=exp(C-D/T)*101325.0*P; // Saturation pressure.
	  ns=Ps/(kb*T); // Saturation concentration of monomers using ideal gas law.
	  S=N[1]/ns; // Saturation ratio.
	  theta = (s1*sigma)/(kb*T);// Dimensionless surface tension theta.
	  a=(2*sigma)/(pi*m1); //Temporary variable.
	  b= theta - (4*pow(theta,3))/(27*pow(log(S),2));//Temporary variable.
	  Jk = pow(ns,2)*S*v1*pow(a,0.5)*exp(b); //Nucleation rate using the classical SCC model.
	  c = 0.6666667*theta/log(S);// Temporary variable.
	  kstar = pow(c,3.0);//Calculating critical cluster size that will determine the node at which nucleation occurs.
	  dpstar= 4*sigma*v1/(kb*T*log(S));//Size of particles corresponding to the critical node.
	  vstar=pi*pow(dpstar,3)/6; // Volume of particle corresponding to the critical node.
	  Ntot=0.0; // Initializing total number of particles to zero, as there are no particles prior to nucleation.
	  for(i=2;i<=MAX-1;i++)
	    Ntot = Ntot + N[i];//Total number of particles (does not include monomers).
	  Vtot=0.0;
	  for(i=2;i<=MAX-1;i++)
	    Vtot = Vtot + N[i]*v[i];//Total volume of particles. Note that loop runs from i=2 because i=1 corresponds to monomers, which we do not count as particles.
	  Vtotal = Vtot + N[1]*v[1];//Total volume for mass conservation check (includes monomers) 
	  Vav = Vtot/Ntot; // Average volume of particles ( number average)
	  dpav = pow((6*Vav/pi),0.3333333);//Volume based mean diameter of particles printing after every 50 times the loop runs.
	  counter=counter + 1; 
	  if (counter==50)
	    {
	      printf("\n%e\t%e\t%e\t%e\t%e\t%e\t%e",t,N[1],Jk,S,dpav,kstar,Ntot);
	      counter=0;
	    }
	     //Calculation of collision frequency function beta(i,j)=K(i,j).
	 
	  if(beta_option==1)
	    {
	      for(i=1;i<=MAX-1;i++)
		{
		  for(j=1;j<=MAX-1;j++)
		    {
		      
		      temp1=1/v[i] + 1/v[j];
		      temp2=pow(v[i],0.333333) + pow(v[j],0.333333);
		      K[i][j]= pow(3.0/(4.0*pi),0.1666666)*pow((6.0*kb*T/rho),0.5)*pow(temp1,0.5)*pow(temp2,2.0);	     
		      if (K[i][j]<Kmin)
			Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time. 
		      
		    }
		}
	      }
	   if(beta_option==2)
	  {
	      mu=A_mu*pow(T,1.5)/(B_mu+T);
	      lambda = (mu/(P*101325.0))* sqrt(pi*R*T/(2.0*0.04));
	      for(i=1;i<=MAX-1;i++)
		{
		  for(j=1;j<=MAX-1;j++)
		    {
		      kn1=(2.0*lambda)/dp[i];
		      kn2=(2.0*lambda)/dp[j];
		      D1=(kb*T)/(3.0*pi*mu*dp[i])*((5.0+4.0*kn1+6.0*kn1*kn1+18.0*kn1*kn1*kn1)/(5.0-kn1+(8.0+pi)*kn1*kn1));
		      D2=(kb*T)/(3.0*pi*mu*dp[j])*((5.0+4.0*kn2+6.0*kn2*kn2+18.0*kn2*kn2*kn2)/(5.0-kn2+(8.0+pi)*kn2*kn2));
		      c1=sqrt((8.0*kb*T)/(pi*m[i]));
		      c2=sqrt((8.0*kb*T)/(pi*m[j]));
		      l1=(8.0*D1)/(pi*c1);
		      l2=(8.0*D2)/(pi*c2);
		      g1=(pow((dp[i]+l1),3)-pow((dp[i]*dp[i]+l1*l1),1.5))/(3.0*dp[i]*l1)-dp[i];
		      g2=(pow((dp[j]+l2),3)-pow((dp[j]*dp[j]+l2*l2),1.5))/(3.0*dp[j]*l2)-dp[j];
		      K[i][j]=2.0*pi*(D1+D2)*(dp[i]+dp[j])/((dp[i]+dp[j])/(dp[i]+dp[j]+2.0*sqrt(g1*g1+g2*g2))+(8.0*(D1+D2))/(sqrt(c1*c1+c2*c2)*(dp[i]+dp[j])));
		      if (K[i][j]<Kmin)
			Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time. 
		    }
		}
	      }
    
	  //Operator to put nucleated particles in the bin just higher than k*.
	  for(k=2;k<=MAX-1;k++)
	    {
	      if (vstar<v1)
		zeta[2]=vstar/v[2];//Putting particles formed smaller than monomers in the smallest particle node (node 2). This situation arises when k* falls below monomer size.
	      else if (v[k-1] <=vstar && vstar < v[k]) 
		zeta[k]=vstar/v[k];// Putting particles in node just larger than k*.	      
	      else
		zeta[k]=0.0;
	    }
	  //Calculation of gain and loss terms for Nk due to coagulation.
	  for(k=2;k<=MAX-1;k++)
	    {
	      //Initializing gain and loss terms for coagulation.
	      sum1=0;// Production term due to collision of two smaller particles to form a k size particle.
	      sum2=0;// Loss term due to collision of k size particle with any other particle.
	      for(i=2;i<=MAX-1;i++)
		{
		  sum2=sum2+K[k][i]*N[i]; // k collides with any other particle to get out of node k, thus loss term.
		  for(j=2;j<=k;j++)
		    sum1=sum1+X[i][j][k]*K[i][j]*N[i]*N[j];// i and j collide to form particle of size k which is then multiplied by the size splitting operator to put the particles in the adjacent nodes after adjusting the volume.
		}
	      // Change in PSD due to nucleation + coagulation.
	      N[k]= N[k] + step*(0.5*sum1 - N[k]*sum2 + Jk*zeta[k]);
	    }
	  //Monomer balance. accounting for the monomers loss due to nucleation.
	  N[1] = N[1] - Jk*kstar*step;
	  T=T-step*coolrate;// Temperature would decrease as time progresses, due to cooling.
	  t=t+step;//Time increment.
	}
      printf("\n\nFinal Particle Size Distribution\n\n");
      printf("\nvolume\t\tnumber");
      for(i=2;i<=MAX-1;i++)
	printf("\n%e\t%e",v[i],N[i]);// Printing the final number distribution obtained due to nucleation and coagulation. (number concentration or particles against nodes).
    }


  //*****************************************
  //SURFACE GROWTH + NUCLEATION + COAGULATION
  //*****************************************
  else if(choice==3) // THE COMPLETE GDE
    {
      v[1]=v1; // Setting volume of the first node equal to the monomer volume, i.e. volume of a molecule, which is v1.
      for(i=1;i<=MAX-1;i++)
	{
	  N[i]=0.0; // Setting number of particles = 0 at all the nodes initially. 
	  N1s[i]=0.0; // Setting number concentration of monomers over different size particles = 0 initially.
	}
      for(i=2;i<=MAX;i++)
	{
	v[i]=v[1]*pow(q,(i-1)); // Calculating volume of all larger nodes using a geometric factor of 2. Note that volume of node 1 is the volume of a molecule.
	dp[i]=pow((6.0*v[i]/pi),0.333333);
	m[i]=rho*v[i];
	}
      	dp[1]=pow((6.0*v[1]/pi),0.333333);
	m[1]=rho*v[1];
      //Size splitting operator for coagulation ( Same as in "coagulation + nucleation" section of the code).
      for(k=1;k<=MAX-1;k++)
	{
	  for(i=1;i<=MAX-1;i++)
	    {
	      for(j=1;j<=MAX-1;j++)
		{
		  if(v[k]<= (v[i]+v[j]) && (v[i]+v[j])<v[k+1])
		    {
		      X[i][j][k]=(v[k+1]-v[i]-v[j])/(v[k+1]-v[k]);
		    }	     
		  else
		    {
		      if(v[k-1]<= (v[i]+v[j]) && (v[i]+v[j])<v[k])
			{
			  X[i][j][k]=(v[i]+v[j]-v[k-1])/(v[k]-v[k-1]);		      
			}
		      else 
			X[i][j][k]=0;
		    }      
		}
	    }
	}
      t=0.0;// Initializing time.
      step=1e-4; // Initial time step for integration.
      S=1.001; // Initial saturation ratio. Assuming that the vapor is just saturated initially.
      temp1 = 6*v1/pi; // Temporary variable.
      s1 = pi*pow(temp1,0.66666666667);// Surface area of monomer.
      m1=rho*v1;// Mass of a monomer.
      Ps=exp(C-D/T)*101325.0*P;// Saturation vapor pressure at temperature T.
      ns=Ps/(kb*T);// Saturation concentration of monomers using ideal gas law.
      N[1]=S*ns;// Monomer concentration. 
      printf("\nTime\t\tN[1]\t\tJk\t\tS\t\tdiameter\tParticle Volume\tTotal Number");// Printing result heading.
      
      // Counters below are temporary variables that determine number of time the loop runs before every print step. Generally time steps are small and printing results is done after every 4000 steps.
      counter=1;
      counter2 = 1;
      while(T>300) // Observing changes until temperature drops down to room temperature ~ 27 C.
	{
	  // Terms below have been explained in "nucleation + coagulation" section of the code above.
	  sigma=(A-B*T)*1e-3;
	  Ps=exp(C-D/T)*101325.0*P;
	  ns=Ps/(kb*T);
	  S=N[1]/ns;
	  theta = (s1*sigma)/(kb*T);
	  a=(2*sigma)/(pi*m1);
	  b= theta - (4*pow(theta,3))/(27*pow(log(S),2));
	  Jk = pow(ns,2)*S*v1*pow(a,0.5)*exp(b); 
	  c = 0.6666667*theta/log(S);
	  kstar = pow(c,3.0);
	  dpstar= 4*sigma*v1/(kb*T*log(S));
	  vstar=pi*pow(dpstar,3)/6;
	  Ntot=0.0;
	  for(i=2;i<=MAX-1;i++)
	    Ntot = Ntot + N[i];//Total number of particles.
	  Vtot=0.0;
	  Vtotal=0.0;
	  for(i=2;i<=MAX-1;i++)
	    Vtot = Vtot + N[i]*v[i];//Total volume of particles . This volume is used to evaluate the mean volume and thus mean diameter of particles.
	  Vtotal = Vtot + N[1]*v[1];//Total mass volume (including monomers) This volume is used to check for mass conservation.
	  Vav = Vtot/Ntot; // Mean volume of particles.
	  dpav = pow((6*Vav/pi),0.3333333);//Mean diameter of particles.
	  if(Ntot>100.0)//Reducing the timestep for integration when there are sufficient number of particles such that surface growth takes over and the code can no longer run for large time-steps.
	    step=1e-5;
	  if(t>0.14)
	    step=5e-8;
	  counter = counter +1;
	  counter2 = counter2 +1;
	  if(counter == 400)// Printing mean diameter etc. after every 4,000 steps.
	    {
	      printf("\n%e\t%e\t%e\t%e\t%e\t%e\t%e",t,N[1],Jk,S,dpav,Vtot,Ntot);
	      counter = 1;
	    }
	  if (counter2 == 40000)// Printing PSD after every 400,000 steps.
	    {
	      for (i=1;i<=MAX-1;i++)
		printf("\nN%d\t%e",i,N[i]);
	      counter2=1;
	    }
	  //Collision frequency function for coagulation.
	  if(beta_option==1)
	    {
	      for(i=1;i<=MAX-1;i++)
		{
		  for(j=1;j<=MAX-1;j++)
		    {
		      
		      temp1=1/v[i] + 1/v[j];
		      temp2=pow(v[i],0.333333) + pow(v[j],0.333333);
		      K[i][j]= pow(3.0/(4.0*pi),0.1666666)*pow((6.0*kb*T/rho),0.5)*pow(temp1,0.5)*pow(temp2,2.0);	     
		      if (K[i][j]<Kmin)
			Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time. 
		      
		    }
		}
	      }
	   if(beta_option==2)
	  {
	      mu=A_mu*pow(T,1.5)/(B_mu+T);
	      lambda = (mu/(P*101325.0))* sqrt(pi*R*T/(2.0*0.04));
	      for(i=1;i<=MAX-1;i++)
		{
		  for(j=1;j<=MAX-1;j++)
		    {
		      kn1=(2.0*lambda)/dp[i];
		      kn2=(2.0*lambda)/dp[j];
		      D1=(kb*T)/(3.0*pi*mu*dp[i])*((5.0+4.0*kn1+6.0*kn1*kn1+18.0*kn1*kn1*kn1)/(5.0-kn1+(8.0+pi)*kn1*kn1));
		      D2=(kb*T)/(3.0*pi*mu*dp[j])*((5.0+4.0*kn2+6.0*kn2*kn2+18.0*kn2*kn2*kn2)/(5.0-kn2+(8.0+pi)*kn2*kn2));
		      c1=sqrt((8.0*kb*T)/(pi*m[i]));
		      c2=sqrt((8.0*kb*T)/(pi*m[j]));
		      l1=(8.0*D1)/(pi*c1);
		      l2=(8.0*D2)/(pi*c2);
		      g1=(pow((dp[i]+l1),3)-pow((dp[i]*dp[i]+l1*l1),1.5))/(3.0*dp[i]*l1)-dp[i];
		      g2=(pow((dp[j]+l2),3)-pow((dp[j]*dp[j]+l2*l2),1.5))/(3.0*dp[j]*l2)-dp[j];
		      K[i][j]=2.0*pi*(D1+D2)*(dp[i]+dp[j])/((dp[i]+dp[j])/(dp[i]+dp[j]+2.0*sqrt(g1*g1+g2*g2))+(8.0*(D1+D2))/(sqrt(c1*c1+c2*c2)*(dp[i]+dp[j])));
		      if (K[i][j]<Kmin)
			Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time. 
		    }
		}
	      }
	  //Putting nucleated particles in bins just larger than k* ( explained in section " coagulation + nucleation" above).
	  for(k=2;k<=MAX-1;k++)
	    {
	      if (vstar<v1)
		zeta[2]=vstar/v[2];	   
	      else if (v[k-1] <=vstar && vstar < v[k]) 
		zeta[k]=vstar/v[k];	      
	      else
		zeta[k]=0.0;
	    }
	  //Calculating the saturation monomer concentration over an i sized particle that has a diameter dp[i].
	  for(i=1;i<=MAX-1;i++)
	    {
	      dp[i]=pow((6*v[i]/pi),0.333333); // Diameter of particle at the ith node.
	      N1s[i]=ns*exp(4*sigma*MW/(R*T*rho*dp[i]));// Calculating saturation monomer concentration over a i size particle including the kelvin effect.
	    }
	  // t1 through t4 are terms used for monomer balance. These are loss and gain terms due to condensation and evaporation respectively. Each of these terms has been explained below as I use them. Initially setting these terms = 0.
	  t1=0;
	  t2=0;
	  t3=0;
	  t4=0;
	  // Note that in the following part "addterm" and "subterm" are addition and subtraction term in particle size distribution due to condensation or evaporation (surface growth), where as terms t1 through t4 are addition and subtraction term for monomer balance.
	  for(k=2;k<=MAX-1;k++)
	    {
	      sum1=0;//Addition term due to coagulation of size i+j to form k.
	      sum2=0;//Loss term due to coagulation of k with any other particle.
	      if(step!=1e-4)//Activating surface growth when sufficient number of particles have nucleated. At smaller times there are very small number of particles and surface growth hardly changes PSD. Doing so makes the code run faster.
		{
		  // The following four "if" statements check if the vapor pressure of monomer in the system is larger or smaller than saturation vapor pressure over the k size particle. The difference in the two vapor pressures is the driving force for evaporation or condensation. 
		  if(N[1]>N1s[k-1]) // i.e. if actual monomer concentration in the aerosol is greater than saturation concentration over a k-1 size particle.
		    {
		      if (k==2) 
			addterm=0.0; // There are no particles smaller than 2nd node. so particles in 2nd node can not be added due to growth of smaller particles.
		      else
			{
			  addterm = (v[1]/(v[k]-v[k-1]))*K[1][k-1]*(N[1]-N1s[k-1])*N[k-1];//Growth of k due to condensation of monomers on k-1.
			  t1 = t1 + K[1][k-1]*(N[1]-N1s[k-1])*N[k-1];//Loss of monomers that have condensed.
			}
		    }
		  // Similarly, we can explain the following terms below.
		  if(N[1]<N1s[k+1])
		    {
		      addterm = -(v[1]/(v[k+1]-v[k]))*K[1][k+1]*(N[1]-N1s[k+1])*N[k+1];//Growth of k due to evaporation of monomers from k+1.
		      t2 = t2 + K[1][k+1]*(-N[1]+N1s[k+1])*N[k+1];//Gain of monomers that have evaporated.
		    }
		  if(N[1]<N1s[k])
		    {
		      subterm = -(v[1]/(v[k]-v[k-1]))*K[1][k]*(N[1]-N1s[k])*N[k];//Loss of k due to evaporation of monomers from k.
		      t3= t3 + K[1][k]*(-N[1]+N1s[k])*N[k];//Gain of monomers that have evaporated.
		    }
		  if(N[1]>N1s[k])
		    {
		      subterm = (v[1]/(v[k+1]-v[k]))*K[1][k]*(N[1]-N1s[k])*N[k];//Loss of k due to condensation of monomers on k.	      
		      t4 = t4 + K[1][k]*(N[1]-N1s[k])*N[k];//Loss of monomers that have condensed.
		    }
		}
	      // Terms due to coagulation.
	      for(i=2;i<=MAX-1;i++)
		{
		  sum2=sum2+K[k][i]*N[i];
		  for(j=2;j<=k;j++)
		    sum1=sum1+X[i][j][k]*K[i][j]*N[i]*N[j];
		}
	      // Change in particle size due to nucleation, coagulation and surface growth.
	      N[k]= N[k] + step*(0.5*sum1 - N[k]*sum2 + Jk*zeta[k]+ addterm - subterm);//PSD when all the three phenomena are activated.
	    }
	  // Change in the concentration of monomers as they are used up in nucleation and surface growth.
	  N[1] = N[1] - Jk*kstar*step - step*(t1 +t4 - t3-t2);//Monomer balance equation.
	  T=T-step*coolrate;//Reducing the temperature at the given cooling rate.
	  t=t+step;// Time increment.
	}
    }

  //*****************************************
  //         PURE SURFACE GROWTH
  //*****************************************
  
  else if (choice ==4)
    { 	  
      for(i=1;i<=MAX-1;i++)
	{
	  N[i]=0.0;
	  N1s[i]=0.0;
	}
      v[1]=v1;

	  std::cout << "Node 1: " << v[1] << std::endl;

	  for (i = 2; i <= MAX; i++)
	  {
		  v[i] = v[1] * pow(q, (i - 1));
		  
		  if (i <= MAX - 1)
		  {
			  std::cout << "Node: " << i << " Vol: " << v[i] << std::endl;
		  }
	  }

	  std::cout << "inital psd created" << std::endl;

      if((fptr2=fopen("grow.inp","r+"))==NULL)//Opening the input file for surface growth problems
	{
	  printf("Warning: Could not open file");
	  return(-1);
	}
      //Scanning the initial size and number concentration of particles whose surface growth is to be studied.
	  fscanf_s(fptr2,"Enter the node in which you want to put in particles:%d Enter the particle number concentration in #/m3 initially:%lf",&i,&temp);
      fclose(fptr2);
      N[i]=temp;
      t=0.0;
      step=1e-8;
      S=1.001;
      temp1 = 6*v1/pi;
      s1 = pi*pow(temp1,0.66666666667);
      m1=rho*v1;
      Ps=exp(C-D/T)*101325.0*P;
      ns=Ps/(kb*T);
      N[1]=S*ns;
      counter=0;
	  while(t < time_end)
     // while(T>300)
	{
	  sigma=(A-B*T)*1e-3;
	  Ps=exp(C-D/T)*101325.0*P;
	  ns=Ps/(kb*T);
	  S=N[1]/ns;
	  theta = (s1*sigma)/(kb*T);
	  a=(2*sigma)/(pi*m1);
	  b= theta - (4*pow(theta,3))/(27*pow(log(S),2));
	  Jk = pow(ns,2)*S*v1*pow(a,0.5)*exp(b); 
	  c = 0.6666667*theta/log(S);
	  kstar = pow(c,3.0);
	 if(counter==10000)//printing after every 1000 steps
	    { 
	     printf("\n\ntime=%es\n",t);
	  //   for(i=1;i<=MAX-2;i++)
	//	printf("\nN[%d]\t%e",i,N[i]);
	    counter=0;
	   }
	  counter=counter +1;
	      //collision frequency functions for the monomers
	      for(j=1;j<=MAX-1;j++)
		{
		  temp1=1/v[1] + 1/v[j];
		  temp2=pow(v[1],0.333333) + pow(v[j],0.333333);
		  K[1][j]= pow(3.0/(4.0*pi),0.1666666)*pow((6.0*kb*T/rho),0.5)*pow(temp1,0.5)*pow(temp2,2.0);	     
		}
	      //calculation of saturation concentration of monomers over an i sized particle
	      for(i=1;i<=MAX-1;i++)
		{
		  dp[i]=pow((6*v[i]/pi),0.333333);
		  N1s[i]=ns*exp(4*sigma*MW/(R*T*rho*dp[i]));
		}
	      for(k=2;k<=MAX-2;k++)
		{
		  if(step!=1e-4)
		    {
		      if(N[1]>N1s[k-1])
			{
			  if (k==2)
			    addterm=0.0;
			  else
			    addterm = (v[1]/(v[k]-v[k-1]))*K[1][k-1]*(N[1]-N1s[k-1])*N[k-1];//growth of k due to condensation of monomers on k-1
			}		
		      if(N[1]<N1s[k+1])
			addterm = -(v[1]/(v[k+1]-v[k]))*K[1][k+1]*(N[1]-N1s[k+1])*N[k+1];//growth of k due to evaporation of monomers from k+1


		      if(N[1]<N1s[k])
			subterm = -(v[1]/(v[k]-v[k-1]))*K[1][k]*(N[1]-N1s[k])*N[k];//loss of k due to evaporation of monomers from k

		      if(N[1]>N1s[k])
			subterm = (v[1]/(v[k+1]-v[k]))*K[1][k]*(N[1]-N1s[k])*N[k];//loss of k due to condensation of monomers on k	      
		    }
		  N[k]= N[k] + step*(addterm - subterm);//changing PSD due to surface growth/evaporation only
		}
	      T=T-step*coolrate;
	      t=t+step;
	}      
    } 

	for (i = 1; i <= MAX - 1; i++) {

		printf("\n%e\t%e", v[i], N[i]);
		resultPSD << i << "\t" << v[i] << "\t" << N[i] << std::endl;
	}
	std::cout << "results printed" << std::endl;

	resultPSD.close();
	return 0;
}

