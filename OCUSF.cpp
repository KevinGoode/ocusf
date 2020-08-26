#include <stdio.h>
#include <math.h>
/*
	This program converted from BASIC to 'C' by K.Goode from original model in 
	"Water and Wasterwater Engineering Hydraulics"

	NB. Conversion sacrifices programminhg elegance for consistency with original code eg Using many global variables
	NB. However, all gotos removed .
*/
#define TRUE (1)
#define FALSE (0)
#define SIZE_OF_LINE (64)
#define SIZE_OF_ARRAY (100)
#define SMALL_NUM (0.000001)
char line[SIZE_OF_LINE];
enum XAREA{CIRC=0,RECT,TRAP};
double Y[SIZE_OF_ARRAY],V[SIZE_OF_ARRAY],YP[SIZE_OF_ARRAY],VP[SIZE_OF_ARRAY],C[SIZE_OF_ARRAY],Q[SIZE_OF_ARRAY],SF[SIZE_OF_ARRAY];
//WB=StoredAtStart+Qin-StoredAtEnd-Qout;
double FlowIn=0.0;
double CurrentStored=0.0;
double CurrentWaterBalance=0.0;
double CumulativeWaterBalance=0.0;
void InitWaterBalance(double dx);
double CalculateVolumeStored(double dx);
double CalculateWaterBalance(double dx, double dt);
double CalculateNetFlowIn();
struct ChannelShape
{
	double Diameter;//D (m)
	double Width;//B (m)
	double SideWallAngle;//FI (deg)
	double Length;//XL (m)
	double ManningsValue;//MN
	double BedSlope;//S0 (sin theta);
	XAREA  CrossSectionShape;
};
enum BoundaryConditionVariable{FLOW=1,DEPTH};
struct BoundaryConditions
{
	BoundaryConditionVariable Variable;
	double RateOfChange;
	double FinalValue;
};
struct ChannelBoundaryConditions
{
	BoundaryConditions UpStream;
	BoundaryConditions DownStream;
};
struct ChannelInitialConditions
{
	double Flow; //(m3/s)
	double USDepth; //(m)
};
struct ChannelShape Channel;
struct ChannelBoundaryConditions BoundaryConditions;
struct ChannelInitialConditions InitConditions;

int N=0;//Number of reaches into which the channel length is divided.
int NITER=0;//Number of computational iterations
double G=9.806; //Acceleration due to gravity m3/s
int NS=0;//Number of computational nodes (N+1)
double DX=0.0;//Distance step (m)
double DT=0.0;//Time step (s)
double T=0.0; //Current time (s)
double TM=0.0; //Current time (mins)
double DXX=0.0; //(m)
double DXI=0.0; //(m)
double ZETA=0.0; //(dimensionless)


//GOSUB routines
void CalculateSteadyFlowDepth(double Q,double *A, double *P,double *TW,double *SteadyDepth,double *SteadyFlow); //a.k.a. GOSUB 2060 in original code
void CalculateSectionFlowParameters(double YY,double *A, double *P,double *TW);//a.k.a. GOSUB 1860 in original code
void Calculate_Q_C_SF(double VV,double A,double P,double TW,double *QCur,double *CCur,double *SFCur);//a.k.a. GOSUB 2010 in original code
void UStreamBoundaryFinalFlowCheck(double *QUS);//a.k.a. GOSUB 2210 in original code
void UStreamBoundaryFinalDepthCheck();//a.k.a. GOSUB 2270 in original code
void DStreamBoundaryFinalFlowCheck(double *QDS);//a.k.a. GOSUB 2340 in original code
void DStreamBoundaryFinalDepthCheck();//a.k.a. GOSUB 2400 in original code

//New routines
void ReadInputParameters(int numArgs, char **argv);
void Initialise(double *A, double *P,double *TW);
void InitialiseFlowAndDepth(double *A, double *P,double *TW,double SteadyDepth,double SteadyFlow);
void OutputFlowAndDepthHeaders();
void OutputFlowAndDepthRecords();
void CheckBoundaryConditions(double *QUS,double *QDS);
//Extra utilitiy methods
void ArraySet(double *pArray,double value,int size);
/* PROGRAM OCUSF*/
int main(int numArgs, char **argv)
{   
	double QUS=0.0;//upstream q
	double QDS=0.0;//Downstream q
	double SFCur=0.0;//Steady SF in steady state (NB in original code this also is called SF which conflict with array of these values)
	double CCur=0.0;///Steady C in steady state (NB in original code this also is called C which conflict with array of these values)
	double QCur=0.0;//Current Q (m3/s)
    double A=0.0;//Area (m2)
    double TW=0.0;
    double P=0.0;//Wetted Perimeter (m)

	/*Option of passing data from input file*/
	ReadInputParameters(numArgs, argv);

	G=9.809;NS=N+1;
	
	Initialise(&A, &P,&TW);
	DX=Channel.Length/(double)N;
	DT=DX/(V[NS-1]+C[NS-1]);
	T=0.0;
    
	OutputFlowAndDepthHeaders();
	OutputFlowAndDepthRecords();
	
    int COUNT=0;
    for(COUNT=0;COUNT<NITER;COUNT++)
	{
		CheckBoundaryConditions(&QUS,&QDS);

		TM=T/60.0;
		DXX=0.0;
		int j=0;
		for(j=0;j<NS;j++)
		{
			DXI=(fabs(V[j])+C[j])*DT;
			if(DXI>DXX) DXX=DXI;
		}
		ZETA=DXX/DX;
		DT=DT/ZETA;
		T=T+DT;TM=T/60.0;
		double TH=0.0;
		double CA=0.0,VR=0.0,CR=0.0,YR=0.0,SR=0.0;
		double CB=0.0,VS=0.0,CS=0.0,YS=0.0,SS=0.0;
		TH=DT/DX;

         InitWaterBalance(DX);
        
		//Interior points
		int k=0;
		for(k=1;k<N;k++)
		{
			CA=C[k]-C[k-1];
			VR=(V[k]+TH*(C[k]*V[k-1]-V[k]*C[k-1]))/(1.00+TH*(V[k]-V[k-1]+CA));
			CR=(C[k]-VR*TH*CA)/(1.0+TH*CA);
			YR=Y[k]-TH*(VR+CR)*(Y[k]-Y[k-1]);
			CalculateSectionFlowParameters(YR,&A, &P,&TW);
			Calculate_Q_C_SF(VR,A,P,TW,&QCur,&CCur,&SFCur);
			CR=CCur;SR=SFCur;

			CB=C[k]-C[k+1];
			VS=(V[k]-TH*(V[k]*C[k+1]-C[k]*V[k+1]))/(1.00-TH*(V[k]-V[k+1]-CB));
			CS=(C[k]+VS*TH*CB)/(1.0+TH*CB);
			YS=Y[k]+TH*(VS-CS)*(Y[k]-Y[k+1]);
			CalculateSectionFlowParameters(YS,&A, &P,&TW);
			Calculate_Q_C_SF(VS,A,P,TW,&QCur,&CCur,&SFCur);
			CS=CCur;SS=SFCur;

			YP[k]=(YS*CR+YR*CS+CR*CS*((VR-VS)/G-DT*(SR-SS)))/(CR+CS);
			VP[k]=VR-G*((YP[k]-YR)/CR+DT*(SR-Channel.BedSlope));
		}


		//Upstream boundary conditions
		CB=C[0]-C[1];
		VS=(V[0]-TH*(V[0]*C[1]-C[0]*V[1]))/(1.0-TH*(V[0]-V[1]-CB));
		CS=(C[0]+VS*TH*CB)/(1.0+TH*CB);
		YS=Y[0]+TH*(VS-CS)*(Y[0]-Y[1]);
        CalculateSectionFlowParameters(YS,&A, &P,&TW);
		Calculate_Q_C_SF(VS,A, P,TW,&QCur,&CCur,&SFCur);
		double C2=0.0,CM=0.0;
		double FY=0.0,FDY=0.0,DELY=0.0;
		C2=G/CCur;
		CM=VS-C2*YS-G*DT*(SFCur-Channel.BedSlope);
		if(BoundaryConditions.UpStream.Variable != DEPTH)
		{
			//Specified variation in Q at US boundary
			double YY;
			YY=Y[0];
			bool bNotConverged=TRUE;
			while(bNotConverged)
			{
				CalculateSectionFlowParameters(YY,&A, &P,&TW);
				FY=QUS/A-C2*YY-CM;
				FDY=-(QUS/pow(A,2))*TW-C2;
				DELY=-FY/FDY;
				YY+=DELY;
				bNotConverged=(fabs(DELY)>0.001);
			}
			YP[0]=YY;VP[0]=C2*YP[0]+CM;
		}
		//Specified variation in Y at US boundary
		VP[0]=C2*YP[0]+CM;
		//Downstream boundary
		CA=C[NS-1]-C[N-1];
		VR=(V[NS-1]+TH*(C[NS-1]*V[N-1]-V[NS-1]*C[N-1]))/(1.0+TH*(V[NS-1]-V[N-1]+CA));
		CR=(C[NS-1]-VR*TH*CA)/(1.0+TH*CA);
		YR=Y[NS-1]-TH*(VR+CR)*(Y[NS-1]-Y[N-1]);
		CalculateSectionFlowParameters(YR,&A, &P,&TW);
		Calculate_Q_C_SF(VR,A, P,TW,&QCur,&CCur,&SFCur);
		double C4=G/CCur;
		double CP=VR+C4*YR-G*DT*(SFCur-Channel.BedSlope);
		if( BoundaryConditions.DownStream.Variable != DEPTH)
		{
			//Specified variation in DS Discharge Q
			double YY=Y[NS-1];
			bool bNotConverged=TRUE;
			while(bNotConverged)
			{
				CalculateSectionFlowParameters(YY,&A, &P,&TW);
				FY=QDS/A+C4*YY-CP;
				FDY=-(QDS/pow(A,2))*TW+C4;
				DELY=-FY/FDY;
				YY+=DELY;
				bNotConverged=(fabs(DELY)>0.001);
			}
			YP[NS-1]=YY;
		}
		//Specified variation in downstream Y
		VP[NS-1]=CP-C4*YP[NS-1];
		//Update variable values for current time step
		int m=0;
		for(m=0;m<NS;m++)
		{
			CalculateSectionFlowParameters(YP[m],&A,&P,&TW);
			Calculate_Q_C_SF(VP[m],A,P,TW,&QCur,&CCur,&SFCur);
			V[m]=VP[m];
			Q[m]=V[m]*A;
			C[m]=CCur;
			Y[m]=YP[m];
			SF[m]=SFCur;
		}
	
        CurrentWaterBalance=CalculateWaterBalance(DX, DT);
        CumulativeWaterBalance+=CurrentWaterBalance;

        if(fabs(ceil((double)COUNT/2.0) -floor((double)COUNT/2.0))>SMALL_NUM)
		{
			//Even COUNT numbers only
			OutputFlowAndDepthRecords();
		}
	
	}

	return 0;
}
//a.k.a GOSUB 2060
void CalculateSteadyFlowDepth( double QCur, double *A, double *P,double *TW,double *YN, double *VN)
{
	double FSR=pow(Channel.BedSlope,0.5);
	double HII=40.0,LOO=0.001,RH=0.0,FSH=0.0,WW=0.0,Z=0.0;
	bool bNotConverged=TRUE;
	if(Channel.CrossSectionShape==CIRC) HII=Channel.Diameter;
	double YY=0.0;
	while(bNotConverged)
	{
		YY=(HII+LOO)/2.0;
		CalculateSectionFlowParameters(YY,A, P,TW);
		RH=*A/(*P);
		FSH=Channel.ManningsValue*QCur/((*A)*pow(RH,0.67));
		WW=FSH-FSR;
		if(WW >0)
		{
			LOO=YY;
		}
		else
		{
			HII=YY;
		}
		Z=(HII+LOO)/2.0;
		bNotConverged=(fabs(Z-YY)>.0002);
	}
	*YN=YY;
	*VN=InitConditions.Flow/(*A);
	return;
}
//a.k.a. GOSUB 1860
void CalculateSectionFlowParameters(double YY,double *A, double *P,double *TW)
{
	double D=Channel.Diameter;
	double B=Channel.Width;
    double FI=Channel.SideWallAngle*3.142/180.0;//Side wall angle in radians

	switch(Channel.CrossSectionShape)
	{
		case CIRC:
		{
			double HI=3.1416,LO=0.000,EST=0.0,XR=0.0,Z=0.0;
			bool bNotConverged=TRUE;
			while(bNotConverged)
			{
				EST =(HI+LO)/2.0;
				XR=1.0-(2.0*YY)/Channel.Diameter-cos(EST);
				if(XR<0)
				{
					LO=EST;
				}
				else
				{
					HI=EST;
				}
				Z=(HI+LO)/2.0;
				bNotConverged=(fabs(Z-EST)>0.001);
			}
			*P=D*EST;
            *A=0.25*D*D*(EST-0.5*sin(2.0*EST));
			*TW=D*sin(EST);
			break;
		}
		case RECT:
		{
			*A=B*YY;
            *P=(B+2.0*YY);
            *TW=B;
			break;
		}
		case TRAP:
		{
			*A=YY*(B+YY/tan(FI));
			*P=B+2*YY/(sin(FI));
            *TW=B+2.0*YY/tan(FI);
			break;
		}
		default:
		{
			printf("Error\n");
			break;
		}
	}
	return;
}
//a.k.a. GOSUB 2010
void Calculate_Q_C_SF(double VV,double A,double P,double TW, double *QCur,double *CCur,double *SFCur)
{   double MN=Channel.ManningsValue;
	*QCur=VV*A;
	*CCur=pow((G*A/TW),0.5);
	*SFCur=MN*MN*pow((P/A),1.33333)*VV*fabs(VV);
}
//TODO both flow and depth cheacks can be consolidated from 4 routines into 2
//a.ka. GOSUB 2210
void UStreamBoundaryFinalFlowCheck(double *QUS)
{
	if(BoundaryConditions.UpStream.RateOfChange > 0.0)
	{
		if(*QUS>=BoundaryConditions.UpStream.FinalValue)
		{
			*QUS=BoundaryConditions.UpStream.FinalValue;
		}
	}
	else
	{
		if(*QUS<=BoundaryConditions.UpStream.FinalValue)
		{
			*QUS=BoundaryConditions.UpStream.FinalValue;
		}
	}
}
//a.ka. GOSUB 2270
void UStreamBoundaryFinalDepthCheck()
{
	if(fabs(BoundaryConditions.UpStream.RateOfChange) <SMALL_NUM)
	{
		YP[0]=Y[0];
	}
	else if(BoundaryConditions.UpStream.RateOfChange >0)
	{
		if(YP[0]>=BoundaryConditions.UpStream.FinalValue)
		{
			YP[0]=BoundaryConditions.UpStream.FinalValue;
		}
	}
	else
	{
		if(YP[0]<=BoundaryConditions.UpStream.FinalValue)
		{
			YP[0]=BoundaryConditions.UpStream.FinalValue;
		}
	}
}
//a.ka. GOSUB 2340
void DStreamBoundaryFinalFlowCheck(double *QDS)
{
	if(BoundaryConditions.DownStream.RateOfChange> 0.0)
	{
		if(*QDS>=BoundaryConditions.DownStream.FinalValue)
		{
			*QDS=BoundaryConditions.DownStream.FinalValue;
		}
	}
	else
	{
		if(*QDS<=BoundaryConditions.DownStream.FinalValue)
		{
			*QDS=BoundaryConditions.DownStream.FinalValue;
		}
	}
}
//a.ka. GOSUB 2400
void DStreamBoundaryFinalDepthCheck()
{
	if(fabs(BoundaryConditions.DownStream.RateOfChange) <SMALL_NUM)
	{
		YP[NS-1]=Y[0];
	}
	else if(BoundaryConditions.DownStream.RateOfChange>0)
	{
		if(YP[NS-1]>=BoundaryConditions.DownStream.FinalValue)
		{
			YP[NS-1]=BoundaryConditions.DownStream.FinalValue;
		}
	}
	else
	{
		if(YP[NS-1]<=BoundaryConditions.DownStream.FinalValue)
		{
			YP[NS-1]=BoundaryConditions.DownStream.FinalValue;
		}
	}
}
void ArraySet(double *pArray,double value,int size)
{
	int i=0;
	for(i=0;i<size;i++)
	{
		pArray[i]=value;
	}
}
void ReadInputParameters(int numArgs, char **argv)
{
	int UN=0; //Upstream boundary condition number  1.) Linear variation of discharge with time 2.) Linear variation of depth with time
	int DN=0; //Downstream boundary condition number  1.) Linear variation of discharge with time 2.) Linear variation of depth with time
	FILE* fInputFile=NULL;
	bool bInputFile=FALSE;
	int NO=0;
	if(numArgs == 2)
	{
		fInputFile =fopen(argv[1],"r");
		if(fInputFile > 0)bInputFile=TRUE;
	}
	if(!bInputFile) fInputFile=stdin;
	printf("---------------------------PROGRAM OCUSF---------------------------\n");
	printf("This program computes the transient flow and water depth in\n");
	printf("open channnels of rectangular, trapezoidal and circular cross-\n");
	printf("sections, using a numerical computation procedure based on the\n");
	printf("method of characterisitics. The computation of frictaional resitance\n");
	printf("is based on the Manning equation.\n");
	printf("\n");
	printf("Computation starts from a specified steady state at time zero.\n");
	printf("The program offers a choice of two initial steady states viz.\n");
	printf("steady uniform flow and zero flow\n");
	printf("\n");
	printf("Boundary conditions: The program caters for the following\n");
	printf("parameter variations at both ends of the channel:\n");
	printf(" a.) linear variation of flow depth with time\n");
	printf(" b.) linear variation of discharge rate with time\n");
	printf("\n");
	printf("A constant value for either boundary paramter is obtained\n");
	printf("by specifying a zero rate for the parameter variation\n");
	printf("\n");
	printf("---------------------PRESS THE SPACE BAR TO CONTINUE---------------\n");
	fgets(line,SIZE_OF_LINE,fInputFile);
	printf("Enter Channel Data:\n");
	printf("1 CIRCULAR 2 RECTANGULAR 3 TRAPEZOIDAL:"); fscanf(fInputFile,"%d",&NO);
	if(NO ==1){ Channel.CrossSectionShape=CIRC;printf("Diameter (m):"); fscanf(fInputFile,"%lf",&Channel.Diameter);}
	if(NO ==2){ Channel.CrossSectionShape=RECT;printf("Channel width (m):"); fscanf(fInputFile,"%lf",&Channel.Width);}
	if(NO ==3){ Channel.CrossSectionShape=TRAP;printf("Bottom width (m):"); fscanf(fInputFile,"%lf",&Channel.Width);}
	if(NO ==3){ printf("Angle of side wall to horl (deg):"); fscanf(fInputFile,"%lf",&Channel.SideWallAngle);}
	printf("Enter Channel Length (m):"); fscanf(fInputFile,"%lf",&Channel.Length);
	printf("Enter Mannings n-value:"); fscanf(fInputFile,"%lf",&Channel.ManningsValue);
	printf("Enter Channel Bed Slope(sin theta):"); fscanf(fInputFile,"%lf",&Channel.BedSlope);
	printf("Enter Initial steadyflow (m3/s):");fscanf(fInputFile,"%lf",&InitConditions.Flow);
	if(fabs(InitConditions.Flow-0.0)<SMALL_NUM){printf("Enter water depth at upstream end of channel(m):");fscanf(fInputFile,"%lf",&InitConditions.USDepth); }
    printf("Select upstream boundary condition\n");
	printf(" 1.) Linear variation of discharge with time\n");
	printf(" 2.) Linear variation of depth with time\n");
	printf("Enter 1 or 2 as appropriate:");fscanf(fInputFile,"%d",&UN);
	BoundaryConditions.UpStream.Variable=(BoundaryConditionVariable)UN;
	if(UN==1){printf("Enter rate of discharge variations (m3/s/s):");fscanf(fInputFile,"%lf",&BoundaryConditions.UpStream.RateOfChange);}
    if(UN==1){printf("Enter final upstream discharge rate (m3/s):");fscanf(fInputFile,"%lf",&BoundaryConditions.UpStream.FinalValue);}
    if(UN==2){printf("Enter rate of depth variation with time (m/s):");fscanf(fInputFile,"%lf",&BoundaryConditions.UpStream.RateOfChange);}
	if(!(fabs(BoundaryConditions.UpStream.RateOfChange-0.0)<SMALL_NUM))
	{
		if(UN==2){printf("Enter final upstream depth (m):");fscanf(fInputFile,"%lf",&BoundaryConditions.UpStream.FinalValue);}
	}
	
	printf("Select downstream boundary condition\n");
	printf(" 1.) Linear variation of discharge with time\n");
	printf(" 2.) Linear variation of depth with time\n");
	printf("Enter 1 or 2 as appropriate:");fscanf(fInputFile,"%d",&DN);
	BoundaryConditions.DownStream.Variable=(BoundaryConditionVariable)DN;
	if(DN==1){printf("Enter rate of discharge variations (m3/s/s):");fscanf(fInputFile,"%lf",&BoundaryConditions.DownStream.RateOfChange);}
    if(DN==1){printf("Enter final downstream discharge rate (m3/s):");fscanf(fInputFile,"%lf",&BoundaryConditions.DownStream.FinalValue);}
    if(DN==2){printf("Enter rate of depth variation with time (m/s):");fscanf(fInputFile,"%lf",BoundaryConditions.DownStream.RateOfChange);}
	if(!(fabs(BoundaryConditions.DownStream.RateOfChange-0.0)<SMALL_NUM))
	{
		if(DN==2){printf("Enter final downstream depth (m):");fscanf(fInputFile,"%lf",&BoundaryConditions.DownStream.FinalValue);}
	}
	printf("\n");
	printf("Enter number of reaches into which the channel length is divided.\n");
	printf("for computational purposes (multiple of 10):");fscanf(fInputFile,"%d",&N);
	printf("Enter number of computational iterations");fscanf(fInputFile,"%d",&NITER);
    printf("\n");
	printf("DATA INPUT COMPLETED; COMPUTATION IN PROGRESS.\n");
	
	if(bInputFile==TRUE) fclose(fInputFile);
}
void OutputFlowAndDepthHeaders()
{
	printf("\n");
	printf("-------------------TABULATION OF COMPUTED VALUES FOLLOWS----------\n");
	printf("-TIME---------------------DISTANCE ALONG CHANNEL------------------\n");
	printf("-(MIN)---0.0L     0.2L     0.4L     0.6L     0.8L     1.0L--------\n");
}
void OutputFlowAndDepthRecords()
{
	int i=0;
	printf(" %02.2lf ",TM);
	for(i=0;i<N;i+=N/5) printf("  %02.3lf  ",Q[i]); printf("  %02.3lf  ",Q[NS-1]);
	printf("   Q (m^3/s)\n      ");
	for(i=0;i<N;i+=N/5) printf("  %02.3lf  ",Y[i]);
	printf("  %02.3lf  [%02.3lf  %02.3lf  ]",Y[NS-1],CurrentWaterBalance,CumulativeWaterBalance);printf("   Depth (m)[Balance (m3), CumBalance(m3)]\n");
}
void CheckBoundaryConditions(double *QUS,double *QDS)
{
		double USRate=BoundaryConditions.UpStream.RateOfChange;
		double DSRate=BoundaryConditions.DownStream.RateOfChange;
	    if(BoundaryConditions.UpStream.Variable==FLOW){*QUS=Q[0]+DT*USRate;UStreamBoundaryFinalFlowCheck(QUS);}
		if(BoundaryConditions.UpStream.Variable==DEPTH){YP[0]=Y[0]+DT*USRate;UStreamBoundaryFinalDepthCheck();}
		if(BoundaryConditions.DownStream.Variable==FLOW){*QDS=Q[NS-1]+DT*DSRate;DStreamBoundaryFinalFlowCheck(QDS);}
		if(BoundaryConditions.DownStream.Variable==DEPTH){YP[NS-1]=Y[NS-1]+DT*DSRate;DStreamBoundaryFinalDepthCheck();}
}
void Initialise(double *A, double *P,double *TW)
{
	/*Computation of steady flow depth YN and VELOCITY VN*/
	double SteadyVelocity=0.0; //(m/s)
	double SteadyDepth=0.0; //(m)
	if (fabs(InitConditions.Flow-0.0)>SMALL_NUM)
	{ 
		double QCur=fabs(InitConditions.Flow);
		CalculateSteadyFlowDepth(QCur,A,P,TW,&SteadyDepth,&SteadyVelocity);
	}
    //Assignment of initial values to variables

	InitialiseFlowAndDepth(A,P,TW,SteadyDepth,SteadyVelocity);
}
void InitialiseFlowAndDepth(double *A, double *P,double *TW,double SteadyDepth, double SteadyVelocity)
{
	double XL=Channel.Length;
	double QCur=0,CCur=0.0,SFCur=0.0;
	int i=0;
	for (i=0;i<NS;i++)
	{
		if(fabs(InitConditions.Flow-0.0)<SMALL_NUM) 
		{
			V[i]=0;
			Y[i]=InitConditions.USDepth+((double)(i+1)-1.0)*XL/(double)N*Channel.BedSlope;
		}
		else
		{
			V[i]=SteadyVelocity;
			Y[i]=SteadyDepth;
		}
		CalculateSectionFlowParameters(Y[i],A,P,TW);
		Calculate_Q_C_SF(V[i],*A,*P,*TW,&QCur,&CCur,&SFCur);
		SF[i]=SFCur;
		Q[i]=InitConditions.Flow;
		C[i]=CCur;
	}
}
void InitWaterBalance(double dx)
{
    CurrentStored =CalculateVolumeStored( dx);
    FlowIn=CalculateNetFlowIn();
}
double CalculateWaterBalance(double dx, double dt)
{
    double stored=CalculateVolumeStored(dx);
    double deltaStored=stored-CurrentStored;
    double wb=deltaStored-FlowIn*dt;
    return wb;
}
double CalculateNetFlowIn()
{
    return Q[0]-Q[NS-1];
}
double CalculateVolumeStored(double dx)
{
    double stored=0.0;
    int i=0;
	for(i=1;i<N;i++)
    {
        double A,P,TW;
        CalculateSectionFlowParameters((Y[i]+Y[i+1])/2.0,&A, &P,&TW);
        stored+=dx*A;
    }
    return stored;
}