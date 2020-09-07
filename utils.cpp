
#include <math.h>
#include <stdio.h>
#include "utils.hpp"

//a.k.a GOSUB 2060
void CalculateSteadyFlowDepth(ChannelInitialConditions initConditions, ChannelShape channel, double QCur, double *A, double *P,double *TW,double *YN, double *VN)
{
	double FSR=pow(channel.BedSlope,0.5);
	double HII=40.0,LOO=0.001,RH=0.0,FSH=0.0,WW=0.0,Z=0.0;
	bool bNotConverged=TRUE;
	if(channel.CrossSectionShape==CIRC) HII=channel.Diameter;
	double YY=0.0;
	while(bNotConverged)
	{
		YY=(HII+LOO)/2.0;
		CalculateSectionFlowParameters(channel, YY,A, P,TW);
		RH=*A/(*P);
		FSH=channel.ManningsValue*QCur/((*A)*pow(RH,0.67));
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
	*VN=initConditions.Flow/(*A);
	return;
}
//a.k.a. GOSUB 1860
void CalculateSectionFlowParameters(ChannelShape channel, double YY,double *A, double *P,double *TW)
{   
	/*
	P is wetted perimeter (IE water perimeter minus surface exposed to air)
	TW is width of water exposed to air ( Hence P + TW = Perimeter of water cross section)
	*/
	double D=channel.Diameter;
	double B=channel.Width;
    double FI=channel.SideWallAngle*3.142/180.0;//Side wall angle in radians

	switch(channel.CrossSectionShape)
	{
		case CIRC:
		{   
			// EST is angle in radians from vertical of line joining centre of circle and water surface 
			// touching perimter of circle. IE wetted perimeter, P = Circum *2*EST/2Pi = 2Pi.r*EST/Pi = D*EST
			// and top water width TW = 2*r*sin(EST) = D*sin(EST) 
			double HI=3.1416,LO=0.000,EST=0.0,XR=0.0,Z=0.0;
			bool bNotConverged=TRUE;
			while(bNotConverged)
			{
				EST =(HI+LO)/2.0;
				XR=1.0-(2.0*YY)/channel.Diameter-cos(EST);
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
            *A=0.25*D*D*(EST-0.5*sin(2.0*EST)); // This looks wrong- TODO CHECK
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
void Calculate_Q_C_SF(ChannelShape channel, double VV,double A,double P,double TW, double *QCur,double *CCur,double *SFCur)
{   double MN=channel.ManningsValue;
	*QCur=VV*A;
	*CCur=pow((G*A/TW),0.5);
	*SFCur=MN*MN*pow((P/A),1.33333)*VV*fabs(VV);
}
