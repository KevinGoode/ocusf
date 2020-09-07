#define TRUE (1)
#define FALSE (0)
const double G=9.806; //Acceleration due to gravity m3/s
enum XAREA{CIRC=0,RECT,TRAP};
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


void CalculateSteadyFlowDepth(ChannelInitialConditions initConditions, ChannelShape channel, double Q,double *A, double *P,double *TW,double *SteadyDepth,double *SteadyFlow); //a.k.a. GOSUB 2060 in original code
void CalculateSectionFlowParameters(ChannelShape channel, double YY,double *A, double *P,double *TW);//a.k.a. GOSUB 1860 in original code
void Calculate_Q_C_SF(ChannelShape channel, double VV,double A,double P,double TW,double *QCur,double *CCur,double *SFCur);//a.k.a. GOSUB 2010 in original code
