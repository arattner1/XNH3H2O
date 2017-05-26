/*******************************************************
NH3H2O.h									02/28/2011
Alexander Rattner							STSL
Function definitions for the thermophysical properties
of Ammonia-Water fluid mixtures.
*******************************************************/

#include "NH3H2O.h"


//==================================================================================================
//TOP LEVEL FUNCTION BELOW
//==================================================================================================
void SetProps(char *Mode, Propset *P)
{
	//TPX mode
	if ( strcmp(Mode, "TPX") == 0 )
	{
		//Set initial guess values for x and y
		P->x	= 0.01;
		P->y	= 0.99;
		//Get missing properties
		Qxy_TPX ( &(P->Q), &(P->x), &(P->y), P->T, P->P, P->psi );
		//Ensure quality is in allowed range
		P->Q	= (P->Q > 1.0) ? 1.0 : P->Q;
		P->Q	= (P->Q < 0.0) ? 0.0 : P->Q;
		//Now get the rest of the properties
		P->MW	= MW_m( P->psi );
		P->rho	= 1/( (P->Q)*V_vm(P->T, P->P, P->y) + (1 - P->Q)*V_lm(P->T, P->P, P->x) );
		//Get cp, only exists for single phase
		if ( P->Q == 0)
		{	P->cp = Cp_lm( P->T, P->P, P->x ); P->x = P->psi; P->y = 0;}
		else if (P->Q == 1)
		{	P->cp = Cp_vm( P->T, P->P, P->y ); P->x = 0; P->y = P->psi;}
		//Get h, u, s
		P->h	= (P->Q)*H_vm(P->T, P->P, P->y) + (1 - P->Q)*H_lm(P->T, P->P, P->x);
		P->u	= (P->Q)*U_vm(P->T, P->P, P->y) + (1 - P->Q)*U_lm(P->T, P->P, P->x);
		P->s	= (P->Q)*S_vm(P->T, P->P, P->y) + (1 - P->Q)*S_lm(P->T, P->P, P->x);
	}
}

//Like above, but searches from a tree for initial guesses
void SetPropsSearch(char *Mode, Propset *P, Trees *TS, Proplib *PL)
{	
	//Flag for if a tree was searched
	int searched = 0;

	//TPX mode
	if      ( strcmp(Mode, "TPX") == 0 )
	{	searched = getTPX(P, TS->TPX); }
	//PXH mode
	else if ( strcmp(Mode, "PXH") == 0 )
	{	searched = getPXH(P, TS->PXH); }
	//TPQ mode
	else if ( strcmp(Mode, "TPQ") == 0 )
	{	searched = getTPQ(P, TS->TPQ); }
	//PXQ mode
	else if ( strcmp(Mode, "PXQ") == 0 )
	{	searched = getPXQ(P, TS->PXQ); }
	//XVU mode
	else if ( strcmp(Mode, "XVU") == 0 )
	{	searched = getXVU2(P, TS, PL); }
	//TXV mode
	else if ( strcmp(Mode, "TXV") == 0 )
	{	searched = getTXV(P, TS->TXV); }
    
	//Check if we had to search
	if (searched == 1)		//We searched a tree, should add the new point to the database
	{	AddPropPoint(Mode, P, PL, TS); }
}

//Searches tree to find closest value to input triple
void PropTreeSearch(double a, double b, double c, Propset *P, void *Tree)
{
	//Define some variables
	void		*kdres;
	
	//Get results from the tree out the tree
	kdres = kd_nearest3(Tree, a, b, c);
	*P = *( (Propset*) kd_res_item_data(kdres) );
	
	//bound the concentrations and quality
	P->x = median( 1E-6, P->x, (1 - 1E-6) );
	P->y = median( 1E-6, P->y, (1 - 1E-6) );
	P->psi = median( 1E-6, P->psi, (1 - 1E-6) );
	P->Q = median( 1E-6, P->Q, (1 - 1E-6) );
	
	//clean up
	kd_res_free(kdres);
}

void PrintProps(Propset *P)
{
	mexPrintf("\nTemperature:\t\t\t%10.4lf\t[K]\n", P->T);
	mexPrintf("Pressure:\t\t\t%10.4lf\t[Pa]\n", P->P);
	mexPrintf("Mixture NH3 mass fraction:\t%10.4lf\t[-]\n", P->psi);
	mexPrintf("Enthalpy of mixture:\t\t%10.4lf\t[J/kg]\n", P->h);
	mexPrintf("Energy of mixture:\t\t%10.4lf\t[J/kg]\n", P->u);
	mexPrintf("Entropy of mixture:\t\t%10.4lf\t[J/kg-K]\n", P->s);
	mexPrintf("Specific heat of mixture:\t%10.4lf\t[J/kg-K]\n", P->cp);
	mexPrintf("Density of mixture:\t\t%10.4lf\t[kg/m^3]\n", P->rho);
	mexPrintf("Specific volume of mixture:\t%14.4lf\t[m^3/kg]\n", P->v);
	mexPrintf("Liquid NH3 mass fraction:\t%10.4lf\t[-]\n", P->x);
	mexPrintf("Vapor NH3 mass fraction:\t%10.4lf\t[-]\n", P->y);
	mexPrintf("Quality:\t\t\t%10.4lf\t[-]\n\n", P->Q);

}

//Reads property data from a file, returns the file length
//File should follow the format: TPXQxyhsuv
void LoadPropData (char *fName, Proplib *L)
{
	//Declare variables
	FILE	*fp;			//Pointer to data file
	char	fLine[256];		//line of data from file
	size_t	dataSize;		//Size of a block of data
	Propset	*PD;		//Tree data array
	int		n = 0;			//Number of datapoints read
	int		nTotal;			//Total number of data points
	int		NBlock = 1000;	//Increment size for data point allocation
	
	//initial allocation
	dataSize		= NBlock*sizeof(Propset);
	//*PropDataPtr	= (Propset*)malloc(dataSize);
	//PropData		= *PropDataPtr;
	L->D			= (Propset*)mxMalloc(dataSize);
	mexMakeMemoryPersistent(L->D);
	PD				= L->D;
	nTotal			= NBlock;
	
	//First open the file and read data
	fp = fopen (fName, "rt");
	while ( fgets(fLine, 256, fp) != NULL )
	{

		sscanf(fLine, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &(PD[n].T), &(PD[n].P), &(PD[n].psi), &(PD[n].Q), &(PD[n].x), &(PD[n].y), &(PD[n].h), &(PD[n].s), &(PD[n].u), &(PD[n].v) );
		//Rho is calculated explicitly
		PD[n].rho = 1.0 / PD[n].v;
		n++;
		
		//Check if we need to allocate more memory
		if (n == nTotal)
		{
			nTotal 		+= NBlock;
			dataSize 	= nTotal*sizeof(Propset);
			L->D 		= (Propset*) mxRealloc( L->D, dataSize );
			mexMakeMemoryPersistent(L->D);
			PD			= L->D;
		}
	}

	//File cleanup
	fclose(fp);
	
	//Put counts in property library
	L->n		= n;
	L->N		= nTotal;
	L->NBlock	= NBlock;
}

//Writes property data to a file
//File will follow the format: TPXQxyhsuv
void SavePropData (char *fName, Proplib *L)
{
	//Declare variables
	FILE	*fp;			//Pointer to data file
	Propset	*P;				//Current property set
	int		i;				//counter
	
	//Open file
	fp = fopen (fName, "wt");
	
	//Loop through property data
	P	= &(L->D[0]);
	fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", P->T, P->P, P->psi, P->Q, P->x, P->y, P->h, P->s, P->u, P->v);
	for (i = 1; i < L->n; i++)
	{
		P	= &(L->D[i]);
		fprintf(fp, "\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", P->T, P->P, P->psi, P->Q, P->x, P->y, P->h, P->s, P->u, P->v);
	}
	
	//Close file
	fclose(fp);
}

//Builds a tree of the data in TreeDataPtr, indexed by the mode
void BuildPropTree(char *Mode, Proplib *L, void **Tree)
{
	//Variables
	int i;								//Counter
	Propset *TreeData = L->D;			//Property data array
	
	//now make the appropriate tree
	if (strcmp(Mode, "TPX") == 0 )	//Make a TPX indexed oct tree
	{	
		//Initialize the tree
		*Tree = kd_create(3);
		//Start feeding in data
		for (i = 0; i < L->n; i++)
			kd_insert3(*Tree, TreeData[i].T/300, TreeData[i].P/1E6, TreeData[i].psi, &(TreeData[i]) );
		
			//Un-normalized
			//kd_insert3(*Tree, TreeData[i].T, TreeData[i].P, TreeData[i].psi, &(TreeData[i]) );
	}
	if (strcmp(Mode, "TPQ") == 0 )	//Make a TPQ indexed oct tree
	{	
		//Initialize the tree
		*Tree = kd_create(3);
		//Start feeding in data
		for (i = 0; i < L->n; i++)
			kd_insert3(*Tree, TreeData[i].T/300, TreeData[i].P/1E6, TreeData[i].Q, &(TreeData[i]) );
			
			//Un-normalized
			//kd_insert3(*Tree, TreeData[i].T, TreeData[i].P, TreeData[i].Q, &(TreeData[i]) );
	}
	else if (strcmp(Mode, "PXH") == 0 )	//Make a PXH indexed oct tree
	{	
		//Initialize the tree
		*Tree = kd_create(3);
		//Start feeding in data
		for (i = 0; i < L->n; i++)
			kd_insert3(*Tree, TreeData[i].P/1E6, TreeData[i].psi, TreeData[i].h/1E6, &(TreeData[i]) );
			//Un-normalized
			//kd_insert3(*Tree, TreeData[i].P, TreeData[i].psi, TreeData[i].h, &(TreeData[i]) );
	}
	else if (strcmp(Mode, "PXQ") == 0 )	//Make a PXQ indexed oct tree
	{	
		//Initialize the tree
		*Tree = kd_create(3);
		//Start feeding in data
		for (i = 0; i < L->n; i++)
			kd_insert3(*Tree, TreeData[i].P/1E6, TreeData[i].psi, TreeData[i].Q, &(TreeData[i]) );
			//Un-normalized
			//kd_insert3(*Tree, TreeData[i].P, TreeData[i].psi, TreeData[i].Q, &(TreeData[i]) );
	}
	else if (strcmp(Mode, "XVU") == 0 )	//Make a XVU indexed oct tree
	{	
		//Initialize the tree
		*Tree = kd_create(3);
		//Start feeding in data
		for (i = 0; i < L->n; i++)
			kd_insert3(*Tree, TreeData[i].psi, TreeData[i].v, TreeData[i].u/1E6, &(TreeData[i]) );
			//Un-normalized
			//kd_insert3(*Tree, TreeData[i].psi, TreeData[i].v, TreeData[i].u, &(TreeData[i]) );
	}
	else if (strcmp(Mode, "TXV") == 0 )	//Make a TXV indexed oct tree
	{	
		//Initialize the tree
		*Tree = kd_create(3);
		//Start feeding in data
		for (i = 0; i < L->n; i++)
			kd_insert3(*Tree, TreeData[i].T/300, TreeData[i].psi, TreeData[i].v, &(TreeData[i]) );
	}
}

//Adds a datapoint to a property library
void AddPropPoint(char *Mode, Propset *P, Proplib *L, Trees *TS)
{
	//Define variables
	size_t		szCur;			//Datablock size
	int			i;				//Current index
	Propset		*TreeData;		//Temporary pointer
	Propset		*OldData;		//Pointer to old data
	
	
	//First add the new datapoint to the library
	//Check if we need to allocate more memory
	if (L->n == L->N)
	{
		L->N 		+= L->NBlock;
		szCur 		= L->N*sizeof(Propset);
		OldData		= L->D;
		L->D 		= (Propset*) mxRealloc( L->D, szCur );
		mexMakeMemoryPersistent(L->D);
		TreeData	= L->D;
		
		//we also need to rebuild the trees if the old pointer moved
		if (OldData != TreeData)
		{
			if (TS->hasTPX == 1)
			{
				kd_clear(TS->TPX);
				for (i=0; i < L->n; i++)
					kd_insert3(TS->TPX, TreeData[i].T/300, TreeData[i].P/1E6, TreeData[i].psi, &(TreeData[i]) );
					//kd_insert3(TS->TPX, TreeData[i].T, TreeData[i].P, TreeData[i].psi, &(TreeData[i]) );
			}
			if (TS->hasTPQ == 1)
			{
				kd_clear(TS->TPQ);
				for (i=0; i < L->n; i++)
					kd_insert3(TS->TPQ, TreeData[i].T/300, TreeData[i].P/1E6, TreeData[i].Q, &(TreeData[i]) );
					//kd_insert3(TS->TPQ, TreeData[i].T, TreeData[i].P, TreeData[i].Q, &(TreeData[i]) );
			}
			if (TS->hasPXH == 1)
			{
				kd_clear(TS->PXH);
				for (i=0; i < L->n; i++)
					kd_insert3(TS->PXH, TreeData[i].P/1E6, TreeData[i].psi, TreeData[i].h/1E6, &(TreeData[i]) );
					//kd_insert3(TS->PXH, TreeData[i].P, TreeData[i].psi, TreeData[i].h, &(TreeData[i]) );
			}
			if (TS->hasPXQ == 1)
			{
				kd_clear(TS->PXQ);
				for (i=0; i < L->n; i++)
					kd_insert3(TS->PXQ, TreeData[i].P/1E6, TreeData[i].psi, TreeData[i].Q, &(TreeData[i]) );
					//kd_insert3(TS->PXQ, TreeData[i].P, TreeData[i].psi, TreeData[i].Q, &(TreeData[i]) );
			}
			if (TS->hasXVU == 1)
			{
				kd_clear(TS->XVU);
				for (i=0; i < L->n; i++)
					kd_insert3(TS->XVU, TreeData[i].psi, TreeData[i].v, TreeData[i].u/1E6, &(TreeData[i]) );
					//kd_insert3(TS->PXQ, TreeData[i].psi, TreeData[i].v, TreeData[i].U, &(TreeData[i]) );
			}
			if (TS->hasTXV == 1)
			{
				kd_clear(TS->TXV);
				for (i=0; i < L->n; i++)
					kd_insert3(TS->TXV, TreeData[i].T/300, TreeData[i].psi, TreeData[i].v, &(TreeData[i]) );
			}
		}
	}
	//Add new datapoint
	TreeData		= L->D;
	i				= L->n;
	TreeData[i]		= *P;	//copy over new datapoint
	(L->n)++;				//Increment point counter
	
	//Now put that datapoint in the tree
	if ( TS->hasTPX == 1 )
	{	
		kd_insert3(TS->TPX, TreeData[i].T/300, TreeData[i].P/1E6, TreeData[i].psi, &(TreeData[i]) );
		//kd_insert3(TS->TPX, TreeData[i].T, TreeData[i].P, TreeData[i].psi, &(TreeData[i]) ); 
	}
	if ( TS->hasPXH == 1 )
	{	
		kd_insert3(TS->PXH, TreeData[i].P/1E6, TreeData[i].psi, TreeData[i].h/1E6, &(TreeData[i]) );
		//kd_insert3(TS->PXH, TreeData[i].P, TreeData[i].psi, TreeData[i].h, &(TreeData[i]) );
	}
	if ( TS->hasTPQ == 1 )
	{	
		kd_insert3(TS->TPQ, TreeData[i].T/300, TreeData[i].P/1E6, TreeData[i].Q, &(TreeData[i]) );
		//kd_insert3(TS->TPQ, TreeData[i].T, TreeData[i].P, TreeData[i].Q, &(TreeData[i]) ); 
	}
	if ( TS->hasPXQ == 1 )
	{	
		kd_insert3(TS->PXQ, TreeData[i].P/1E6, TreeData[i].psi, TreeData[i].Q, &(TreeData[i]) ); 
		//kd_insert3(TS->PXQ, TreeData[i].P, TreeData[i].psi, TreeData[i].Q, &(TreeData[i]) ); 
	}
	if ( TS->hasXVU == 1 )
	{	
		kd_insert3(TS->XVU, TreeData[i].psi, TreeData[i].v, TreeData[i].u/1E6, &(TreeData[i]) ); 
		//kd_insert3(TS->XVU, TreeData[i].psi, TreeData[i].v, TreeData[i].u, &(TreeData[i]) ); 
	}
	if ( TS->hasTXV == 1 )
	{	
		kd_insert3(TS->TXV, TreeData[i].T/300, TreeData[i].psi, TreeData[i].v, &(TreeData[i]) ); 
	}
}


//==================================================================================================
//INDIVIDUAL QUERY MODES
//==================================================================================================
//Return mode: 0-no search, 1-had to search a tree (potentially add new point to tree)
int getTPX(Propset *P, void *Tree)
{
	//Back up values that get overwritten by the search
	double	bT, bP, bpsi;
	//double 	Ps_a, Ps_w;
	double 	Tb, Td;
	int 	searched = 0;
	int 	its;		//iteration count
	double	rad;		//distance from closest point in database
	
	//First check if we are definitely in the single pressure regime, and can skip iteration
	//Ps_a	= Psat_a(P->T);
	//Ps_w	= Psat_w(P->T);
	Tb		= Tbub (P->P, P->psi);
	Td		= Tdew (P->P, P->psi);
		
	//if (P->P > Ps_a)	//Definitely in liquid phase
	if (P->T < 0.98*Tb)	//Definitely in liquid phase
	{
		P->x	= P->psi;
		P->y	= 0;
		P->Q	= 0;
		P->MW	= MW_m( P->x );	
		P->v	= V_lm(P->T, P->P, P->x);
		P->rho	= 1/P->v;
		P->cp 	= Cp_lm( P->T, P->P, P->x );
		P->h	= H_lm(P->T, P->P, P->x);
		P->u	= U_lm(P->T, P->P, P->x);
		P->s	= S_lm(P->T, P->P, P->x);
	}
	//else if (P->P < Ps_w)	//Definitely in vapor phase
	else if (P->T > 1.02*Td)	//Definitely in vapor phase
	{
		P->x	= 0;
		P->y	= P->psi;
		P->Q	= 1;
		P->MW	= MW_m( P->y );	
		P->v	= V_vm(P->T, P->P, P->y);
		P->rho	= 1/P->v;
		P->cp 	= Cp_vm( P->T, P->P, P->y );
		P->h	= H_vm(P->T, P->P, P->y);
		P->u	= U_vm(P->T, P->P, P->y);
		P->s	= S_vm(P->T, P->P, P->y);
	}		
	else	//Potentially two-phase mixture, more complicated
	{
		//Back up values that will get overwritten by the search
		bT = P->T; bP = P->P; bpsi = P->psi;
		//Search for initial guess
		PropTreeSearch(P->T/300, P->P/1E6, P->psi, P, Tree);
		//PropTreeSearch(P->T, P->P, P->psi, P, Tree);
		//Check how close we are to a point in the database
		rad	= sqrt( pow( (P->T-bT)/bT, 2 ) + pow( (P->P-bP)/bP, 2 ) + pow( P->psi-bpsi, 2 ) );
		//Replace value that got overwritten by search
		P->T = bT; P->P = bP; P->psi = bpsi;
		//Use linear interpolation between sat points if we are too far from a point in the database
		if( rad > 0.05 )
		{
			P->Q = (P->T - Tb)/(Td - Tb);
			P->Q = median(1E-6, P->Q, 1-1E-6);
			P->x = P->psi *  (Td - P->T)/(Td - Tb) ;
			P->x = median(1E-6, P->x, 1-1E-6);
			P->y = 1 + ( (P->psi - 1)*(P->T - Tb)/(Td - Tb) );
			P->y = median(1E-6, P->y, 1-1E-6);
		}
		
		//Iterate to get solution
		its = Qxy_TPX ( &(P->Q), &(P->x), &(P->y), P->T, P->P, P->psi );
		//Ensure quality is in allowed range, double check if Q should be 0 or 1 when outside of range
		if ( (P->Q < 0) || (P->Q > 1) )
		{
			if ( P->T < ((Tb + Td)/2) )
			{	P->Q	= 0;}
			else
			{	P->Q	= 1;}
		}
		
		//Set searched to 1 if fluid is 2 phase and iteration count is sufficient
		searched = ( (P->Q > 0.0) && (P->Q < 1.0) && (its > 2) ) ? 1 : 0;
		
		//If Quality is 0 or 1, short-circuit the concentrations
		P->x	= (P->Q == 0)  ? P->psi : P->x;
		P->y	= (P->Q == 1)  ? P->psi : P->y;
		
		//Now get the rest of the properties
		P->MW	= MW_m( P->psi );
		P->v    = (P->Q)*V_vm(P->T, P->P, P->y) + (1 - P->Q)*V_lm(P->T, P->P, P->x);
        P->rho	= 1/P->v;
        
		//Get cp, only exists for single phase
		if ( P->Q == 0)
		{	P->cp = Cp_lm( P->T, P->P, P->x ); P->x = P->psi; P->y = 0;}
		else if (P->Q == 1)
		{	P->cp = Cp_vm( P->T, P->P, P->y ); P->x = 0; P->y = P->psi;}
		
		//Get h, u, s
		P->h	= (P->Q)*H_vm(P->T, P->P, P->y) + (1 - P->Q)*H_lm(P->T, P->P, P->x);
		P->u	= (P->Q)*U_vm(P->T, P->P, P->y) + (1 - P->Q)*U_lm(P->T, P->P, P->x);
		P->s	= (P->Q)*S_vm(P->T, P->P, P->y) + (1 - P->Q)*S_lm(P->T, P->P, P->x);
	}
	
	return searched;
}

//Return mode: 0-no search, 1-had to search a tree (potentially add new point to tree)
int getPXH(Propset *P, void *Tree)
{
	//Variables
	double	bP, bpsi, bH;	//Back up values that get overwritten by the search
	double 	Tb, Td;			//Saturation temperatures
	double	hs_l, hs_v;		//Saturation enthalpies
	int 	searched = 0;
	int		its;			//iteration count
	double	rad;			//Radius from point in database
	
	//First check if we are definitely in the single pressure regime, and can skip iteration
	Tb		= Tbub(P->P, P->psi);
	Td		= Tdew(P->P, P->psi);
	hs_l	= H_lm(Tb, P->P, P->psi);	//minimum possible enthalpy for a saturated liquid at this pressure
	hs_v	= H_vm(Td, P->P, P->psi);	//maximum possible enthalpy for a saturated vapor at this pressure
		
	if (P->h < 0.98*hs_l)		//Definitely in liquid phase
	{
		P->x	= P->psi;
		P->y	= 0;
		P->Q	= 0;
		P->MW	= MW_m( P->x );	
		//We have to iterate to get T, start with saturation temp as an initial guess
		P->T	= Tb;
		T_PQHxy ( &(P->T), P->P, P->Q, P->h, P->x, P->y);
		//Now we can get the rest of the parameters
		P->v	= V_lm(P->T, P->P, P->x);
		P->rho	= 1/P->v;
		P->cp 	= Cp_lm( P->T, P->P, P->x );
		P->h	= H_lm(P->T, P->P, P->x);
		P->u	= U_lm(P->T, P->P, P->x);
		P->s	= S_lm(P->T, P->P, P->x);
	}
	else if (P->h > 1.02*hs_v)	//Definitely in vapor phase
	{
		P->x	= 0;
		P->y	= P->psi;
		P->Q	= 1;
		P->MW	= MW_m( P->y );	
		//We have to iterate to get T, start with saturation temp as an initial guess
		P->T	= Td;
		T_PQHxy ( &(P->T), P->P, P->Q, P->h, P->x, P->y);
		//Now we can get the rest of the parameters
		P->v	= V_vm(P->T, P->P, P->y);
		P->rho	= 1/P->v;
		P->cp 	= Cp_vm( P->T, P->P, P->y );
		P->h	= H_vm(P->T, P->P, P->y);
		P->u	= U_vm(P->T, P->P, P->y);
		P->s	= S_vm(P->T, P->P, P->y);
	}
	else	//Potentially two-phase mixture, more complicated
	{
		//Back up values that will get overwritten by the search
		bP = P->P; bpsi = P->psi; bH = P->h;
		//Search for initial guess
		PropTreeSearch(P->P/1E6, P->psi, P->h/1E6, P, Tree);
		//PropTreeSearch(P->P, P->psi, P->h, P, Tree);
		//Check how close we are to a point in the database
		rad	= sqrt( pow( (P->P-bP)/bP, 2 ) + pow( P->psi-bpsi, 2 ) + pow( (P->h-bH)/bH, 2 ) );
		
		//Replace value that got overwritten by search
		P->P = bP; P->psi = bpsi; P->h = bH;
		
		if( rad > 0.05 )
		{
			P->T = Tb + (Td - Tb)*(P->h - hs_l)/(hs_v - hs_l);
			P->Q = (P->h - hs_l)/(hs_v - hs_l);
			P->x = P->psi * (hs_v - P->h)/(hs_v - hs_l);
            P->x = median(1E-6, P->x, 1 - 1E-6);
			P->y = 1 + ( (P->psi - 1)*(P->h - hs_l)/(hs_v - hs_l) );
            P->y = median(1E-6, P->y, 1 - 1E-6);
		}
        
		//Iterate to get solution
		its = TQxy_PXH ( &(P->T), &(P->Q), &(P->x), &(P->y), P->P, P->psi, P->h, 0.1*fabs(Td-Tb) );

		//Check if quality is outside of [0,1], then short circuit to single phase iterative solution
		if (P->Q <= 0)		//Liquid mode
		{
			T_PHx_l (&(P->T), P->P, P->h, P->psi);		//Find temperature for the liquid mixture
			P->Q		= 0;
			P->x		= P->psi;
			P->y		= 0;
		}
		else if (P->Q >= 1)	//Vapor mode
		{
			T_PHy_v (&(P->T), P->P, P->h, P->psi);		//Find temperature for the vapor mixture
			P->Q		= 1;
			P->y		= P->psi;
			P->x		= 0;
		}
		
		//Set searched to 1 if fluid is 2 phase
		searched = ( (P->Q > 0.0) && (P->Q < 1.0) && (its > 3) ) ? 1 : 0;
		
		//Now get the rest of the properties
		P->MW	= MW_m( P->psi );
		P->v    = (P->Q)*V_vm(P->T, P->P, P->y) + (1 - P->Q)*V_lm(P->T, P->P, P->x);
        P->rho	= 1/P->v;
		//Get cp, only exists for single phase
		if ( P->Q == 0)
		{	P->cp = Cp_lm( P->T, P->P, P->x ); P->x = P->psi; P->y = 0;}
		else if (P->Q == 1)
		{	P->cp = Cp_vm( P->T, P->P, P->y ); P->x = 0; P->y = P->psi;}
		
		//Get h, u, s
		P->h	= (P->Q)*H_vm(P->T, P->P, P->y) + (1 - P->Q)*H_lm(P->T, P->P, P->x);
		P->u	= (P->Q)*U_vm(P->T, P->P, P->y) + (1 - P->Q)*U_lm(P->T, P->P, P->x);
		P->s	= (P->Q)*S_vm(P->T, P->P, P->y) + (1 - P->Q)*S_lm(P->T, P->P, P->x);
	}
	
	return searched;
}

//Return mode: 0-no search, 1-had to search a tree (potentially add new point to tree)
int getTPQ(Propset *P, void *Tree)
{
	//Variables
	
	//Back up values that get overwritten by the search
	double	bT, bP, bQ;
	double	xb, yd;
	double	Ps_a, Ps_w;		//Saturation limits
	double	rad;			//Radius from point in database
	int 	searched = 0;
	int 	its = 0;
	
	//Apply limits to P.Q
    if      (P->Q < 0)
    {   P->Q = 0; }
    else if (P->Q > 1)
    {   P->Q = 1; }
	
	//Check if state is within the saturation limits
	Ps_a	= Psat_a(P->T);
	Ps_w	= Psat_w(P->T);
	
	//Get saturation states
	xb		= xbub(P->T, P->P);
	yd		= ydew(P->T, P->P);
	
	//Search a tree for guesses for x, y, psi
	//Back up values that will get overwritten by the search
	bT = P->T; bP = P->P; bQ = P->Q;
	//Search for initial guess
	PropTreeSearch(P->T/300, P->P/1E6, P->Q, P, Tree);
	//Check how close we are to a point in the database
	rad	= sqrt( pow( (P->T-bT)/bT, 2 ) + pow( (P->P-bP)/bP, 2 ) + pow( P->Q-bQ, 2 ) );
			
	//Replace value that got overwritten by search
	P->T = bT; P->P = bP; P->Q = bQ;
		
	//Initial guesses for x, y, psi if no close points in database
	if( rad > 0.05 )
	{
		P->x	= xb;
		P->y	= yd;
		P->psi	= P->x*(1-P->Q) + P->y*P->Q;
	}
    P->x	= xb;
	P->y	= yd;
	
	//Iterative solution for x, y, psi
	P->psi		= (xb + yd)/2;
	its 		= Qxy_TPX ( &(P->Q), &(P->x), &(P->y), P->T, P->P, P->psi );
	P->Q		= bQ;
    P->psi		= (P->Q)*(P->y) + (1-P->Q)*(P->x);
    
	//Set searched to 1 if fluid is 2 phase
	searched = ( (P->Q >= 0.0) && (P->Q <= 1.0) && (its > 3) ) ? 1 : 0;
	
	//Now get the rest of the properties
	P->MW	= MW_m( P->psi );
	P->v    = (P->Q)*V_vm(P->T, P->P, P->y) + (1 - P->Q)*V_lm(P->T, P->P, P->x);
    P->rho	= 1/P->v;
	//Get cp, only exists for single phase
	if ( P->Q == 0)
	{	P->cp = Cp_lm( P->T, P->P, P->x ); P->psi = P->x; P->y = 0;}
	else if (P->Q == 1)
	{	P->cp = Cp_vm( P->T, P->P, P->y ); P->x = 0; P->psi = P->y;}
		
	//Get h, u, s
	P->h	= (P->Q)*H_vm(P->T, P->P, P->y) + (1 - P->Q)*H_lm(P->T, P->P, P->x);
	P->u	= (P->Q)*U_vm(P->T, P->P, P->y) + (1 - P->Q)*U_lm(P->T, P->P, P->x);
	P->s	= (P->Q)*S_vm(P->T, P->P, P->y) + (1 - P->Q)*S_lm(P->T, P->P, P->x);
	
	return searched;
}


//Return mode: 0-no search, 1-had to search a tree (potentially add new point to tree)
int getPXQ(Propset *P, void *Tree)
{
	//Variables
	double	bP, bpsi, bQ;	//Back up values that get overwritten by the search
	double 	Tb, Td;			//Saturation temperatures
	int 	searched = 0;
	int		its;			//iteration count
	double	rad;			//Radius from point in database
	
	//Apply limits to P.Q
	P->Q	= median(1E-6, P->Q, 1-1E-6);
	P->psi	= median(1E-6, P->psi, 1-1E-6);
	
	//Use Tb and Td to estimate initial temperature
	Tb		= Tbub(P->P, P->psi);
	Td		= Tdew(P->P, P->psi);
		
	//Back up values that will get overwritten by the search
	bP = P->P; bpsi = P->psi; bQ = P->Q;
	//Search for initial guess
	PropTreeSearch(P->P/1E6, P->psi, P->Q, P, Tree);
	//PropTreeSearch(P->P, P->psi, P->Q, P, Tree);
	//Check how close we are to a point in the database
	rad	= sqrt( pow( (P->P-bP)/bP, 2 ) + pow( P->psi-bpsi, 2 ) + pow( P->Q - bQ, 2 ) );
		
	//Replace value that got overwritten by search
	P->P = bP; P->psi = bpsi; P->Q = bQ;
		
	if( rad > 0.05 )
	{
		P->T = Tb + (Td - Tb)*P->Q;
		P->x = P->psi * (1 - P->Q);
		P->y = 1 - (1 - P->psi)*P->Q;
	}
	
    //Iterate to get solution
	its = Txy_PXQ ( &(P->T), &(P->x), &(P->y), P->P, P->psi, P->Q, 0.1*fabs(Td-Tb) );

    //mexPrintf("Here, info to Txy_PXQ:  %f  %f  %f  %f  %f  %f  %f\n", P->T, P->x, P->y, P->P, P->psi, P->Q, 0.1*fabs(Td-Tb) );
    
    
	//Something went wrong during iteration:
	if ( its == -1 )
	{
//mexPrintf("Iteration Fail, Psi = %f\n", P->psi);
		//The initial guess from the DB is likely the culprit
		//So reset the initial guess to a manual one
		P->T = Tb + (Td - Tb)*P->Q;
		P->x = P->psi * (1 - P->Q);
		P->y = 1 - (1 - P->psi)*P->Q;
		
		//And rerun the iteration
		its = Txy_PXQ ( &(P->T), &(P->x), &(P->y), P->P, P->psi, P->Q, 0.1*fabs(Td-Tb) );
//mexPrintf("\tRerun its: %d\n", its);
	}

	//Set searched to 1 if fluid is 2 phase
	searched = ( (P->Q > 0.0) && (P->Q < 1.0) && (its > 3) ) ? 1 : 0;
		
	//Now get the rest of the properties
	P->MW	= MW_m( P->psi );
	P->v    = (P->Q)*V_vm(P->T, P->P, P->y) + (1 - P->Q)*V_lm(P->T, P->P, P->x);
    P->rho	= 1/P->v;		
		
	//Get h, u, s
	P->h	= (P->Q)*H_vm(P->T, P->P, P->y) + (1 - P->Q)*H_lm(P->T, P->P, P->x);
	P->u	= (P->Q)*U_vm(P->T, P->P, P->y) + (1 - P->Q)*U_lm(P->T, P->P, P->x);
	P->s	= (P->Q)*S_vm(P->T, P->P, P->y) + (1 - P->Q)*S_lm(P->T, P->P, P->x);
	
	return searched;
}

//Gets TPX from XVU given an old XVU as an initial guess
int getXVU2(Propset *Pin,  Trees *TS, Proplib *PL)
{
    //Variables
	double	bpsi, bV, bU;	//Back up values that get overwritten by the search
	int 	searched = 0;
	//double	rad;			//Radius from point in database
	//FOR Iterative solution
	//Set input values
	Propset		PTPX_n;
    //double      X_i = Pin->psi;
    double      V_i = Pin->v;
    double      U_i = Pin->u;
    //Current values of target properties
    double      V, U;
    double      T, P;
    //Values for estimating derivatives
    double      V_T = 0;
    double      V_P = 0;
    double      U_T = 0;
    double      U_P = 0;
    double      dT = 0.01;
    double      dP = 10;
    //Newton-Rhapson vecs/mats
    double      R[2];                   //Residual
    double      J[2][2], Jinv[2][2];    //Jacobians
    double      dTP[2];                 //change in T, P
    //Iteration parameters
    double      err = 1;
    double      errMax = 1E-5;          //Max error
    int         its = 0;                //Iteration counter
    int         itsMax = 100;           //Maximum iterations
    double      Plim;					//Maximum change in pressure
    double      Tlim;		            //Maximum change in temperature (K)
    double      alphaT, alphaP;         //Limiter factors for T and P
    double      alpha = 1;              //Under-relaxation factor
	
	//Back up values that will get overwritten by the search
	bpsi = Pin->psi; bV = Pin->v; bU = Pin->u;
	//Search for initial guess
	PropTreeSearch(Pin->psi, Pin->v, Pin->u/1E6, &PTPX_n, TS->XVU);
	//Replace value that got overwritten by search
	PTPX_n.psi = bpsi; PTPX_n.v = bV; PTPX_n.u = bU;

	//Set new values of T and P
	T	= PTPX_n.T;
	P	= PTPX_n.P;

	//Set limiters
    Plim = 0.2*PTPX_n.P;
    Tlim = 5;
	
	//Main iteration loop
    while ( (its < itsMax) && (err > errMax) )
    {   
        //Update V and U
        getTPXLite(&PTPX_n, TS->TPX);
        V       = PTPX_n.v;
        U       = PTPX_n.u;
        
        //Get residual
        R[0]    = V_i - V;
        R[1]    = U_i - U;
        
        //Get derivatives
        PTPX_n.T   += dT;
        getTPXLite(&PTPX_n, TS->TPX);
        V_T         = PTPX_n.v;
        U_T         = PTPX_n.u;
        PTPX_n.T   -= dT;
        PTPX_n.P   += dP;
        getTPXLite(&PTPX_n, TS->TPX);
        V_P         = PTPX_n.v;
        U_P         = PTPX_n.u;
        PTPX_n.P   -= dP;
        
        //Get Jacobian
        J[0][0]     = -(V_T - V)/dT;
		J[0][1]     = -(V_P - V)/dP;
        J[1][0]     = -(U_T - U)/dT;
		J[1][1]     = -(U_P - U)/dP;
        
        //Get inverse Jacobian
        Inv2x2(J, Jinv);
        
        //Get change in T and P
        dTP[0]  = -(Jinv[0][0]*R[0] + Jinv[0][1]*R[1]);
        dTP[1]  = -(Jinv[1][0]*R[0] + Jinv[1][1]*R[1]);
        
        //Limiter through alpha
        if ( fabs(dTP[0]) > Tlim )
        {   alphaT = Tlim/fabs(dTP[0]); }
        else
        {    alphaT = 1; }
        if ( fabs(dTP[1]) > Plim )
        {   alphaP = Plim/fabs(dTP[1]); }
        else
        {   alphaP = 1; }
        alpha   = min_val(alphaT, alphaP);
        
        //Update T and P
        T       = T + alpha*dTP[0];
        P       = P + alpha*dTP[1];
        
        //Update property structure
        PTPX_n.T   = T;
        PTPX_n.P   = P;
        
        //Update iteration counter
        its++;
        
        //Update error
        err = sqrt( (dTP[0]*dTP[0])/(T*T) + (dTP[1]*dTP[1])/(P*P) );
    }

	//Set P <- PTPX_n
	Pin->T = PTPX_n.T;
	Pin->P = PTPX_n.P;
	Pin->x = PTPX_n.x;
	Pin->y = PTPX_n.y;
	Pin->Q = PTPX_n.Q;
		
	//Set searched to 1 if fluid is 2 phase
	searched = ( (Pin->Q > 0.0) && (Pin->Q < 1.0) && (its > 3) ) ? 1 : 0;
	
	//Now get the rest of the properties
	Pin->MW	= MW_m( Pin->psi );
	Pin->v	= (Pin->Q)*V_vm(Pin->T, Pin->P, Pin->y) + (1 - Pin->Q)*V_lm(Pin->T, Pin->P, Pin->x);
	Pin->rho	= 1/(Pin->v);
	//Get cp, only exists for single phase
	if ( Pin->Q == 0)
	{	Pin->cp = Cp_lm( Pin->T, Pin->P, Pin->x ); Pin->x = Pin->psi; Pin->y = 0;}
	else if (Pin->Q == 1)
	{	Pin->cp = Cp_vm( Pin->T, Pin->P, Pin->y ); Pin->x = 0; Pin->y = Pin->psi;}
	else
	{	Pin->cp	= 0; }
		
	//Get h, u, s
	Pin->h	= (Pin->Q)*H_vm(Pin->T, Pin->P, Pin->y) + (1 - Pin->Q)*H_lm(Pin->T, Pin->P, Pin->x);
	Pin->u	= (Pin->Q)*U_vm(Pin->T, Pin->P, Pin->y) + (1 - Pin->Q)*U_lm(Pin->T, Pin->P, Pin->x);
	Pin->s	= (Pin->Q)*S_vm(Pin->T, Pin->P, Pin->y) + (1 - Pin->Q)*S_lm(Pin->T, Pin->P, Pin->x);
	
	return searched;
}

//Return mode: 0-no search, 1-had to search a tree (potentially add new point to tree)
int getXVU(Propset *P, void *Tree)
{
	//Variables
	double	bpsi, bV, bU;	//Back up values that get overwritten by the search
	double	dpsi, dV, dU;	//Closest values stored in the database
	double 	Tb, Td;			//Saturation temperatures
	double	Pb, Pd;			//Saturation pressures
	double	Ub, Ud;			//Saturation energies
	double	vb, vd;			//Saturation 
	int 	searched = 0;
	int		its;			//iteration count
	double	rad;			//Radius from point in database
	
	//First check region on the flow map
	Pb		= Pbub (P->psi, P->v);
	Pd		= Pdew (P->psi, P->v);
	Tb		= Tbub(Pb, P->psi);
	Td		= Tdew(Pd, P->psi);
	//Get saturation energies
	Ub		= U_lm(Tb, Pb, P->psi);
	Ud		= U_vm(Td, Pd, P->psi);
	//Get saturation volumes
	vb		= V_lm(Tb, Pb, P->psi);
	vd		= V_vm(Td, Pd, P->psi);
	
	if (P->v < 0.95*vb)		//Definitely in liquid phase
	{
		P->x	= P->psi;
		P->y	= 0;
		P->Q	= 0;
		P->MW	= MW_m( P->x );	
		//We have to iterate to get T and P, start with saturation temp as an initial guess
		P->T	= Tb;
		P->P	= Pb;
		//iterate to find P,T (function doesn't exist yet)
		TP_QVUxy ( &(P->T), &(P->P), P->Q, P->v, P->u, P->x, P->y );
		//Now we can get the rest of the parameters
		P->v	= V_lm(P->T, P->P, P->x);
		P->rho	= 1/P->v;
		P->cp 	= Cp_lm( P->T, P->P, P->x );
		P->h	= H_lm(P->T, P->P, P->x);
		P->u	= U_lm(P->T, P->P, P->x);
		P->s	= S_lm(P->T, P->P, P->x);
	}
	else if (P->v > 1.05*vd )
	{
		P->x	= 0;
		P->y	= P->psi;
		P->Q	= 1;
		P->MW	= MW_m( P->y );	
		//We have to iterate to get T and P, start with saturation temp as an initial guess
		P->T	= Td;
		P->P	= Pd;
		//iterate to find P,T (function doesn't exist yet)
		TP_QVUxy ( &(P->T), &(P->P), P->Q, P->v, P->u, P->x, P->y );
		//Now we can get the rest of the parameters
		P->v	= V_vm(P->T, P->P, P->y);
		P->rho	= 1/P->v;
		P->cp 	= Cp_vm( P->T, P->P, P->y );
		P->h	= H_vm(P->T, P->P, P->y);
		P->u	= U_vm(P->T, P->P, P->y);
		P->s	= S_vm(P->T, P->P, P->y);
	}
	else	//Potentially two-phase mixture, more complicated
	{
		//Back up values that will get overwritten by the search
		bpsi = P->psi; bV = P->v; bU = P->u;
		//Search for initial guess
		PropTreeSearch(P->psi, P->v, P->u/1E6, P, Tree);
		//Store values retrieved from the database
		dpsi = P->psi; dV = P->v; dU = P->u;
		//Check how close we are to a point in the database
		rad	= sqrt( pow( bpsi - P->psi, 2 ) + pow( (P->v-bV)/(P->v+bV) , 2 ) + pow( (P->u-bU)/(P->u+bU), 2 ) );
		
		//Replace value that got overwritten by search
		P->psi = bpsi; P->v = bV; P->u = bU;
		
		if( rad > 0.25 )
		{
			P->T = Tb + (Td - Tb)*(P->v - vb)/(vd - vb);
			P->P = Pb + (Pd - Pb)*(P->v - vb)/(vd - vb);
			P->Q = (P->v - vb)/(vd - vb);
			P->x = P->psi * (vd - P->v)/(vd - vb);
			P->y = 1 + ( (P->psi - 1)*(P->v - vb)/(vd - vb) );
		}
		
		//Iterate to get solution (function doesn't exist yet
		its = TPQxy_XVU ( &(P->T), &(P->P), &(P->Q), &(P->x), &(P->y), P->psi, P->v, P->u, 0.2*fabs(Td-Tb), 0.2*fabs(Pd-Pb) );
		//Failure to converge, return intermediate property point to evaluate before reattempting original XVU
		if ( its < 0 )
		{
			P->psi	= (bpsi + dpsi)/2;
			P->v	= (bV   + dV)  /2;
			P->u	= (bU   + dU)  /2;
			return (-1);
		}
		
		//Check if quality is outside of [0,1], then short circuit to single phase iterative solution
		if (P->Q <= 0)		//Liquid mode
		{
			P->Q		= 0;
			P->x		= P->psi;
			P->y		= 1;
			TP_QVUxy ( &(P->T), &(P->P), P->Q, P->v, P->u, P->x, P->y );
		}
		else if (P->Q >= 1)	//Vapor mode
		{
			P->Q		= 1;
			P->x		= 0;
			P->y		= P->psi;
			TP_QVUxy ( &(P->T), &(P->P), P->Q, P->v, P->u, P->x, P->y );
		}
		
		//Set searched to 1 if fluid is 2 phase
		searched = ( (P->Q > 0.0) && (P->Q < 1.0) && (its > 3) ) ? 1 : 0;
		
		//Now get the rest of the properties
		P->MW	= MW_m( P->psi );
		P->v    = (P->Q)*V_vm(P->T, P->P, P->y) + (1 - P->Q)*V_lm(P->T, P->P, P->x);
        P->rho	= 1/P->v;
		//Get cp, only exists for single phase
		if ( P->Q == 0)
		{	P->cp = Cp_lm( P->T, P->P, P->x ); P->x = P->psi; P->y = 0;}
		else if (P->Q == 1)
		{	P->cp = Cp_vm( P->T, P->P, P->y ); P->x = 0; P->y = P->psi;}
		
		//Get h, u, s
		P->h	= (P->Q)*H_vm(P->T, P->P, P->y) + (1 - P->Q)*H_lm(P->T, P->P, P->x);
		P->u	= (P->Q)*U_vm(P->T, P->P, P->y) + (1 - P->Q)*U_lm(P->T, P->P, P->x);
		P->s	= (P->Q)*S_vm(P->T, P->P, P->y) + (1 - P->Q)*S_lm(P->T, P->P, P->x);
	}
	
	return searched;
}


//Return mode: 0-no search, 1-had to search a tree (potentially add new point to tree)
int getTXV(Propset *P, void *Tree)
{
	double	Pb, Pd;			//Bubble and dew pressures
	double	Vb, Vd;			//Bubble and dew specific volumes
	double	Vl, Vv;			//Temporary liquid and vapor volumes
	double	bT, bpsi, bV;	//Temporary storage for T,X,V
	int 	searched = 0;	//If an extensive search was necessary
	int		its;			//iteration count
	double	rad;			//Radius from point in database
	
	//First get bubble and dew pressures and temps for TX
	Pb		= Pbub_TX (P->T, P->psi);
	Pd		= Pdew_TY (P->T, P->psi);
	Vb		= V_lm (P->T, Pb, P->psi);
	Vd		= V_vm (P->T, Pd, P->psi);
	
	//Check what regime we are in
	if (P->v < 0.9*Vb)		//Definitely in liquid phase
	{
		P->x	= P->psi;
		P->y	= 0;
		P->Q	= 0;
		P->MW	= MW_m( P->x );	
		//We have to iterate to get P, start with saturation pressure as an initial guess
		P->P	= Pb;
		P_TQVxy ( &(P->P), P->T, P->Q, P->v, P->x, P->y);
		//P_TXQxy ( &(P->P), P->T, P->psi, P->Q, P->x, P->y);
		//Now we can get the rest of the parameters
		P->v	= V_lm(P->T, P->P, P->x);
		P->rho	= 1/P->v;
		P->cp 	= Cp_lm( P->T, P->P, P->x );
		P->h	= H_lm(P->T, P->P, P->x);
		P->u	= U_lm(P->T, P->P, P->x);
		P->s	= S_lm(P->T, P->P, P->x);
	}
	else if (P->v > 1.1*Vd)	//Definitely in vapor phase
	{
		P->x	= 0;
		P->y	= P->psi;
		P->Q	= 1;
		P->MW	= MW_m( P->y );	
		//We have to iterate to get P, start with saturation pressure as an initial guess
		P->P	= Pd;
		P_TQVxy ( &(P->P), P->T, P->Q, P->v, P->x, P->y);
		//Now we can get the rest of the parameters
		P->v	= V_vm(P->T, P->P, P->y);
		P->rho	= 1/P->v;
		P->cp 	= Cp_vm( P->T, P->P, P->y );
		P->h	= H_vm(P->T, P->P, P->y);
		P->u	= U_vm(P->T, P->P, P->y);
		P->s	= S_vm(P->T, P->P, P->y);
	}
	else	//Potentially two-phase mixture, more complicated
	{
		//Back up values that will get overwritten by the search
		bT = P->T; bpsi = P->psi; bV = P->v;
		//Search for initial guess
		PropTreeSearch(P->T/300, P->x, P->v, P, Tree);
		//Check how close we are to a point in the database
		rad	= sqrt( pow( (P->T-bT)/bT, 2 ) + pow( P->psi-bpsi, 2 ) + pow( (P->v-bV)/bV, 2 ) );
		
		//Replace value that got overwritten by search
		P->T = bT; P->psi = bpsi; P->v = bV;
		
		if( rad > 0.05 )
		{
			P->P = Pb + (Pd - Pb)*(P->v - Vb)/(Vd - Vb);
			P->Q = (P->v - Vb)/(Vd - Vb);
			if	 	(P->Q < 0)
			{	P->Q	= 0;}
			else if	(P->Q > 1)
			{	P->Q	= 1;}
			P->x = P->psi * (Vd - P->v)/(Vd - Vb);
			P->x = median(1E-6, P->x, 1 - 1E-6);
			P->y = 1 + ( (P->psi - 1)*(P->v - Vb)/(Vd - Vb) );
			P->y = median(1E-6, P->y, 1 - 1E-6);
		}
		
		//Iterate to get solution
		its = PQxy_TXV ( &(P->P), &(P->Q), &(P->x), &(P->y), P->T, P->psi, P->v );
		
		//Sometimes when Q is outside of (0, 1) it helps to check if Q should be 0 or 1
		if ( ( P->Q <= 0 ) || (P->Q >= 1) )
		{
			Vl		= V_lm(P->T, P->P, P->psi);
			Vv		= V_vm(P->T, P->P, P->psi);
			
			if ( fabs(Vl-P->v) < fabs(Vv-P->v) )
			{	P->Q = 0;	}
			else
			{	P->Q = 1;	}
		}
		
		//Check if quality is outside of [0,1], then short circuit to single phase iterative solution
		if (P->Q <= 0)		//Liquid mode
		{
			P->Q		= 0;
			P->x		= P->psi;
			P->y		= 0;
			P_TQVxy ( &(P->P), P->T, P->Q, P->v, P->x, P->y);
		}
		else if (P->Q >= 1)	//Vapor mode
		{
			P->Q		= 1;
			P->y		= P->psi;
			P->x		= 0;
			P_TQVxy ( &(P->P), P->T, P->Q, P->v, P->x, P->y);
		}
		
		//Set searched to 1 if fluid is 2 phase
		searched = ( (P->Q > 0.0) && (P->Q < 1.0) && (its > 3) ) ? 1 : 0;
		
		//Now get the rest of the properties
		P->MW	= MW_m( P->psi );
		P->v    = (P->Q)*V_vm(P->T, P->P, P->y) + (1 - P->Q)*V_lm(P->T, P->P, P->x);
        P->rho	= 1/P->v;
		//Get cp, only exists for single phase
		if ( P->Q == 0)
		{	P->cp = Cp_lm( P->T, P->P, P->x ); P->x = P->psi; P->y = 0;}
		else if (P->Q == 1)
		{	P->cp = Cp_vm( P->T, P->P, P->y ); P->x = 0; P->y = P->psi;}
		
		//Get h, u, s
		P->h	= (P->Q)*H_vm(P->T, P->P, P->y) + (1 - P->Q)*H_lm(P->T, P->P, P->x);
		P->u	= (P->Q)*U_vm(P->T, P->P, P->y) + (1 - P->Q)*U_lm(P->T, P->P, P->x);
		P->s	= (P->Q)*S_vm(P->T, P->P, P->y) + (1 - P->Q)*S_lm(P->T, P->P, P->x);
	}
	
	return searched;
}

//Return mode: 0-no search, 1-had to search a tree (potentially add new point to tree)
int getTPXLite(Propset *P, void *Tree)
{
	//double 	Ps_a, Ps_w;
	double 	Tb, Td;
	int 	its;		//iteration count
	
	//First check if we are definitely in the single pressure regime, and can skip iteration
	Tb		= Tbub (P->P, P->psi);
	Td		= Tdew (P->P, P->psi);
		
	//if (P->P > Ps_a)	//Definitely in liquid phase
	if (P->T < 0.98*Tb)	//Definitely in liquid phase
	{
		P->x	= P->psi;
		P->y	= 0;
		P->Q	= 0;
		P->v	= V_lm(P->T, P->P, P->x);
		P->u	= U_lm(P->T, P->P, P->x);
	}
	//else if (P->P < Ps_w)	//Definitely in vapor phase
	else if (P->T > 1.02*Td)	//Definitely in vapor phase
	{
		P->x	= 0;
		P->y	= P->psi;
		P->Q	= 1;
		P->v	= V_vm(P->T, P->P, P->y);
		P->u	= U_vm(P->T, P->P, P->y);
	}		
	else	//Potentially two-phase mixture, more complicated
	{	
		//Iterate to get solution
		its = Qxy_TPX ( &(P->Q), &(P->x), &(P->y), P->T, P->P, P->psi );
		//Ensure quality is in allowed range, double check if Q should be 0 or 1 when outside of range
		if ( (P->Q < 0) || (P->Q > 1) )
		{
			if ( P->T < ((Tb + Td)/2) )
			{	P->Q	= 0;}
			else
			{	P->Q	= 1;}
		}
		
		//If Quality is 0 or 1, short-circuit the concentrations
		P->x	= (P->Q == 0)  ? P->psi : P->x;
		P->y	= (P->Q == 1)  ? P->psi : P->y;
		
		//Now get the rest of the properties
		P->v	= (P->Q)*V_vm(P->T, P->P, P->y) + (1 - P->Q)*V_lm(P->T, P->P, P->x);
		//Get cp, only exists for single phase
		if ( P->Q == 0)
		{	P->x = P->psi; P->y = 0;}
		else if (P->Q == 1)
		{	P->x = 0; P->y = P->psi;}
		
		//Get u
		P->u	= (P->Q)*U_vm(P->T, P->P, P->y) + (1 - P->Q)*U_lm(P->T, P->P, P->x);
	}
	
	return 0;
}

//==================================================================================================
//Special modes
//==================================================================================================

//returns a 9 long vector with the values: dT/dX, dT/dv, dT/du, dP/dX, dP/dv, dP/du, dX/dX, dX/dv, dX/du
void getDTPX_DXVU(Propset *P, Trees *TS, Proplib *PL, double *Jvec)
{
    //Variables
    double  J[3][3], Jinv[3][3];                    //Jacobian and inverse Jacobian matrix
    //Step sizes in T, P, x
    double  dT = 0.01;                              //Change in K
    double  dP = 1E-4*P->P;                         //Change in Pa
    double  XP = min_val(P->psi+1E-4, 1-1E-6);      //Values of psi
    double  XM = max_val(P->psi-1E-4, 1E-6);        //Values of psi
    //Stored quantities
    double  X_TP, X_TM, V_TP, V_TM, U_TP, U_TM;     //Values from changing T
    double  X_PP, X_PM, V_PP, V_PM, U_PP, U_PM;     //Values from changing P
    double  X_XP, X_XM, V_XP, V_XM, U_XP, U_XM;     //Values from changing X
    
    //Evaluate values at T+dT/2, T-dT/2
    P->T += dT/2;
    SetPropsSearch("TPX", P, TS, PL);
    X_TP = P->psi; V_TP = P->v; U_TP = P->u;
    P->T -= dT;
    SetPropsSearch("TPX", P, TS, PL);
    X_TM = P->psi; V_TM = P->v; U_TM = P->u;
    P->T += dT/2;
    
    //Evaluate values at P+dP/2, P-dP/2
    P->P += dP/2;
    SetPropsSearch("TPX", P, TS, PL);
    X_PP = P->psi; V_PP = P->v; U_PP = P->u;
    P->P -= dP;
    SetPropsSearch("TPX", P, TS, PL);
    X_PM = P->psi; V_PM = P->v; U_PM = P->u;
    P->P += dP/2;
    
    //Evaluate values at XP, XM
    P->psi = XP;
    SetPropsSearch("TPX", P, TS, PL);
    X_XP = P->psi; V_XP = P->v; U_XP = P->u;
    P->psi = XM;
    SetPropsSearch("TPX", P, TS, PL);
    X_XM = P->psi; V_XM = P->v; U_XM = P->u;
    
    //Fill in approximate Jacobian
    J[0][0] = (X_TP - X_TM)/dT;
    J[0][1] = (X_PP - X_PM)/dP;
    J[0][2] = (X_XP - X_XM)/(XP - XM);
    J[1][0] = (V_TP - V_TM)/dT;
    J[1][1] = (V_PP - V_PM)/dP;
    J[1][2] = (V_XP - V_XM)/(XP - XM);
    J[2][0] = (U_TP - U_TM)/dT;
    J[2][1] = (U_PP - U_PM)/dP;
    J[2][2] = (U_XP - U_XM)/(XP - XM);
    //Get inverse of Jacobian
    Inv3x3(J, Jinv);
    
    //Return values as a 1x9 vector
    Jvec[0] = Jinv[0][0]; Jvec[1] = Jinv[0][1]; Jvec[2] = Jinv[0][2];
    Jvec[3] = Jinv[1][0]; Jvec[4] = Jinv[1][1]; Jvec[5] = Jinv[1][2];
    Jvec[6] = Jinv[2][0]; Jvec[7] = Jinv[2][1]; Jvec[8] = Jinv[2][2];
}

//Gets TPX from XVU given an old XVU as an initial guess
void getTPX_XVU(Propset *PTPX_n, Propset *PXVU_n, Propset *PTPX_p,  Trees *TS, Proplib *PL)
{
    //Set input values
    double      X_i = PXVU_n->psi;
    double      V_i = PXVU_n->v;
    double      U_i = PXVU_n->u;
    //Current values of target properties
    double      V, U;
    double      T, P;
    //Values for estimating derivatives
    double      V_T = 0;
    double      V_P = 0;
    double      U_T = 0;
    double      U_P = 0;
    double      dT = 0.01;
    double      dP = 10;
    //Newton-Rhapson vecs/mats
    double      R[2];                   //Residual
    double      J[2][2], Jinv[2][2];    //Jacobians
    double      dTP[2];                 //change in T, P
    //Iteration parameters
    double      err = 1;
    double      errMax = 1E-5;          //Max error
    int         its = 0;                //Iteration counter
    int         itsMax = 100;           //Maximum iterations
    double      Plim = 0.05*PTPX_p->P;  //Maximum change in pressure
    double      Tlim = 0.4;             //Maximum change in temperature (K)
    double      alphaT, alphaP;         //Limiter factors for T and P
    double      alpha = 1;              //Under-relaxation factor
    
    //Initialize PTPX_n
    *PTPX_n     = *PTPX_p;
    PTPX_n->psi = PXVU_n->psi;
    
    //Set initial values of T and P
    T       = PTPX_p->T;
    P       = PTPX_p->P;
    
    //Main iteration loop
    while ( (its < itsMax) && (err > errMax) )
    {   
        //Update V and U
        SetPropsSearch("TPX", PTPX_n, TS, PL);
        V       = PTPX_n->v;
        U       = PTPX_n->u;
        
        //Get residual
        R[0]    = V_i - V;
        R[1]    = U_i - U;
        
        //Get derivatives
        PTPX_n->T   += dT;
        SetPropsSearch("TPX", PTPX_n, TS, PL);
        V_T         = PTPX_n->v;
        U_T         = PTPX_n->u;
        PTPX_n->T   -= dT;
        PTPX_n->P   += dP;
        SetPropsSearch("TPX", PTPX_n, TS, PL);
        V_P         = PTPX_n->v;
        U_P         = PTPX_n->u;
        PTPX_n->P   -= dP;
        
        //Get Jacobian
        J[0][0]     = -(V_T - V)/dT;
		J[0][1]     = -(V_P - V)/dP;
        J[1][0]     = -(U_T - U)/dT;
		J[1][1]     = -(U_P - U)/dP;
        
        //Get inverse Jacobian
        Inv2x2(J, Jinv);
        
        //Get change in T and P
        dTP[0]  = -(Jinv[0][0]*R[0] + Jinv[0][1]*R[1]);
        dTP[1]  = -(Jinv[1][0]*R[0] + Jinv[1][1]*R[1]);
        
        //Limiter through alpha
        if ( fabs(dTP[0]) > Tlim )
        {   alphaT = Tlim/fabs(dTP[0]); }
        else
        {    alphaT = 1; }
        if ( fabs(dTP[1]) > Plim )
        {   alphaP = Plim/fabs(dTP[1]); }
        else
        {   alphaP = 1; }
        alpha   = min_val(alphaT, alphaP);
        
        //Update T and P
        T       = T + alpha*dTP[0];
        P       = P + alpha*dTP[1];
        
        //Update property structure
        PTPX_n->T   = T;
        PTPX_n->P   = P;
        
        //Update iteration counter
        its++;
        
        //Update error
        err = sqrt( (dTP[0]*dTP[0])/(T*T) + (dTP[1]*dTP[1])/(P*P) );
    }
}


//==================================================================================================
//WATER PROPERTIES BELOW
//==================================================================================================

//Evaluates the Gibbs energy of liquid water
double G_lw (double T, double P)
{
	//Define variable
	double		G_r_lw;		//Reduced Gibbs energy of liquid water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculation of reduced Gibbs energy
	G_r_lw		= H_R_L_O_W - T_r*S_R_L_O_W + B_1_W*(T_r - T_R_O_W) + (B_2_W/2)*( pow(T_r,2) - pow(T_R_O_W,2) ) + (B_3_W/3)*( pow(T_r,3) - pow(T_R_O_W,3) ) - B_1_W*T_r*log(T_r/T_R_O_W) - B_2_W*T_r*(T_r-T_R_O_W) - (B_3_W/2)*T_r*(  pow(T_r,2) - pow(T_R_O_W,2) ) + (A_1_W + A_3_W*T_r + A_4_W*pow(T_r,2) )*(P_r - P_R_O_W) + (A_2_W/2)*( pow(P_r,2) - pow(P_R_O_W,2) );
	
	//Return Gibbs energy
	return ((G_r_lw * R_U * T_B / MW_W) - G_REF_L_W);
}

//Evaluates the enthalpy of liquid water
double H_lw (double T, double P)
{
	//Define variable
	double		H_r_lw;		//Reduced enthalpy of liquid water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced enthalpy of liquid water
	H_r_lw	= H_R_L_O_W + B_1_W*(T_r - T_R_O_W) + (B_2_W/2)*( pow(T_r,2) - pow(T_R_O_W,2) ) + (B_3_W/3)*( pow(T_r,3) - pow(T_R_O_W,3) ) + A_1_W*(P_r - P_R_O_W) + (A_2_W/2)*( pow(P_r,2) - pow(P_R_O_W,2) ) - A_4_W*pow(T_r,2)*(P_r - P_R_O_W);

	//Return enthalpy
	return ((H_r_lw * R_U * T_B / MW_W) - H_REF_L_W);
}

//Evaluates the specific volume of liquid water
double V_lw (double T, double P)
{
	//Define variable
	double		V_r_lw;		//Reduced volume of liquid water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced volume of liquid water
	V_r_lw	= A_1_W + A_3_W*T_r + A_4_W*pow(T_r,2) + A_2_W*P_r;

	//Return specific volume
	return ( V_r_lw * R_U * T_B / (P_B * MW_W) );
}

//Evaluates the internal energy of liquid water
double U_lw (double T, double P)
{
	//Define variable
	double		U_r_lw;		//Reduced energy of liquid water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced internal energy of liquid water
	U_r_lw	= H_R_L_O_W + B_1_W*(T_r - T_R_O_W) + (B_2_W/2)*( pow(T_r,2) - pow(T_R_O_W,2) ) + (B_3_W/3)*( pow(T_r,3) - pow(T_R_O_W,3) ) - A_1_W*P_R_O_W - (A_2_W/2)*( pow(P_r,2) + pow(P_R_O_W,2) ) - A_3_W*T_r*P_r - A_4_W*pow(T_r,2)*(2*P_r - P_R_O_W);

	//Return internal energy
	return ((U_r_lw * R_U * T_B / MW_W) - U_REF_L_W);
}

//Evaluates the density of liquid water
double Rho_lw (double T, double P)
{
	return (1.0/V_lw(T,P));
}

//Evaluates the entropy of liquid water
double S_lw (double T, double P)
{
	//Define variable
	double		S_r_lw;		//Reduced entropy of liquid water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced entropy of liquid water
	S_r_lw	= S_R_L_O_W + B_1_W*log(T_r / T_R_O_W) + B_2_W*( T_r - T_R_O_W ) + (B_3_W/2)*( pow(T_r,2) - pow(T_R_O_W,2) ) - A_3_W*(P_r - P_R_O_W) - 2*A_4_W*T_r*(P_r - P_R_O_W);

	//Return entropy
	return ((S_r_lw * R_U / MW_W) - S_REF_L_W);
}

//Evaluates the specific heat at constant pressure of liquid water
double Cp_lw (double T, double P)
{
	//Define variable
	double		Cp_r_lw;	//Reduced Cp of liquid water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced Cp of liquid water
	Cp_r_lw	= B_1_W + B_2_W*T_r + B_3_W*pow(T_r,2) - 2*A_4_W*T_r*(P_r - P_R_O_W);

	//Return Cp
	return ( Cp_r_lw * R_U / MW_W );
}

//Chemical potential of liquid water
double ChemPo_lw (double T, double P, double x)
{
	double g_l;					//Gibbs energy of liquid mixture
	double T_r, P_r;			//Reduced quantities
	double x_mo;				//Molar ammonia fraction
	double MW_M;				//Molecular weight of the mixture
	double f_1, f_2, f_3;		//Components of excess enthalpy
	double h_la, h_lw;			//Enthalpies of each component
	double s_la, s_lw;			//Entropy of the components
	double dhedx, dsedx, dsmdx;	//Derivatives of the excess enthalpy and entropy, and entropy of mixing
	double mu;					//Chemical potential
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculate molar quantities
	x_mo	= (x/MW_A)/ ( (x/MW_A) + ((1-x)/MW_W) );
	MW_M	= x*MW_A + (1-x)*MW_W;
	
	//Gibbs energy of the mixture
	g_l		= G_lm(T, P, x)*MW_M;
	
	//Enthalpies of the components
	h_lw	= H_lw(T, P)*MW_W;
	h_la	= H_la(T, P)*MW_A;
	
	//Calculate excess enthalpy derivative
	f_1		= E_1 + E_2*P_r + 2*E_5/T_r + 3*E_6/pow(T_r,2);
	f_2		= E_7 + E_8*P_r + 2*E_11/T_r + 3*E_12/pow(T_r,2);
	f_3		= E_13 + E_14*P_r + 2*E_15/T_r + 3*E_16/pow(T_r, 2);
	dhedx	= R_U*T_B*( f_1*(-2*x_mo + 1) + f_2*(-6*pow(x_mo,2) + 6*x_mo - 1) + f_3*(-16*pow(x_mo,3) + 24*pow(x_mo,2) - 10*x_mo + 1) );
	
	//Calculate entropies
	s_lw	= S_lw(T, P)*MW_W;
	s_la	= S_la(T, P)*MW_A;
	
	//Calculate excess entropy derivative
	f_1		= -(E_3 + E_4*P_r) +E_5/pow(T_r,2) + 2*E_6/pow(T_r,3);
	f_2		= -(E_9 + E_10*P_r) + E_11/pow(T_r,2) + 2*E_12/pow(T_r,3);
	f_3		= E_15/pow(T_r,2) + 2*E_16/pow(T_r,3);
	dsedx	= R_U*( f_1*(-2*x_mo + 1) + f_2*(-6*pow(x_mo,2) + 6*x_mo - 1) + f_3*(-16*pow(x_mo,3) + 24*pow(x_mo,2) - 10*x_mo + 1) );
	
	//Entropy of mixing
	if (x_mo == 0)
	{	dsmdx = 1E6;}	//Some big number
	else if (x_mo == 1)
	{	dsmdx = -1E6;}	//Some small number
	else
	{	dsmdx	= -R_U*log(x_mo / (1-x_mo));}
	
	
	//return chemical potential
	mu		= g_l - x_mo*(h_la - h_lw + dhedx - T*( s_la - s_lw + dsedx + dsmdx ) );
	return (mu/MW_W);
}

//Evaluates the Gibbs energy of vapor water
double G_vw (double T, double P)
{
	//Define variable
	double		G_r_vw;		//Reduced Gibbs energy of vapor water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculation of reduced Gibbs energy
	G_r_vw		= H_R_V_O_W - T_r*S_R_V_O_W + D_1_W*(T_r - T_R_O_W) + (D_2_W/2)*( pow(T_r,2) - pow(T_R_O_W,2) ) + (D_3_W/3)*( pow(T_r,3) - pow(T_R_O_W,3) ) - D_1_W*T_r*log(T_r/T_R_O_W) - D_2_W*T_r*(T_r - T_R_O_W) - (D_3_W/2)*T_r*( pow(T_r,2) - pow(T_R_O_W,2) ) + T_r*log(P_r/P_R_O_W) + C_1_W*(P_r-P_R_O_W) + C_2_W*( P_r/pow(T_r,3) - 4*P_R_O_W/pow(T_R_O_W,3) + 3*P_R_O_W*T_r/pow(T_R_O_W,4) ) + C_3_W*(P_r/pow(T_r,11) - 12*P_R_O_W/pow(T_R_O_W,11) + 11*P_R_O_W*T_r/pow(T_R_O_W,12) ) + (C_4_W/3)*( pow(P_r,3)/pow(T_r,11) - 12*pow(P_R_O_W,3)/pow(T_R_O_W,11) + 11*pow(P_R_O_W,3)*T_r/pow(T_R_O_W,12) );	
	
	//Return Gibbs energy
	return ((G_r_vw * R_U * T_B / MW_W) - G_REF_V_W);
}


//Evaluates the Enthalpy of vapor water
double H_vw (double T, double P)
{
	//Define variable
	double		H_r_vw;		//Reduced Enthalpy of vapor water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculation of reduced Enthalpy
	H_r_vw		= H_R_V_O_W + (4*C_2_W)*(P_r/pow(T_r,3) - P_R_O_W/pow(T_R_O_W,3) ) + (D_2_W/2)*( pow(T_r,2) - pow(T_R_O_W,2) ) + (D_3_W/3)*( pow(T_r,3) - pow(T_R_O_W,3) ) + 12*C_3_W*(P_r/pow(T_r,11) - P_R_O_W/pow(T_R_O_W,11) ) + 4*C_4_W*( pow(P_r,3)/pow(T_r,11) - pow(P_R_O_W,3)/pow(T_R_O_W,11) )  + D_1_W*(T_r - T_R_O_W) + C_1_W*(P_r - P_R_O_W);
	
	
	//Return Enthalpy
	return ((H_r_vw * R_U * T_B / MW_W) - H_REF_V_W);
}

//Evaluates the specific volume of vapor water
double V_vw (double T, double P)
{
	//Define variable
	double		V_r_vw;		//Reduced volume of vapor water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced volume of vapor water
	V_r_vw	= T_r/P_r + C_1_W + C_2_W/pow(T_r,3) + C_3_W/pow(T_r,11) + C_4_W*pow(P_r,2)/pow(T_r,11);

	//Return specific volume
	return ( V_r_vw * R_U * T_B / (P_B * MW_W) );
}

//Evaluates the density of vapor water
double Rho_vw (double T, double P)
{
	return (1.0/V_vw(T,P));
}

//Evaluates the entropy of vapor water
double S_vw (double T, double P)
{
	//Define variable
	double		S_r_vw;		//Reduced entropy of vapor water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced entropy of vapor water
	S_r_vw	= S_R_V_O_W + D_1_W*log(T_r/T_R_O_W) + D_2_W*(T_r - T_R_O_W) + (D_3_W/2)*( pow(T_r,2) - pow(T_R_O_W,2) ) - log( P_r / P_R_O_W ) + (3*C_2_W)*( P_r/pow(T_r,4) - P_R_O_W/pow(T_R_O_W,4) ) + (11*C_3_W)*(P_r/pow(T_r,12) - P_R_O_W/pow(T_R_O_W,12) ) + (11/3)*C_4_W*( pow(P_r,3)/pow(T_r,12) - pow(P_R_O_W,3)/pow(T_R_O_W,12) );	
	
	
	//Return vapor
	return ((S_r_vw * R_U / MW_W) - S_REF_V_W);
}

//Evaluates the internal energy of vapor water
double U_vw (double T, double P)
{
	//Define variable
	double		U_r_vw;		//Reduced energy of vapor water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced internal energy of vapor water
	U_r_vw	= H_R_V_O_W + (C_2_W )*(3*P_r/pow(T_r,3) - 4*P_R_O_W/pow(T_R_O_W,3) ) + (D_2_W/2)*( pow(T_r,2) - pow(T_R_O_W,2) ) + (D_3_W/3)*( pow(T_r,3) - pow(T_R_O_W,3) ) + C_3_W*(11*P_r/pow(T_r,11) - 12*P_R_O_W/pow(T_R_O_W,11) ) + C_4_W*( 3*pow(P_r,3)/pow(T_r,11) - 4*pow(P_R_O_W,3)/pow(T_R_O_W,11) ) + D_1_W*(T_r - T_R_O_W) - C_1_W*P_R_O_W - T_r;
	
	
	//Return internal energy
	return ((U_r_vw * R_U * T_B / MW_W) - U_REF_V_W);
}

//Evaluates the specific heat at constant pressure of vapor water
double Cp_vw (double T, double P)
{
	//Define variable
	double		Cp_r_vw;	//Reduced Cp of liquid water
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced Cp of vapor water
	Cp_r_vw	= D_1_W + D_2_W*T_r + D_3_W*pow(T_r,2) - 12*C_2_W*P_r/pow(T_r,4) - (11*P_r/pow(T_r,12) )*( 12*C_3_W + 4*C_4_W*pow(P_r,2) );

	//Return Cp
	return ( Cp_r_vw * R_U / MW_W );
}

//Chemical potential of vapor water
double ChemPo_vw (double T, double P, double y)
{
	double g_v;					//Gibbs energy of vapor mixture
	double y_mo;				//Molar ammonia fraction
	double MW_M;				//Molecular weight of the mixture
	double h_va, h_vw;			//Enthalpies of each component
	double s_va, s_vw;			//Entropy of the components
	double dsmdx;	//Derivatives of the excess enthalpy and entropy, and entropy of mixing
	double mu;					//Chemical potential
	
	//Calculate molar quantities
	y_mo	= (y/MW_A)/ ( (y/MW_A) + ((1-y)/MW_W) );
	MW_M	= y*MW_A + (1-y)*MW_W;
	
	//Gibbs energy of the mixture
	g_v		= G_vm(T, P, y)*MW_M;
	
	//Enthalpies of the components
	h_vw	= H_vw(T, P)*MW_W;
	h_va	= H_va(T, P)*MW_A;
	
	//Calculate entropies
	s_vw	= S_vw(T, P)*MW_W;
	s_va	= S_va(T, P)*MW_A;
	
	//Entropy of mixing
	if (y_mo == 0)
	{	dsmdx = 1E6;}	//Some big number
	else if (y_mo == 1)
	{	dsmdx = -1E6;}	//Some small number
	else
	{	dsmdx	= -R_U*log(y_mo / (1-y_mo));}
	
	
	//return chemical potential
	mu		= g_v - y_mo*(h_va - h_vw - T*( s_va - s_vw + dsmdx ) );
	//return mu;
	return (mu/MW_W);
}

//Saturation pressure of water (T in K, Psat in Pa)
double Psat_w (double T)
{
	//Define some normalize variables
	double Tr	= T/TC_W;						//Reduced temperature
	double Ttr	= TT_W/TC_W;					//Reduced triple point temperature of water
	double t	= (Tr - Ttr)/(1 - Ttr);			//Normalized temperature
	double Ptr	= PT_W/PC_W;					//Reduced triple pressure
	double tb	= ((TB_W/TC_W) - Ttr)/(1-Ttr);	//Normalized boiling temperature of water
	double Psatr;									//Reduced saturation pressure
	
	//Calculate reduced saturation pressure
	Psatr = exp( (Ttr*log(Ptr))*(1-t) /( (t*(1-Ttr) + Ttr)*(1+F_1_W*t)*(1 + F_2_W*t*(t-tb)) ) );

	//Return saturation pressure
	return (Psatr * PC_W);
}

//Saturation temperature of water (P in Pa, Tsat in K)
double Tsat_w (double P)
{
	//Define some normalize variables
	double Tr;									//Reduced temperature
	double Ttr	= TT_W/TC_W;					//Reduced triple point temperature
	double t;									//Normalized temperature
	double Ptr	= PT_W/PC_W;					//Reduced triple pressure
	double tb	= ((TB_W/TC_W) - Ttr)/(1-Ttr);	//Normalized boiling temperature
	double Psatr = P/PC_W;						//Reduced saturation pressure
	
	//iteration stuff
	double	eps		= 1E-6;						//convergence criterion
	double	err		= 1;						//Current error
	double	R;									//Current residual
	double	J;									//Current derivative of the residual
	double	dt;									//Change in normalized temperature
	
	//Initial guess for t
	t	= Psatr;
	
	//Iteration
	while ( err > eps )
	{
		//Calculate residual
		R	= log(Psatr) - 	(Ttr*log(Ptr))*(1-t) /( (t*(1-Ttr) + Ttr)*(1+F_1_W*t)*(1 + F_2_W*t*(t-tb)) );
		
		//Calculate derivative of residual
		J	= Ttr*log(Ptr)*( 1 + 2*F_1_W*t + 2*F_2_W*t - F_2_W*tb - 2*Ttr*F_1_W*t + F_1_W*F_2_W*tb*pow(t,2) + 2*F_1_W*F_2_W*Ttr*pow(t,3) - 4*F_1_W*F_2_W*Ttr*pow(t,2) - 2*F_1_W*F_2_W*tb*t + 2*F_1_W*F_2_W*Ttr*t - F_1_W*F_2_W*Ttr*tb - F_1_W*pow(t,2) - F_2_W*pow(t,2) + (F_1_W+F_2_W)*Ttr + F_1_W*Ttr*pow(t,2) - 2*F_1_W*F_2_W*pow(t,3) + F_2_W*Ttr*pow(t,2) - 2*F_2_W*Ttr*t + 3*F_1_W*F_2_W*pow(t,2) - F_1_W*F_2_W*Ttr*tb*pow(t,2) + 2*F_1_W*F_2_W*Ttr*tb*t )/ pow( (-t + t*Ttr - Ttr)*(1 + F_1_W*t)*(1 + F_2_W*t - F_2_W*tb), 2);
	
		//calculate change in t
		dt	= -R/J;
		t	+= dt;
	
		//update error
		err = dt;
	}

	//Un-normalize T
	Tr		= t*(1-Ttr) + Ttr;

	//Return saturation temperature
	return (Tr * TC_W);
}

//==================================================================================================
//AMMONIA PROPERTIES BELOW
//==================================================================================================

//Evaluates the Gibbs energy of liquid ammonia
double G_la (double T, double P)
{
	//Define variable
	double		G_r_la;		//Reduced Gibbs energy of liquid ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculation of reduced Gibbs energy
	G_r_la		= H_R_L_O_A - T_r*S_R_L_O_A + B_1_A*(T_r - T_R_O_A) + (B_2_A/2)*( pow(T_r,2) - pow(T_R_O_A,2) ) + (B_3_A/3)*( pow(T_r,3) - pow(T_R_O_A,3) ) - B_1_A*T_r*log(T_r/T_R_O_A) - B_2_A*T_r*(T_r-T_R_O_A) - (B_3_A/2)*T_r*(  pow(T_r,2) - pow(T_R_O_A,2) ) + (A_1_A + A_3_A*T_r + A_4_A*pow(T_r,2) )*(P_r - P_R_O_A) + (A_2_A/2)*( pow(P_r,2) - pow(P_R_O_A,2) );
	
	//Return Gibbs energy
	return ((G_r_la * R_U * T_B / MW_A) - G_REF_L_A);
}

//Evaluates the enthalpy of liquid ammonia
double H_la (double T, double P)
{
	//Define variable
	double		H_r_la;		//Reduced enthalpy of liquid ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced enthalpy of liquid ammonia
	H_r_la	= H_R_L_O_A + B_1_A*(T_r - T_R_O_A) + (B_2_A/2)*( pow(T_r,2) - pow(T_R_O_A,2) ) + (B_3_A/3)*( pow(T_r,3) - pow(T_R_O_A,3) ) + A_1_A*(P_r - P_R_O_A) + (A_2_A/2)*( pow(P_r,2) - pow(P_R_O_A,2) ) - A_4_A*pow(T_r,2)*(P_r - P_R_O_A);

	//Return enthalpy
	return ((H_r_la * R_U * T_B / MW_A) - H_REF_L_A);
}

//Evaluates the specific volume of liquid ammonia
double V_la (double T, double P)
{
	//Define variable
	double		V_r_la;		//Reduced volume of liquid ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced volume of liquid ammonia
	V_r_la	= A_1_A + A_3_A*T_r + A_4_A*pow(T_r,2) + A_2_A*P_r;

	//Return specific volume
	return ( V_r_la * R_U * T_B / (P_B * MW_A) );
}

//Evaluates the internal energy of liquid ammonia
double U_la (double T, double P)
{
	//Define variable
	double		U_r_la;		//Reduced energy of liquid ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced internal energy of liquid ammonia
	U_r_la	= H_R_L_O_A + B_1_A*(T_r - T_R_O_A) + (B_2_A/2)*( pow(T_r,2) - pow(T_R_O_A,2) ) + (B_3_A/3)*( pow(T_r,3) - pow(T_R_O_A,3) ) - A_1_A*P_R_O_A - (A_2_A/2)*( pow(P_r,2) + pow(P_R_O_A,2) ) - A_3_A*T_r*P_r - A_4_A*pow(T_r,2)*(2*P_r - P_R_O_A);

	//Return internal energy
	return ((U_r_la * R_U * T_B / MW_A) - U_REF_L_A);
}

//Evaluates the density of liquid ammonia
double Rho_la (double T, double P)
{
	return (1.0/V_la(T,P));
}

//Evaluates the entropy of liquid ammonia
double S_la (double T, double P)
{
	//Define variable
	double		S_r_la;		//Reduced entropy of liquid ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced entropy of liquid ammonia
	S_r_la	= S_R_L_O_A + B_1_A*log(T_r / T_R_O_A) + B_2_A*( T_r - T_R_O_A ) + (B_3_A/2)*( pow(T_r,2) - pow(T_R_O_A,2) ) - A_3_A*(P_r - P_R_O_A) - 2*A_4_A*T_r*(P_r - P_R_O_A);

	//Return entropy
	return ((S_r_la * R_U / MW_A) - S_REF_L_A);
}

//Evaluates the specific heat at constant pressure of liquid ammonia
double Cp_la (double T, double P)
{
	//Define variable
	double		Cp_r_la;	//Reduced Cp of liquid ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced Cp of liquid ammonia
	Cp_r_la	= B_1_A + B_2_A*T_r + B_3_A*pow(T_r,2) - 2*A_4_A*T_r*(P_r - P_R_O_A);

	//Return Cp
	return ( Cp_r_la * R_U / MW_A );
}

//Chemical potential of liquid ammonia
double ChemPo_la (double T, double P, double x)
{
	double g_l;					//Gibbs energy of liquid mixture
	double T_r, P_r;			//Reduced quantities
	double x_mo;				//Molar ammonia fraction
	double MW_M;				//Molecular weight of the mixture
	double f_1, f_2, f_3;		//Components of excess enthalpy
	double h_la, h_lw;			//Enthalpies of each component
	double s_la, s_lw;			//Entropy of the components
	double dhedx, dsedx, dsmdx;	//Derivatives of the excess enthalpy and entropy, and entropy of mixing
	double mu;					//Chemical potential
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculate molar quantities
	x_mo	= (x/MW_A)/ ( (x/MW_A) + ((1-x)/MW_W) );
	MW_M	= x*MW_A + (1-x)*MW_W;
	
	//Gibbs energy of the mixture
	g_l		= G_lm(T, P, x)*MW_M;
	
	//Enthalpies of the components
	h_lw	= H_lw(T, P)*MW_W;
	h_la	= H_la(T, P)*MW_A;
	
	//Calculate excess enthalpy derivative
	f_1		= E_1 + E_2*P_r + 2*E_5/T_r + 3*E_6/pow(T_r,2);
	f_2		= E_7 + E_8*P_r + 2*E_11/T_r + 3*E_12/pow(T_r,2);
	f_3		= E_13 + E_14*P_r + 2*E_15/T_r + 3*E_16/pow(T_r, 2);
	dhedx	= R_U*T_B*( f_1*(-2*x_mo + 1) + f_2*(-6*pow(x_mo,2) + 6*x_mo - 1) + f_3*(-16*pow(x_mo,3) + 24*pow(x_mo,2) - 10*x_mo + 1) );
	
	//Calculate entropies
	s_lw	= S_lw(T, P)*MW_W;
	s_la	= S_la(T, P)*MW_A;
	
	//Calculate excess entropy derivative
	f_1		= -(E_3 + E_4*P_r) +E_5/pow(T_r,2) + 2*E_6/pow(T_r,3);
	f_2		= -(E_9 + E_10*P_r) + E_11/pow(T_r,2) + 2*E_12/pow(T_r,3);
	f_3		= E_15/pow(T_r,2) + 2*E_16/pow(T_r,3);
	dsedx	= R_U*( f_1*(-2*x_mo + 1) + f_2*(-6*pow(x_mo,2) + 6*x_mo - 1) + f_3*(-16*pow(x_mo,3) + 24*pow(x_mo,2) - 10*x_mo + 1) );
	
	//Entropy of mixing
	if (x_mo == 0)
	{	dsmdx = 1E6;}	//Some big number
	else if (x_mo == 1)
	{	dsmdx = -1E6;}	//Some small number
	else
	{	dsmdx	= -R_U*log(x_mo / (1-x_mo));}
	
	
	//return chemical potential
	mu		= g_l + (1-x_mo)*(h_la - h_lw + dhedx - T*( s_la - s_lw + dsedx + dsmdx ) );
	return (mu/MW_A);
}


//Evaluates the Gibbs energy of vapor ammonia
double G_va (double T, double P)
{
	//Define variable
	double		G_r_va;		//Reduced Gibbs energy of vapor ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculation of reduced Gibbs energy
	G_r_va		= H_R_V_O_A - T_r*S_R_V_O_A + D_1_A*(T_r - T_R_O_A) + (D_2_A/2)*( pow(T_r,2) - pow(T_R_O_A,2) ) + (D_3_A/3)*( pow(T_r,3) - pow(T_R_O_A,3) ) - D_1_A*T_r*log(T_r/T_R_O_A) - D_2_A*T_r*(T_r - T_R_O_A) - (D_3_A/2)*T_r*( pow(T_r,2) - pow(T_R_O_A,2) ) + T_r*log(P_r/P_R_O_A) + C_1_A*(P_r-P_R_O_A) + C_2_A*( P_r/pow(T_r,3) - 4*P_R_O_A/pow(T_R_O_A,3) + 3*P_R_O_A*T_r/pow(T_R_O_A,4) ) + C_3_A*(P_r/pow(T_r,11) - 12*P_R_O_A/pow(T_R_O_A,11) + 11*P_R_O_A*T_r/pow(T_R_O_A,12) ) + (C_4_A/3)*( pow(P_r,3)/pow(T_r,11) - 12*pow(P_R_O_A,3)/pow(T_R_O_A,11) + 11*pow(P_R_O_A,3)*T_r/pow(T_R_O_A,12) );	
	
	//Return Gibbs energy
	return ((G_r_va * R_U * T_B / MW_A) - G_REF_V_A);
}


//Evaluates the Enthalpy of vapor ammonia
double H_va (double T, double P)
{
	//Define variable
	double		H_r_va;		//Reduced Enthalpy of vapor ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculation of reduced Enthalpy
	H_r_va		= H_R_V_O_A + (4*C_2_A)*(P_r/pow(T_r,3) - P_R_O_A/pow(T_R_O_A,3) ) + (D_2_A/2)*( pow(T_r,2) - pow(T_R_O_A,2) ) + (D_3_A/3)*( pow(T_r,3) - pow(T_R_O_A,3) ) + 12*C_3_A*(P_r/pow(T_r,11) - P_R_O_A/pow(T_R_O_A,11) ) + 4*C_4_A*( pow(P_r,3)/pow(T_r,11) - pow(P_R_O_A,3)/pow(T_R_O_A,11) )  + D_1_A*(T_r - T_R_O_A) + C_1_A*(P_r - P_R_O_A);
	
	
	//Return Enthalpy
	return ((H_r_va * R_U * T_B / MW_A) - H_REF_V_A);
}

//Evaluates the specific volume of vapor ammonia
double V_va (double T, double P)
{
	//Define variable
	double		V_r_va;		//Reduced volume of vapor ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced volume of vapor ammonia
	V_r_va	= T_r/P_r + C_1_A + C_2_A/pow(T_r,3) + C_3_A/pow(T_r,11) + C_4_A*pow(P_r,2)/pow(T_r,11);

	//Return specific volume
	return ( V_r_va * R_U * T_B / (P_B * MW_A) );
}

//Evaluates the density of vapor ammonia
double Rho_va (double T, double P)
{
	return (1.0/V_va(T,P));
}

//Evaluates the entropy of vapor ammonia
double S_va (double T, double P)
{
	//Define variable
	double		S_r_va;		//Reduced entropy of vapor ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced entropy of vapor ammonia
	S_r_va	= S_R_V_O_A + D_1_A*log(T_r/T_R_O_A) + D_2_A*(T_r - T_R_O_A) + (D_3_A/2)*( pow(T_r,2) - pow(T_R_O_A,2) ) - log( P_r / P_R_O_A ) + (3*C_2_A)*( P_r/pow(T_r,4) - P_R_O_A/pow(T_R_O_A,4) ) + (11*C_3_A)*(P_r/pow(T_r,12) - P_R_O_A/pow(T_R_O_A,12) ) + (11/3)*C_4_A*( pow(P_r,3)/pow(T_r,12) - pow(P_R_O_A,3)/pow(T_R_O_A,12) );	
	
	
	//Return vapor
	return ((S_r_va * R_U / MW_A) - S_REF_V_A);
}

//Evaluates the internal energy of vapor ammonia
double U_va (double T, double P)
{
	//Define variable
	double		U_r_va;		//Reduced energy of vapor ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced internal energy of vapor ammonia
	U_r_va	= H_R_V_O_A + (C_2_A )*(3*P_r/pow(T_r,3) - 4*P_R_O_A/pow(T_R_O_A,3) ) + (D_2_A/2)*( pow(T_r,2) - pow(T_R_O_A,2) ) + (D_3_A/3)*( pow(T_r,3) - pow(T_R_O_A,3) ) + C_3_A*(11*P_r/pow(T_r,11) - 12*P_R_O_A/pow(T_R_O_A,11) ) + C_4_A*( 3*pow(P_r,3)/pow(T_r,11) - 4*pow(P_R_O_A,3)/pow(T_R_O_A,11) ) + D_1_A*(T_r - T_R_O_A) - C_1_A*P_R_O_A - T_r;
	
	
	//Return internal energy
	return ((U_r_va * R_U * T_B / MW_A) - U_REF_V_A);
}

//Evaluates the specific heat at constant pressure of vapor ammonia
double Cp_va (double T, double P)
{
	//Define variable
	double		Cp_r_va;	//Reduced Cp of liquid ammonia
	double		T_r, P_r;	//Reduced quantities
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Reduced Cp of vapor ammonia
	Cp_r_va	= D_1_A + D_2_A*T_r + D_3_A*pow(T_r,2) - 12*C_2_A*P_r/pow(T_r,4) - (11*P_r/pow(T_r,12) )*( 12*C_3_A + 4*C_4_A*pow(P_r,2) );

	//Return Cp
	return ( Cp_r_va * R_U / MW_A );
}

//Chemical potential of vapor ammonia
double ChemPo_va (double T, double P, double y)
{
	double g_v;					//Gibbs energy of vapor mixture
	double y_mo;				//Molar ammonia fraction
	double MW_M;				//Molecular weight of the mixture
	double h_va, h_vw;			//Enthalpies of each component
	double s_va, s_vw;			//Entropy of the components
	double dsmdx;	//Derivatives of the excess enthalpy and entropy, and entropy of mixing
	double mu;					//Chemical potential
	
	//Calculate molar quantities
	y_mo	= (y/MW_A)/ ( (y/MW_A) + ((1-y)/MW_W) );
	MW_M	= y*MW_A + (1-y)*MW_W;
	
	//Gibbs energy of the mixture
	g_v		= G_vm(T, P, y)*MW_M;
	
	//Enthalpies of the components
	h_vw	= H_vw(T, P)*MW_W;
	h_va	= H_va(T, P)*MW_A;
	
	//Calculate entropies
	s_vw	= S_vw(T, P)*MW_W;
	s_va	= S_va(T, P)*MW_A;
	
	//Entropy of mixing
	if (y_mo == 0)
	{	dsmdx = 1E6;}	//Some big number
	else if (y_mo == 1)
	{	dsmdx = -1E6;}	//Some small number
	else
	{	dsmdx	= -R_U*log(y_mo / (1-y_mo));}
	
	
	//return chemical potential
	mu		= g_v + (1-y_mo)*(h_va - h_vw - T*( s_va - s_vw + dsmdx ) );
	//return mu;
	return (mu/MW_A);
}

//Saturation pressure of ammonia (T in K, Psat in Pa)
double Psat_a (double T)
{
	//Define some normalize variables
	double Tr	= T/TC_A;						//Reduced temperature
	double Ttr	= TT_A/TC_A;					//Reduced triple point temperature
	double t	= (Tr - Ttr)/(1 - Ttr);			//Normalized temperature
	double Ptr	= PT_A/PC_A;					//Reduced triple pressure
	double tb	= ((TB_A/TC_A) - Ttr)/(1-Ttr);	//Normalized boiling temperature
	double Psatr;								//Reduced saturation pressure
	
	//Calculate reduced saturation pressure
	Psatr = exp( (Ttr*log(Ptr))*(1-t) /( (t*(1-Ttr) + Ttr)*(1+F_1_A*t)*(1 + F_2_A*t*(t-tb)) ) );

	//Return saturation pressure
	return (Psatr * PC_A);
}

//==================================================================================================
//VAPOR MIXTURES PROPERTIES BELOW
//==================================================================================================
//Enthalpy of the vapor mixture
double H_vm (double T, double P, double y)
{
	double h_va, h_vw;
	//get these quantities from the functions already defined
	h_va	= H_va(T, P);
	h_vw	= H_vw(T, P);
	//return ideal mixture
	return (y*h_va + (1-y)*h_vw);
}

//Specific volume of the vapor mixture
double V_vm (double T, double P, double y)
{
	double v_va, v_vw;
	//get these quantities from the functions already defined
	v_va	= V_va(T, P);
	v_vw	= V_vw(T, P);
	//return ideal mixture
	return (y*v_va + (1-y)*v_vw);
}

//Density of the vapor mixture
double Rho_vm (double T, double P, double y)
{
	//comes directly from specific volume
	return (1/V_vm(T,P,y));
}

//Entropy of the vapor mixture
double S_vm (double T, double P, double y)
{
	double s_va, s_vw, s_mix;
	double y_mo;		//Molar fraction of ammonia
	double MW_M;		//Molar weight of the mixture
	
	//get these quantities from the functions already defined
	s_va	= S_va(T, P);
	s_vw	= S_vw(T, P);
	
	//Evaluate the entropy of mixing, double check this formulation
	y_mo	= (y/MW_A)/ ( (y/MW_A) + ((1-y)/MW_W) );
	MW_M	= y*MW_A + (1-y)*MW_W;
	if ( (y > 0) && (y < 1) )
	{	s_mix	= -R_U*(y_mo*log(y_mo) + (1-y_mo)*log(1-y_mo))/MW_M; }
	else
	{	s_mix 	= 0;}
	
	//return ideal mixture
	return (y*s_va + (1-y)*s_vw + s_mix);
}

//Internal energy of the vapor mixture
double U_vm (double T, double P, double y)
{
	double h_vm, v_vm;
	//get these quantities from the functions already defined
	h_vm	= H_vm(T, P, y);
	v_vm	= V_vm(T, P, y);
	//return ideal mixture
	return (h_vm - P*v_vm);
}

//GIbbs energy of the vapor mixture
double G_vm (double T, double P, double y)
{
	double h_vm, s_vm;
	//get these quantities from the functions already defined
	h_vm	= H_vm(T, P, y);
	s_vm	= S_vm(T, P, y);
	//return ideal mixture
	return (h_vm - T*s_vm);
}

//Specific heat of the vapor mixture
double Cp_vm (double T, double P, double y)
{
	double cp_va, cp_vw;
	//get these quantities from the functions already defined
	cp_va	= Cp_va(T, P);
	cp_vw	= Cp_vw(T, P);
	//return ideal mixture
	return (y*cp_va + (1-y)*cp_vw);
}

//Saturation temperature of ammonia (P in Pa, Tsat in K)
double Tsat_a (double P)
{
	//Define some normalize variables
	double Tr;									//Reduced temperature
	double Ttr	= TT_A/TC_A;					//Reduced triple point temperature
	double t;									//Normalized temperature
	double Ptr	= PT_A/PC_A;					//Reduced triple pressure
	double tb	= ((TB_A/TC_A) - Ttr)/(1-Ttr);	//Normalized boiling temperature
	double Psatr = P/PC_A;						//Reduced saturation pressure
	
	//iteration stuff
	double	eps		= 1E-6;						//convergence criterion
	double	err		= 1;						//Current error
	double	R;									//Current residual
	double	J;									//Current derivative of the residual
	double	dt;									//Change in normalized temperature
	
	//Initial guess for t
	t	= Psatr;
	
	//Iteration
	while ( err > eps )
	{
		//Calculate residual
		R	= log(Psatr) - 	(Ttr*log(Ptr))*(1-t) /( (t*(1-Ttr) + Ttr)*(1+F_1_A*t)*(1 + F_2_A*t*(t-tb)) );
		
		//Calculate derivative of residual
		J	= Ttr*log(Ptr)*( 1 + 2*F_1_A*t + 2*F_2_A*t - F_2_A*tb - 2*Ttr*F_1_A*t + F_1_A*F_2_A*tb*pow(t,2) + 2*F_1_A*F_2_A*Ttr*pow(t,3) - 4*F_1_A*F_2_A*Ttr*pow(t,2) - 2*F_1_A*F_2_A*tb*t + 2*F_1_A*F_2_A*Ttr*t - F_1_A*F_2_A*Ttr*tb - F_1_A*pow(t,2) - F_2_A*pow(t,2) + (F_1_A+F_2_A)*Ttr + F_1_A*Ttr*pow(t,2) - 2*F_1_A*F_2_A*pow(t,3) + F_2_A*Ttr*pow(t,2) - 2*F_2_A*Ttr*t + 3*F_1_A*F_2_A*pow(t,2) - F_1_A*F_2_A*Ttr*tb*pow(t,2) + 2*F_1_A*F_2_A*Ttr*tb*t )/ pow( (-t + t*Ttr - Ttr)*(1 + F_1_A*t)*(1 + F_2_A*t - F_2_A*tb), 2);
	
		//calculate change in t
		dt	= -R/J;
		t	+= dt;
	
		//update error
		err = dt;
	}

	//Un-normalize T
	Tr		= t*(1-Ttr) + Ttr;

	//Return saturation temperature
	return (Tr * TC_A);
}

//Dew temperature for a vapor mixture, from Patek and Klomfar, 1995
double Tdew	(double P, double y)
{
	//Variables
	int		i;			//Counter
	double	T = 0;		//Temperature
	double	y_mo;		//Molar concentration
	
	//Constants from Patek
	double	T0	= 100;	//Ref temperature [K]
	double	P0	= 2E6;	//Reference pressure [MPa]
	double	m[]	= {0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.50, 0.50, 0.75, 0.75, 1.00, 1.00, 1.25, 1.25, 1.50, 1.75};
	int		n[]	= {0, 1, 2, 3, 0, 1, 2, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2};
	double 	a[]	= {0.324004E1, -0.395920E0, 0.435624E-1, -0.218943E-2, -0.143526E1, 0.105256E1, -0.719281E-1, 0.122362E2, -0.224368E1, -0.201780E2, 0.110834E1, 0.145399E2, 0.644312E0, -0.221246E1, -0.756266E0, -0.135529E1, 0.183541E0};

	//Get molar concentration
	y_mo	= (y/MW_A)/ ( (y/MW_A) + ((1-y)/MW_W) );
	
	//High order approximation for the saturation curve
	for (i = 0; i < 17; i++)
		T	+= a[i]*pow( 1-y_mo, m[i] )*pow( log(P0/P), n[i] );
	
	T	= T*T0;		//Un-normalize temperature
	return T;		//Return saturation temperature
}

//Gets vapor concentration at dew point
double ydew		(double T, double P)
{
	//Declare variables
	double 	Ps_a, Ps_w;		//Saturation temperatures for Ammonia + Water
	double 	y;				//Dew saturation concentration
	//Iteration
	double	err = 1;		//Current error
	double	eps = 1E-6;		//Target error
	int		i = 0;			//Iteration counter
	double 	deltay;			//Change in y for initial calculation of derivative
	double	dy;				//Change in y each iteration
	double	Td;				//current dew temperature
	double	yp, Tp;			//Previous concentration and dew temperature
	double	J, Jinv, R;		//Jacobian, its inverse, and the current residual
	
	//Get single component saturation points.
	Ps_a		= Psat_a (T);
	Ps_w		= Psat_w (T);
	//Initial guess for saturation concentration
	y			= (Ps_w - P)/(Ps_w - Ps_a);
	//enforce y in [0,1]
	y	= median(1E-4, y, 1-1E-4);
	
	//Iterate to find dew concentration
	while ( (err > eps) && ( i < 10 ) ) 
	{
		//Calculate new dew temperature
		Td		= Tdew(P, y);
		//Update residual
		R		= T - Td;		
		//Get Jacobian and inverse
		if (i == 0)		//First run, can't use previous result to estimate derivative
		{	deltay  = min_val( (1E-4)*(1-y), (1E-4)*y );		//for calculating derivative of Residual
			J		= -( Tdew(P, y+deltay) - Td )/deltay;	}
		else			//Use previous result to estimate the derivative
		{	J		= -(Td - Tp)/(y - yp); }
		Jinv		= 1/J;
		//update previous iteration values, for calculating residual
		Tp			= Td;
		yp			= y;
		//Get change in y
		dy			= -Jinv*R;
		//Apply limiter
		if (dy > 0) 
		{	dy = min_val(dy, 0.9*(1-y) ); }
		else
		{	dy = max_val(dy, -0.9*y); }
		//Update y
		y			+= dy;
		//Update error
		err			= fabs(dy);
		//Update iteration counter
		i++;
	}

	//enforce y in [0,1]
	y	= median(0, y, 1);

	//Return saturated vapor concentration
	return y;
}

//Get internal energy at the dew point
double Pdew	 (double y, double v)
{
	//Declare variables
	double 	P;				//Dew saturation pressure
	//Iteration
	double	err = 1;		//Current error
	double	eps = 1E-6;		//Target error
	int		i = 0;			//Iteration counter
	double 	deltaP;			//Change in pressure for initial calculation of derivative
	double	dP;				//Change in pressure each iteration
	double	Td;				//current dew temperature
	double	vd;				//current dew volume
	double	yp, Tp, vp, Pp;	//Previous concentration and dew point temperature
	double	J, Jinv, R;		//Jacobian, its inverse, and the current residual
	
	//Initial guess for saturation pressure
	P			= 1E6;
	//enforce y in [0,1]
	y	= median(1E-6, y, 1-1E-6);
	
	//Iterate to find bubble point concentration
	while ( (err > eps) && (i < 10) )
	{
		//Calculate new bubble point temperature and volume
		Td		= Tdew(P, y);
		vd		= V_vm(Td, P, y);
		//Update residual
		R		= v - vd;		
		//Get Jacobian and inverse
		if (i == 0)		//First run, can't use previous result to estimate derivative
		{	deltaP  = 1;		//for calculating derivative of Residual
			J		= -( V_vm( Tdew(P+deltaP, y), P+deltaP, y ) - vd )/deltaP;	}
		else			//Use previous result to estimate the derivative
		{	J		= -(v - vp)/(P - Pp); }
		Jinv		= 1/J;
		//update previous iteration values, for calculating residual
		Pp			= P;
		vp			= vd;
		Tp			= Td;
		yp			= y;
		//Get change in P
		dP			= -Jinv*R;
		//Apply limiter
		if (dP > 0) 
		{	dP = min_val(dP, 0.5*P ); }
		else
		{	dP = max_val(dP, -0.5*P); }
		//Update P
		P			+= dP;
		//Update error
		err			= fabs(dP/P);
		//Update iteration counter
		i++;
	}

	//Return saturated liquid concentration
	return P;
}

//Gets dew point pressure given T, Y
double Pdew_TY	(double T_i, double Y_i)
{
	double	P;			//Dew pressure
	double	P_w, P_a;	//For calculating initial guess
	double	T;			//Temporary placeholder
	//Iteration variables
	double	err = 1;		//Current error
	double	eps = 1E-6;		//Target error
	int		i = 0;			//Iteration counter
	double 	deltaP = 1;		//Change in P for initial calculation of derivative
	double	dP;				//Change in P each iteration
	double	Pp, Tp;			//Previous concentration and dew point temperature
	double	J, Jinv, R;		//Jacobian, its inverse, and the current residual
	
	//Get saturation pressures and make initial guess for P
	P_w		= Psat_w (T_i);
	P_a		= Psat_a (T_i);
	P		= Y_i*P_a + (1-Y_i)*P_w;
	
	//begin iteration to find P
	while ( (err > eps) && (i < 10) )
	{
		//Calculate new dew point temperature
		T		= Tdew(P, Y_i);
		//Update residual
		R		= T_i - T;		
		//Get Jacobian and inverse
		if (i == 0)		//First run, can't use previous result to estimate derivative
		{	J		= -( Tdew(P+deltaP, Y_i) - T ) / deltaP;	}
		else			//Use previous result to estimate the derivative
		{	J		= -( T - Tp )/( P - Pp); }
		Jinv		= 1/J;
		//update previous iteration values, for calculating residual
		Tp			= T;
		Pp			= P;
		//Get change in P
		dP			= -Jinv*R;
		//Apply limiter
		if (dP > 0.5*P) 
		{	dP = 0.5*P; }
		else if ( dP < (-0.5*P) )
		{	dP = -0.5*P; }
		//Update P
		P			+= dP;
		//Update error
		err			= fabs(dP/P);
		//Update iteration counter
		i++;
	}


	return P;
}

//==================================================================================================
//LIQUID MIXTURES PROPERTIES BELOW
//==================================================================================================
//Gibbs energy of liquid mixture
double G_lm (double T, double P, double x)
{
	double h_lm, s_lm;			//Components of Gibbs energy

	//get individual component Gibbs energies
	h_lm	= H_lm(T, P, x);
	s_lm	= S_lm(T, P, x);

	//return total Gibbs energy
	return (h_lm - T*s_lm);
}

//Enthalpy of liquid mixture
double H_lm (double T, double P, double x)
{
	double h_la, h_lw, h_em;	//enthalpy of the ammonia, water, and excess mixing
	double T_r, P_r;			//Reduced quantities
	double x_mo;				//Molar ammonia fraction
	double MW_M;				//Molecular weight of the mixture
	double f_1, f_2, f_3;		//Components of excess enthalpy

	//get individual component enthalpies
	h_la	= H_la(T, P);
	h_lw	= H_lw(T, P);

	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculate molar quantities
	x_mo	= (x/MW_A)/ ( (x/MW_A) + ((1-x)/MW_W) );
	MW_M	= x*MW_A + (1-x)*MW_W;
	
	//Calculate excess enthalpy
	f_1		= E_1 + E_2*P_r + 2*E_5/T_r + 3*E_6/pow(T_r,2);
	f_2		= E_7 + E_8*P_r + 2*E_11/T_r + 3*E_12/pow(T_r,2);
	f_3		= E_13 + E_14*P_r + 2*E_15/T_r + 3*E_16/pow(T_r, 2);
	h_em	= (R_U*T_B/MW_M)*( f_1 + f_2*(2*x_mo - 1) + f_3*pow((2*x_mo - 1),2) )*x_mo*( 1-x_mo );
	
	//return total enthalpy
	return (x_mo*h_la + (1-x_mo)*h_lw + h_em);
}

//Specific volume of liquid mixture
double V_lm (double T, double P, double x)
{
	double v_la, v_lw, v_em;	//specific volume of the ammonia, water, and excess mixing
	double T_r, P_r;			//Reduced quantities
	double x_mo;				//Molar ammonia fraction
	double MW_M;				//Molecular weight of the mixture
	double f_1, f_2, f_3;		//Components of excess volume

	//get individual component specific volumes
	v_la	= V_la(T, P);
	v_lw	= V_lw(T, P);

	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculate molar quantities
	x_mo	= (x/MW_A)/ ( (x/MW_A) + ((1-x)/MW_W) );
	MW_M	= x*MW_A + (1-x)*MW_W;
	
	//Calculate excess specific volume
	f_1		= E_2 + E_4*T_r;
	f_2		= E_8 + E_10*T_r;
	f_3		= E_14;
	v_em	= (R_U*T_B/(P_B*MW_M))*( f_1 + f_2*(2*x_mo - 1) + f_3*pow((2*x_mo - 1),2) )*x_mo*( 1-x_mo );
	//return total specific volume
	return (x_mo*v_la + (1-x_mo)*v_lw + v_em);
}

//Density of the liquid mixture
double Rho_lm (double T, double P, double x)
{
	//comes directly from specific volume
	return (1/V_lm(T,P,x));
}

//Entropy of liquid mixture
double S_lm (double T, double P, double x)
{
	double s_la, s_lw, s_mix, s_em;	//specific volume of the ammonia, water, and excess mixing
	double T_r, P_r;				//Reduced quantities
	double x_mo;					//Molar ammonia fraction
	double MW_M;					//Molecular weight of the mixture
	double f_1, f_2, f_3;			//Components of excess volume

	//get individual component entropies
	s_la	= S_la(T, P);
	s_lw	= S_lw(T, P);

	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculate molar quantities
	x_mo	= (x/MW_A)/ ( (x/MW_A) + ((1-x)/MW_W) );
	MW_M	= x*MW_A + (1-x)*MW_W;
	
	//Evaluate the entropy of mixing
	if ( (x > 0) && (x < 1) )
	{	s_mix	= -(R_U/MW_M)*(x_mo*log(x_mo) + (1-x_mo)*log(1-x_mo));}
	else
	{	s_mix	= 0;}
	
	//Calculate excess entropies
	f_1		= -(E_3 + E_4*P_r) +E_5/pow(T_r,2) + 2*E_6/pow(T_r,3);
	f_2		= -(E_9 + E_10*P_r) + E_11/pow(T_r,2) + 2*E_12/pow(T_r,3);
	f_3		= E_15/pow(T_r,2) + 2*E_16/pow(T_r,3);
	s_em	= (R_U/MW_M)*( f_1 + f_2*(2*x_mo - 1) + f_3*pow((2*x_mo - 1),2) )*x_mo*( 1-x_mo );
	
	//return total specific volume
	return (x_mo*s_la + (1-x_mo)*s_lw + s_mix + s_em);
}

//Internal energy of the liquid mixture
double U_lm (double T, double P, double x)
{
	double h_lm, v_lm;
	//get these quantities from the functions already defined
	h_lm	= H_lm(T, P, x);
	v_lm	= V_lm(T, P, x);
	//return internal energy
	return (h_lm - P*v_lm);
}

//Specific heat of the liquid mixture
double Cp_lm (double T, double P, double x)
{
	double cp_la, cp_lw;
	//get these quantities from the functions already defined
	cp_la	= Cp_la(T, P);
	cp_lw	= Cp_lw(T, P);
	//return ideal mixture
	return (x*cp_la + (1-x)*cp_lw);
}

//Bubble temperature for a liquid mixture, from Patek and Klomfar, 1995
double Tbub	(double P, double x)
{
	//Variables
	int		i;			//Counter
	double	T = 0;		//Temperature
	double	x_mo;		//Molar concentration
	
	//Constants from Patek
	double	T0	= 100;	//Ref temperature [K]
	double	P0	= 2E6;	//Reference pressure [MPa]
	int		m[]	= {0,0,0,0,0,1,1,1,2,4,5,5,6,13};
	int		n[]	= {0,1,2,3,4,0,1,2,3,0,0,1,0,1};
	double 	a[]	= {0.322302E1, -0.384206E0, 0.460965E-1, -0.378945E-2, 0.135610E-3, 0.487755E0, -0.120108E0, 0.106154E-1, -0.533589E-3, 0.785041E1, -0.115941E2, -0.523150E-1, 0.489596E1, 0.421059E-1};

	//Get molar concentration
	x_mo	= (x/MW_A)/ ( (x/MW_A) + ((1-x)/MW_W) );
	
	//High order approximation for the saturation curve
	for (i = 0; i < 14; i++)
		T	+= a[i]*pow( 1-x_mo, m[i] )*pow( log(P0/P), n[i] );
	
	T	= T*T0;		//Un-normalize temperature
	return T;		//Return saturation temperature
}

//Gets vapor concentration at dew point
double xbub		(double T, double P)
{
	//Declare variables
	double 	Ps_a, Ps_w;		//Saturation pressures for Ammonia + Water
	double 	x;				//Bubble saturation concentration
	//Iteration
	double	err = 1;		//Current error
	double	eps = 1E-6;		//Target error
	int		i = 0;			//Iteration counter
	double 	deltax;			//Change in x for initial calculation of derivative
	double	dx;				//Change in x each iteration
	double	Tb;				//current bubble temperature
	double	xp, Tp;			//Previous concentration and bubble point temperature
	double	J, Jinv, R;		//Jacobian, its inverse, and the current residual
	
	//Get single component saturation points.
	Ps_a		= Psat_a (T);
	Ps_w		= Psat_w (T);
	//Initial guess for saturation concentration
	x			= (Ps_w - P)/(Ps_w - Ps_a);
	//enforce x in [0,1]
	x	= median(1E-4, x, 1-1E-4);
	
	//Iterate to find bubble point concentration
	while ( (err > eps) && (i < 10) )
	{
		//Calculate new bubble point temperature
		Tb		= Tbub(P, x);
		//Update residual
		R		= T - Tb;		
		//Get Jacobian and inverse
		if (i == 0)		//First run, can't use previous result to estimate derivative
		{	deltax  = min_val( (1E-4)*(1-x), (1E-4)*x );		//for calculating derivative of Residual
			J		= -( Tbub(P, x+deltax) - Tb )/deltax;	}
		else			//Use previous result to estimate the derivative
		{	J		= -(Tb - Tp)/(x - xp); }
		Jinv		= 1/J;
		//update previous iteration values, for calculating residual
		Tp			= Tb;
		xp			= x;
		//Get change in x
		dx			= -Jinv*R;
		//Apply limiter
		if (dx > 0) 
		{	dx = min_val(dx, 0.9*(1-x) ); }
		else
		{	dx = max_val(dx, -0.9*x); }
		//Update x
		x			+= dx;
		//Update error
		err			= fabs(dx);
		//Update iteration counter
		i++;
	}

	//enforce x in [0,1]
	x	= median(0, x, 1);

	//Return saturated liquid concentration
	return x;
}

//Get internal energy at the bubble point
double Pbub		(double x, double v)
{
	//Declare variables
	//double 	Ps_a, Ps_w;		//Saturation pressures for Ammonia + Water
	double 	P;				//Bubble saturation pressure
	//Iteration
	double	err = 1;		//Current error
	double	eps = 1E-6;		//Target error
	int		i = 0;			//Iteration counter
	double 	deltaP;			//Change in x for initial calculation of derivative
	double	dP;				//Change in x each iteration
	double	Tb;				//current bubble temperature
	double	vb;				//current bubble volume
	double	xp, Tp, vp, Pp;	//Previous concentration and bubble point temperature
	double	J, Jinv, R;		//Jacobian, its inverse, and the current residual
	
	//Get single component saturation points.
	//Ps_a		= Psat_a (T);
	//Ps_w		= Psat_w (T);
	//Initial guess for saturation pressure
	P			= 1E6;
	//enforce x in [0,1]
	x	= median(1E-6, x, 1-1E-6);
	
	//Iterate to find bubble point concentration
	while ( (err > eps) && (i < 10) )
	{
		//Calculate new bubble point temperature and volume
		Tb		= Tbub(P, x);
		vb		= V_lm(Tb, P, x);
		//Update residual
		R		= v - vb;		
		//Get Jacobian and inverse
		if (i == 0)		//First run, can't use previous result to estimate derivative
		{	deltaP  = 1;		//for calculating derivative of Residual
			J		= -( V_lm( Tbub(P+deltaP, x), P+deltaP, x ) - vb )/deltaP;	}
		else			//Use previous result to estimate the derivative
		{	J		= -(v - vp)/(P - Pp); }
		Jinv		= 1/J;
		//update previous iteration values, for calculating residual
		Pp			= P;
		vp			= vb;
		Tp			= Tb;
		xp			= x;
		//Get change in x
		dP			= -Jinv*R;
		//Apply limiter
		if (dP > 0) 
		{	dP = min_val(dP, 0.5*P ); }
		else
		{	dP = max_val(dP, -0.5*P); }
		//Update P
		P			+= dP;
		//Update error
		err			= fabs(dP/P);
		//Update iteration counter
		i++;
	}

	//Return saturated liquid concentration
	return P;
}

//Gets bubble point pressure given T, X
double Pbub_TX	(double T_i, double X_i)
{
	double	P;			//Bubble pressure
	double	P_w, P_a;	//For calculating initial guess
	double	T;			//Temporary placeholder
	//Iteration variables
	double	err = 1;		//Current error
	double	eps = 1E-6;		//Target error
	int		i = 0;			//Iteration counter
	double 	deltaP = 1;		//Change in P for initial calculation of derivative
	double	dP;				//Change in P each iteration
	double	Pp, Tp;			//Previous concentration and bubble point temperature
	double	J, Jinv, R;		//Jacobian, its inverse, and the current residual
	
	//Get saturation pressures and make initial guess for P
	P_w		= Psat_w (T_i);
	P_a		= Psat_a (T_i);
	P		= X_i*P_a + (1-X_i)*P_w;
	
	//begin iteration to find P
	while ( (err > eps) && (i < 10) )
	{
		//Calculate new bubble point temperature
		T		= Tbub(P, X_i);
		//Update residual
		R		= T_i - T;		
		//Get Jacobian and inverse
		if (i == 0)		//First run, can't use previous result to estimate derivative
		{	J		= -( Tbub(P+deltaP, X_i) - T ) / deltaP;	}
		else			//Use previous result to estimate the derivative
		{	J		= -( T - Tp )/( P - Pp); }
		Jinv		= 1/J;
		//update previous iteration values, for calculating residual
		Tp			= T;
		Pp			= P;
		//Get change in P
		dP			= -Jinv*R;
		//Apply limiter
		if (dP > 0.5*P) 
		{	dP = 0.5*P; }
		else if ( dP < (-0.5*P) )
		{	dP = -0.5*P; }
		//Update P
		P			+= dP;
		//Update error
		err			= fabs(dP/P);
		//Update iteration counter
		i++;
	}


	return P;
}

//==================================================================================================
//DERIVATIVES FOR PERFORMING ITERATIVE SOLUTION
//==================================================================================================

// Derivative of liquid ammonia chemical potential with respect to mass ammonia concentration
double d_mu_la_d_x (double T, double P, double x)
{
	double d_g_d_x;				//Derivative of Gibbs energy of liquid mixture
	double T_r, P_r;			//Reduced quantities
	double x_mo, d_x_d_x;		//Molar ammonia fraction and its derivative
	double MW_M;				//Molecular weight of the mixture
	double f_1, f_2, f_3;		//Components of excess enthalpy
	double h_la, h_lw;			//Enthalpies of each component
	double s_la, s_lw;			//Entropy of the components
	double dhedx, dsedx, dsmdx;	//Derivatives of the excess enthalpy and entropy, and entropy of mixing
	double ddheddx, ddseddx, ddsmddx;	//Second derivatives of the above
	double d_mu_d_x;			//Derivative of chemical potential
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculate molar quantities
	x_mo	= (x/MW_A)/ ( (x/MW_A) + ((1-x)/MW_W) );
	d_x_d_x	= MW_A*MW_W / pow( MW_A*x - MW_A - MW_W*x, 2 );
	MW_M	= x*MW_A + (1-x)*MW_W;
	
	//Enthalpies of the components
	h_lw	= H_lw(T, P)*MW_W;
	h_la	= H_la(T, P)*MW_A;
	
	//Calculate excess enthalpy derivative
	f_1		= E_1 + E_2*P_r + 2*E_5/T_r + 3*E_6/pow(T_r,2);
	f_2		= E_7 + E_8*P_r + 2*E_11/T_r + 3*E_12/pow(T_r,2);
	f_3		= E_13 + E_14*P_r + 2*E_15/T_r + 3*E_16/pow(T_r, 2);
	dhedx	= R_U*T_B*( f_1*(-2*x_mo + 1) + f_2*(-6*pow(x_mo,2) + 6*x_mo - 1) + f_3*(-16*pow(x_mo,3) + 24*pow(x_mo,2) - 10*x_mo + 1) );
	ddheddx	= R_U*T_B*( -2*f_1 + f_2*(-12*x_mo + 6) + f_3*(-48*pow(x_mo,2) + 48*x_mo - 10) );
	
	
	//Calculate entropies
	s_lw	= S_lw(T, P)*MW_W;
	s_la	= S_la(T, P)*MW_A;
	
	//Calculate excess entropy derivative
	f_1		= -(E_3 + E_4*P_r) +E_5/pow(T_r,2) + 2*E_6/pow(T_r,3);
	f_2		= -(E_9 + E_10*P_r) + E_11/pow(T_r,2) + 2*E_12/pow(T_r,3);
	f_3		= E_15/pow(T_r,2) + 2*E_16/pow(T_r,3);
	dsedx	= R_U*( f_1*(-2*x_mo + 1) + f_2*(-6*pow(x_mo,2) + 6*x_mo - 1) + f_3*(-16*pow(x_mo,3) + 24*pow(x_mo,2) - 10*x_mo + 1) );
	ddseddx	= R_U*( -2*f_1 + f_2*(-12*x_mo + 6) + f_3*(-48*pow(x_mo,2) + 48*x_mo - 10) );
	
	//Entropy of mixing
	if (x_mo == 0)
	{	dsmdx = 1E6;	//Some big number
		ddsmddx	= -1E6;	
		d_g_d_x	= -1E6; }
	else if (x_mo == 1)
	{	dsmdx = -1E6;	//Some small number
		ddsmddx	= -1E6; 
		d_g_d_x = 1E6; }
	else
	{	dsmdx	= -R_U*log(x_mo / (1-x_mo));
		ddsmddx	= -R_U/( x_mo*(1 - x_mo) );
		d_g_d_x	= h_la - h_lw + dhedx - T*( s_la - s_lw - R_U*( log(x_mo) - log(1-x_mo) ) + dsedx );
	}
	
	//return chemical potential
	d_mu_d_x	= d_g_d_x - ( h_la - h_lw + dhedx - T*(s_la - s_lw + dsedx + dsmdx) ) + (1-x_mo)*( ddheddx - T*(ddseddx + ddsmddx ));
	//Convert to per mass basis and return
	d_mu_d_x	*= d_x_d_x;
	return (d_mu_d_x/MW_A);
}

// Derivative of liquid water chemical potential with respect to mass ammonia concentration
double d_mu_lw_d_x (double T, double P, double x)
{
	double d_g_d_x;				//Derivative of Gibbs energy of liquid mixture
	double T_r, P_r;			//Reduced quantities
	double x_mo, d_x_d_x;		//Molar ammonia fraction and its derivative
	double MW_M;				//Molecular weight of the mixture
	double f_1, f_2, f_3;		//Components of excess enthalpy
	double h_la, h_lw;			//Enthalpies of each component
	double s_la, s_lw;			//Entropy of the components
	double dhedx, dsedx, dsmdx;	//Derivatives of the excess enthalpy and entropy, and entropy of mixing
	double ddheddx, ddseddx, ddsmddx;	//Second derivatives of the above
	double d_mu_d_x;			//Derivative of chemical potential
	
	//Calculate reduced variables
	T_r		= T / T_B;
	P_r		= P / P_B;
	
	//Calculate molar quantities
	x_mo	= (x/MW_A)/ ( (x/MW_A) + ((1-x)/MW_W) );
	d_x_d_x	= MW_A*MW_W / pow( MW_A*x - MW_A - MW_W*x, 2 );
	MW_M	= x*MW_A + (1-x)*MW_W;
	
	//Enthalpies of the components
	h_lw	= H_lw(T, P)*MW_W;
	h_la	= H_la(T, P)*MW_A;
	
	//Calculate excess enthalpy derivative
	f_1		= E_1 + E_2*P_r + 2*E_5/T_r + 3*E_6/pow(T_r,2);
	f_2		= E_7 + E_8*P_r + 2*E_11/T_r + 3*E_12/pow(T_r,2);
	f_3		= E_13 + E_14*P_r + 2*E_15/T_r + 3*E_16/pow(T_r, 2);
	dhedx	= R_U*T_B*( f_1*(-2*x_mo + 1) + f_2*(-6*pow(x_mo,2) + 6*x_mo - 1) + f_3*(-16*pow(x_mo,3) + 24*pow(x_mo,2) - 10*x_mo + 1) );
	ddheddx	= R_U*T_B*( -2*f_1 + f_2*(-12*x_mo + 6) + f_3*(-48*pow(x_mo,2) + 48*x_mo - 10) );
	
	
	//Calculate entropies
	s_lw	= S_lw(T, P)*MW_W;
	s_la	= S_la(T, P)*MW_A;
	
	//Calculate excess entropy derivative
	f_1		= -(E_3 + E_4*P_r) +E_5/pow(T_r,2) + 2*E_6/pow(T_r,3);
	f_2		= -(E_9 + E_10*P_r) + E_11/pow(T_r,2) + 2*E_12/pow(T_r,3);
	f_3		= E_15/pow(T_r,2) + 2*E_16/pow(T_r,3);
	dsedx	= R_U*( f_1*(-2*x_mo + 1) + f_2*(-6*pow(x_mo,2) + 6*x_mo - 1) + f_3*(-16*pow(x_mo,3) + 24*pow(x_mo,2) - 10*x_mo + 1) );
	ddseddx	= R_U*( -2*f_1 + f_2*(-12*x_mo + 6) + f_3*(-48*pow(x_mo,2) + 48*x_mo - 10) );
	
	//Entropy of mixing
	if (x_mo == 0)
	{	dsmdx = 1E6;	//Some big number
		ddsmddx	= -1E6;	
		d_g_d_x	= -1E6; }
	else if (x_mo == 1)
	{	dsmdx = -1E6;	//Some small number
		ddsmddx	= -1E6; 
		d_g_d_x = 1E6; }
	else
	{	dsmdx	= -R_U*log(x_mo / (1-x_mo));
		ddsmddx	= -R_U/( x_mo*(1 - x_mo) );
		d_g_d_x	= h_la - h_lw + dhedx - T*( s_la - s_lw - R_U*( log(x_mo) - log(1-x_mo) ) + dsedx );
	}
	
	//return chemical potential
	d_mu_d_x	= d_g_d_x - ( h_la - h_lw + dhedx - T*(s_la - s_lw + dsedx + dsmdx) ) - x_mo*( ddheddx - T*(ddseddx + ddsmddx ));
	//Convert to per mass basis and return
	d_mu_d_x	*= d_x_d_x;
	return (d_mu_d_x/MW_W);
}

// Derivative of vapor ammonia chemical potential with respect to mass ammonia concentration
double d_mu_va_d_y (double T, double P, double y)
{
	double y_mo;				//Molar ammonia fraction
	double MW_M;				//Molecular weight of the mixture
	double h_va, h_vw;			//Enthalpies of each component
	double s_va, s_vw;			//Entropy of the components
	double dsmdy, ddsmddy;		//Derivatives of the excess enthalpy and entropy, and entropy of mixing
	double d_g_d_y;				//Derivative of the Gibbs energy
	double d_y_d_y;				//Derivative of the molar concentration wrt the mass concentration
	double d_mu_d_y;			//Derivative of the chemical potential
	
	//Calculate molar quantities
	y_mo	= (y/MW_A)/ ( (y/MW_A) + ((1-y)/MW_W) );
	MW_M	= y*MW_A + (1-y)*MW_W;
	d_y_d_y	= MW_A*MW_W / pow( MW_A*y - MW_A - MW_W*y, 2 );
	
	//Enthalpies of the components
	h_vw	= H_vw(T, P)*MW_W;
	h_va	= H_va(T, P)*MW_A;
	
	//Calculate entropies
	s_vw	= S_vw(T, P)*MW_W;
	s_va	= S_va(T, P)*MW_A;
	
	//Entropy of mixing, and derivatives
	if (y_mo == 0)
	{	dsmdy 	=  1E6;	//Some big number
		ddsmddy	= -1E6;	
		d_g_d_y	= -1E6; }
	else if (y_mo == 1)
	{	dsmdy 	= -1E6;	//Some small number
		ddsmddy	= -1E6; 
		d_g_d_y = 1E6; }
	else
	{	dsmdy	= -R_U*log(y_mo / (1-y_mo));
		ddsmddy	= -R_U/( y_mo*(1 - y_mo) );
		d_g_d_y	= h_va - h_vw - T*( s_va - s_vw - R_U*( log(y_mo) - log(1-y_mo) ) ); }
	
	
	//return chemical potential derivative
	d_mu_d_y	= d_g_d_y - ( h_va - h_vw - T*(s_va - s_vw + dsmdy) ) - (1-y_mo)*T*ddsmddy;
	
	//return derivative, convert back to per mass basis
	d_mu_d_y	*= d_y_d_y;
	return (d_mu_d_y/MW_A);
}

// Derivative of vapor water chemical potential with respect to mass ammonia concentration
double d_mu_vw_d_y (double T, double P, double y)
{
	double y_mo;				//Molar ammonia fraction
	double MW_M;				//Molecular weight of the mixture
	double h_va, h_vw;			//Enthalpies of each component
	double s_va, s_vw;			//Entropy of the components
	double dsmdx, ddsmddx;		//Derivatives of the excess enthalpy and entropy, and entropy of mixing
	double d_g_d_y;				//Derivative of the Gibbs energy
	double d_y_d_y;				//Derivative of the molar concentration wrt the mass concentration
	double d_mu_d_y;			//Derivative of the chemical potential
	
	//Calculate molar quantities
	y_mo	= (y/MW_A)/ ( (y/MW_A) + ((1-y)/MW_W) );
	MW_M	= y*MW_A + (1-y)*MW_W;
	d_y_d_y	= MW_A*MW_W / pow( MW_A*y - MW_A - MW_W*y, 2 );
	
	//Enthalpies of the components
	h_vw	= H_vw(T, P)*MW_W;
	h_va	= H_va(T, P)*MW_A;
	
	//Calculate entropies
	s_vw	= S_vw(T, P)*MW_W;
	s_va	= S_va(T, P)*MW_A;
	
	//Entropy of mixing, and derivatives
	if (y_mo == 0)
	{	dsmdx 	=  1E6;	//Some big number
		ddsmddx	= -1E6;	
		d_g_d_y	= -1E6; }
	else if (y_mo == 1)
	{	dsmdx 	= -1E6;	//Some small number
		ddsmddx	= -1E6; 
		d_g_d_y = 1E6; }
	else
	{	dsmdx	= -R_U*log(y_mo / (1-y_mo));
		ddsmddx	= -R_U/( y_mo*(1 - y_mo) );
		d_g_d_y	= h_va - h_vw - T*( s_va - s_vw - R_U*d_y_d_y*( log(y_mo) - log(1-y_mo) ) ); }
	
	
	//return chemical potential derivative
	d_mu_d_y	= d_g_d_y - ( h_va - h_vw - T*(s_va - s_vw + dsmdx) ) + y_mo*T*ddsmddx;
	
	//return derivative, convert to per mass basis
	d_mu_d_y	*= d_y_d_y;
	return (d_mu_d_y/MW_W);
}


//==================================================================================================
//ITERATIVE SOLUTION SCHEMES
//==================================================================================================

//Solve for state given T, P, X
int Qxy_TPX (double *Q_o, double *x_o, double *y_o, double T_i, double P_i, double X_i)
{
	double	x, y;							//Liquid and vapor ammonia mass concentrations
	double	dx, dy;							//Increments of the above
	double	mu_lw, mu_la, mu_vw, mu_va;		//Chemical potentials
	double	J[2][2], Jinv[2][2];			//Matrices for the Jacobian, and inverse of Jacobian
	double	DetInv;							//Inverse determinant of J
	double	R[2];							//Residual at the instant
	double	eps = 1E-6;						//Convergence criterion
	double	err	= 1;						//Error at this iteration
	double	alpha;							//under-relaxation factor
	int		i = 0;							//iteration counter
	int		itsMax = 75;					//Iteration count limit
	
	//Initialize x and y
	//x	= 0.01;
	//y	= 0.99;
	x	= *x_o;
	y	= *y_o;
	
	//Get potentials and check initial error
	mu_lw		= ChemPo_lw (T_i, P_i, x);
	mu_la		= ChemPo_la (T_i, P_i, x);
	mu_vw		= ChemPo_vw (T_i, P_i, y);
	mu_va		= ChemPo_va (T_i, P_i, y);
	//We use the relative error in the chemical potentials as the uncertainty
	//Enforce 1 solver step
    //err			= sqrt( pow( mu_lw-mu_vw, 2) + pow( mu_la-mu_va, 2) ) / sqrt( pow( 0.5*(mu_lw+mu_vw), 2) + pow( 0.5*(mu_la+mu_va), 2) );
	err = 1;
    
	
	while ( (err > eps) && (i >= 0) )
	{
		//Calculate residual
		R[0]	= mu_lw - mu_vw;
		R[1]	= mu_la - mu_va;
		
		//Calculate Jacobian
		J[0][0]	=  d_mu_lw_d_x (T_i, P_i, x);
		J[0][1]	= -d_mu_vw_d_y (T_i, P_i, y);
		J[1][0]	=  d_mu_la_d_x (T_i, P_i, x);
		J[1][1]	= -d_mu_va_d_y (T_i, P_i, y);
		
		//Calculate inverse
		DetInv		= 1/ ( J[0][0]*J[1][1] - J[0][1]*J[1][0] );
		Jinv[0][0]	=  DetInv*J[1][1];
		Jinv[0][1]	= -DetInv*J[0][1];
		Jinv[1][0]	= -DetInv*J[1][0];
		Jinv[1][1]	=  DetInv*J[0][0];
		
		//Calculate increments for x,y
		dx	= -( Jinv[0][0]*R[0] + Jinv[0][1]*R[1] );
		dy	= -( Jinv[1][0]*R[0] + Jinv[1][1]*R[1] );
		
		//under-relax increments to stay well behaved
		alpha	= fabs( (double)( 0.5*(1-x)/dx ) );
		alpha	= min_val( alpha, fabs( 0.5*(1-y)/dy ) );
		alpha	= min_val( alpha, fabs( 0.5*x/dx ) );
		alpha	= min_val( alpha, fabs( 0.5*y/dy ) );
		alpha	= min_val( alpha, 1);
		
		//Update x and y
		x	+= alpha*dx;
		y	+= alpha*dy;
		
		//Update chemical potentials and error
		mu_lw		= ChemPo_lw (T_i, P_i, x);
		mu_la		= ChemPo_la (T_i, P_i, x);
		mu_vw		= ChemPo_vw (T_i, P_i, y);
		mu_va		= ChemPo_va (T_i, P_i, y);
		//err			= sqrt( pow( mu_lw-mu_vw, 2) + pow( mu_la-mu_va, 2) ) / sqrt( pow( 0.5*(mu_lw+mu_vw), 2) + pow( 0.5*(mu_la+mu_va), 2) );
		err	= alpha*sqrt(dx*dx + dy*dy);
		
		//increment iteration count
		i++;
		
		//Check iteration count limit
		if (i > itsMax)
			i = -1;
	}

    *Q_o	= ( (X_i-x)/(y-x) );
	*x_o	= x;
	*y_o	= y;
	
	//return iteration count
	return i;
}

//Solve for T, Q, x, y given P, X, H
int TQxy_PXH (double *T_o, double *Q_o, double *x_o, double *y_o, double P_i, double X_i, double H_i, double Tlim)
{
	double	T, Tp;							//Temperature and updated temperature
	double	Q;								//Quality
	double	x, y;							//Liquid and vapor ammonia mass concentrations
	double	dT;								//Temperature increment
	double	mu_lw, mu_la, mu_vw, mu_va;		//Chemical potentials
	double	H;								//Current enthalpy
	double	eps = 1E-6;						//Convergence criterion
	double	err	= 1;						//Error at this iteration
	double	errA, errB;						//Two different errors to calculate (chemical potential and enthalpy)
	double	alpha = 1;						//under-relaxation factor
	double	deltaT	= 1E-4;					//step in T for calculating derivatives
	double	Qp;								//Temporary quality
	double	errP = err;
	int	i = 0;								//iteration counter
	
	//new iteration scheme parameters
	double	dH_dT;
	double  Hp, xp, yp;
	double 	RH, JH, JinvH;
	
	double	Hm = 0;
	double	Tm = 0;							//Previous enthalpy, temperature
	
	//Initialize T, Q, x and y
	T	= *T_o;
	Tp	= T;
	x	= *x_o;
	y	= *y_o;
	Q	= ( (X_i-x)/(y-x) );
	
	//Get potentials for checking initial error
	mu_lw		= ChemPo_lw (T, P_i, x);
	mu_la		= ChemPo_la (T, P_i, x);
	mu_vw		= ChemPo_vw (T, P_i, y);
	mu_va		= ChemPo_va (T, P_i, y);
	H			= Q*H_vm(T, P_i, y) + (1-Q)*H_lm(T, P_i, x);
	//We use the relative error in the chemical potentials as the uncertainty
	errA		= sqrt( pow( mu_lw-mu_vw, 2) + pow( mu_la-mu_va, 2) ) / sqrt( pow( 0.5*(mu_lw+mu_vw), 2) + pow( 0.5*(mu_la+mu_va), 2) );
	errB		= fabs(H_i - H)/fabs(H_i + H);
	err			= max_val( errA, errB);
	errP 		= err;
	
	//New iteration scheme
	while ( fabs(err) > eps)
	{
		//First get residual
		Qxy_TPX (&Q, &x, &y, T, P_i, X_i);		//Get new Q,x,y
		Q	   = median( 0, Q, 1);
		x	   = (Q == 0) ? X_i : x;
		y	   = (Q == 1) ? X_i : y;
		//New enthalpy
		H = Q*H_vm(T, P_i, y) + (1-Q)*H_lm(T, P_i, x);
		//Residual for enthalpy
		RH = H_i - H;
		
		//Now get Jacobian (-dH/dT)
		if (i == 0)	//No previous data to use to estimate dH/dT
		{
			Tp     = T+deltaT; 											//New T for estimating dH_dT
			Qp     = Q; xp = x; yp = y;  								//Set guesses for new Q,x,y
			Qxy_TPX (&Qp, &xp, &yp, Tp, P_i, X_i);						//Get new Q,x,y
			Qp	   = median( 0, Qp, 1);
			xp	   = (Qp == 0) ? X_i : xp;
			yp	   = (Qp == 1) ? X_i : yp;
			Hp     = Qp*H_vm(Tp, P_i, yp) + (1-Qp)*H_lm(Tp, P_i, xp);	//Get Hp
			dH_dT  = (Hp - H)/deltaT;
		}
		else	//Use previous H, T pair to estimate derivative
		{
			dH_dT  = (H - Hm)/(T - Tm);
		}
		//Evaluate Jacobian
		JH     = -dH_dT;
		JinvH  = 1/JH;
		
		//Set previous values for future estimations of dH_dT
		Tm	   = T;
		Hm     = H;
		
		//Get delta T
		dT     = -JinvH*RH;
		dT     = (fabs(dT) > fabs(Tlim) ) ? fabs(Tlim)*sign(dT) : dT;				//Limiter
		
		//Update T
		T	   += alpha*dT;
		
		//Calculate error
		errP   = err;
		err    = dT / T;									//Update error
		
		//Avoid oscillations
		if ( sign(errP) != sign(err) )
		{	alpha	*= 0.8;	}
		
		//iteration counter
		i++;
	}
	
	//Return quality, x, and y
	*T_o	= T;
	*Q_o	= ( (X_i-x)/(y-x) );
	*x_o	= x;
	*y_o	= y;
	
	//return iteration count
	return i;
}

//Gets temperature and concentrations from PXQ
int Txy_PXQ (double *T_o, double *x_o, double *y_o, double P_i, double X_i, double Q_i, double Tlim)
{
    double  Tl, Tm, Th;
    double  Ql, Qm, Qh;
    double  xl = X_i/2.0;
    double  xm = X_i/2.0;
    double  xh = X_i/2.0;
    double  yl = 1.0 - (1.0 - X_i)/2.0;
    double  ym = 1.0 - (1.0 - X_i)/2.0;
    double  yh = 1.0 - (1.0 - X_i)/2.0;
    
    
	double	T, Tp;							//Concentration and updated concentration
	double	Q, Qp;							//Quality and previous quality
	double	x, y;							//Phasic concentrations
	double	dT;								//Concentration increment
	double	mu_lw, mu_la, mu_vw, mu_va;		//Chemical potentials
	double	eps = 1E-6;						//Convergence criterion
	double	err	= 1;						//Error at this iteration
	double	alpha = 1;						//under-relaxation factor
	double	deltaT	= 1.0;//1E-4;			//step in T for calculating derivatives
	double	errP = err;
	int	    i = 0;							//iteration counter
   	int	    itsMax = 50;					//Iteration count limit
    int     flag = 0;
    
	//new iteration scheme parameters
	double	dQ_dT;
	double  xp, yp;
	double 	R, J, Jinv;
	double  T_b, T_d;
	
//	double	Tm = 0;
//	double	Qm = 0;							//Previous enthalpy, temperature
	
	//Initialize T, Q, x and y
	T	= *T_o;
	Tp	= T;
	x	= *x_o;
	y	= *y_o;
	Q	= Q_i;
	Qp	= Q;
	
	//Apply limits to Q
    Q			= median(1E-8, Q, 1-1E-8);


    //Start with Bisection method
    mexPrintf("\n\nIterative solver start\n");

	//Get bubble and dew point temperature
    T_b = Tbub(P_i, X_i) - 2;
	T_d = Tdew(P_i, X_i) + 2;
	
	Tl = T_b;
	Ql = 0; xl = X_i; yl = 1.0 - 1E-6;
	
	Th  = T_d;
    Qh = 1; xh = 1E-6; yh = X_i;
        
    while ( (fabs(err) > eps) && (i < itsMax) )
    {
//mexPrintf("Tl = %f, Ql = %f.        Th = %f, Qh = %f\n", Tl, Ql, Th, Qh);
//mexPrintf("xl = %f, yl = %f.        xh = %f, yh = %f\n", xl, yl, xh, yh);
        //Tm = Tl + (Th - Tl)*(Q_i - Ql)/(Qh - Ql);
        Tm = (Tl+Th)/2.0;
        
		xm = (xh+xl)/2; ym = (yh+yl)/2;
		Qxy_TPX (&Qm, &xm, &ym, Tm, P_i, X_i);
//mexPrintf("    Qm = %f, xm = %f, ym = %f\n", Qm, xm, ym);
		if ((Qm < 0.0) || (Qm > 1.0))
		{
			if (Tm < (T_b+T_d)/2.0)
			{  Qm = 0.0; xm = X_i; ym = 1.0 - 1E-6;}
			else
			{  Qm = 1.0;  xm = 1E-6; ym = X_i;}
		}
		
		
//mexPrintf("    Qm = %f, xm = %f, ym = %f\n", Qm, xm, ym);
        if (Q >= Qm)
        {   Tl = Tm; xl = xm; yl = ym; Ql = Qm; }
        else if (Q < Qm)
        {   Th = Tm; xh = xm; yh = ym; Qh = Qm; }


        err = fabs( Qh - Ql );
        i++;  
    }

    if ( fabs(Qh - Q) < fabs(Q - Ql) )
    {  *T_o	= Th;  *x_o = xh;  *y_o = yh; }
    else
    {  *T_o	= Tl;  *x_o = xl;  *y_o = yl; }
    
//Old iteration scheme
/*
	while ( (fabs(err) > eps) && (i < itsMax) )
	{
		//First get residual
		flag = Qxy_TPX (&Q, &x, &y, T, P_i, X_i);		//Get new Q,x,y
//mexPrintf("flag = Qxy_TPX (&Q, &x, &y, T, P_i, X_i); flag = %d, Q = %f, x = %f, y = %f, T = %f, P_i = %f, X_i = %f\n", flag, Q, x, y, T, P_i, X_i);
		x	   = (Q <= 0) ? X_i : x;
		y	   = (Q >= 1) ? X_i : y;
		
		//Residual for concentration
		R = Q_i - Q;
//mexPrintf("R = Q_i - Q;, R = %f, Q_i = %f, Q = %f\n", R, Q_i, Q);
		//Now get Jacobian (-dQ/dX)
		if (i == 0)	//No previous data to use to estimate dH/dT
		{
			Tp     = T+deltaT; 								//New T for estimating dH_dT
			xp	= x; yp = y; Qp = (X_i - x)/(y - x);
//mexPrintf("Qp = (X_i - x)/(y - x); Q = %f\n", Qp);
			flag = Qxy_TPX (&Qp, &xp, &yp, Tp, P_i, X_i);			//Get new Q,x,y
//mexPrintf("flag = Qxy_TPX (&Qp, &xp, &yp, Tp, P_i, X_i); flag = %d, Qp = %f, xp = %f, yp = %f, Tp = %f, P_i = %f, X_i = %f\n", flag, Qp, xp, yp, Tp, P_i, X_i);
			dQ_dT  = (Qp - Q)/deltaT;
//mexPrintf("dQ_dT = (Qp - Q)/deltaT; Qp = %f, Q = %f, dQ_dT = %f\n", Qp, Q, dQ_dT);
		}
		//Sometimes derivative goes to 0 if changes in Q are too small, in this case, use previous value
		else if (Q != Qm)	//Use previous H, T pair to estimate derivative
		{
			dQ_dT  = (Q - Qm)/(T - Tm);
		}
		//Evaluate Jacobian
		J		= -dQ_dT;
		Jinv	= 1/J;
		
		//Set previous values for future estimations of dQ_dX
		Tm	   = T;
		Qm     = Q;
		
		//Get delta T
		dT     = -Jinv*R;
		
		//Apply limiter
		if (dT > 0) 
		{	dT = min_val(dT, Tlim ); }
		else
		{	dT = max_val(dT, -Tlim); }
		
		//Update T
		T	   += alpha*dT;
		//Calculate error
		errP   = err;
		err    = dT/T;									//Update error
		
		//Avoid oscillations
		if ( sign(errP) != sign(err) )
		{	alpha	*= 0.8;	}
		
		//Break if nothing is changing
		if (err == errP)
		{	err = 0; }
		
		//iteration counter
		i++;
	}
	
		//Send an error if the residual is still too high
	if ( fabs(R) > 1E-3 )
	{	i = -1; }

	//Return T, x, and y
	*T_o	= T;
	*x_o	= x;
	*y_o	= y;
*/	
	//return iteration count
	return i;
}

void T_PHx_l (double *T_o, double P_i, double H_i, double x_i)
{

	//General variables
	double	T = *T_o;			//get input temperature
	double	delta	= 1E-4;		//increment for temperature
	
	//For iterative solution
	double	dT;					//Change in T
	double	eps = 1E-6;			//Convergence criterion
	double	err	= 1;			//Current error in H
	double	R;					//Current residual
	double	J;					//Jacobian
	
	//iterative solve for T
	while (err > eps)
	{
		//Get residual
		R	= H_i - H_lm(T, P_i, x_i);
		
		//Get derivative
		J	= -( H_lm(T+(delta/2), P_i, x_i) - H_lm(T-(delta/2), P_i, x_i) )/delta;
		
		//calculate change in T
		dT	= -R/J;
		T	+= dT;
			
		//update error
		err = fabs( dT/T );
	}
	
	//Set T_o
	*T_o	= T;
}

void T_PHy_v (double *T_o, double P_i, double H_i, double y_i)
{

	//General variables
	double	T = *T_o;			//get input temperature
	double	delta	= 1E-4;		//increment for temperature
	
	//For iterative solution
	double	dT;					//Change in T
	double	eps = 1E-6;			//Convergence criterion
	double	err	= 1;			//Current error in H
	double	R;					//Current residual
	double	J;					//Jacobian
	
	//iterative solve for T
	while (err > eps)
	{
		//Get residual
		R	= H_i - H_vm(T, P_i, y_i);
		
		//Get derivative
		J	= -( H_vm(T+(delta/2), P_i, y_i) - H_vm(T-(delta/2), P_i, y_i) )/delta;
		
		//calculate change in T
		dT	= -R/J;
		T	+= dT;
			
		//update error
		err = fabs( dT/T );
	}
	
	//Set T_o
	*T_o	= T;
}

//Solve for T given P, Q, H, x, y
void T_PQHxy (double *T_o, double P_i, double Q_i, double H_i, double x_i, double y_i)
{
	//General variables
	double	T = *T_o;			//get input temperature
	double	delta	= 1E-4;		//increment for temperature
	
	//For iterative solution
	double	dT;					//Change in T
	double	eps = 1E-6;			//Convergence criterion
	double	err	= 1;			//Current error in H
	double	R;					//Current residual
	double	J;					//Jacobian
	
	//iterative solve for T
	while (err > eps)
	{
		//Get residual
		R	= H_i	- Q_i*H_vm(T, P_i, y_i) - (1-Q_i)*H_lm(T, P_i, x_i);
		
		//Get derivative
		J	= -Q_i * ( H_vm(T+(delta/2), P_i, y_i) - H_vm(T-(delta/2), P_i, y_i) )/delta - (1-Q_i) * ( H_lm(T+(delta/2), P_i, x_i) - H_lm(T-(delta/2), P_i, x_i) )/delta;
		
		//calculate change in T
		dT	= -R/J;
		T	+= dT;
			
		//update error
		err = fabs( dT/T );
	}
	
	//Set T_o
	*T_o	= T;
}

//Solve for T, P given Q, V, U, x & y
void TP_QVUxy (double *T_o, double *P_o, double Q_i, double V_i, double U_i, double x_i, double y_i)
{
	double	T, P;							//Liquid and vapor ammonia mass concentrations
	double	dT, dP;							//Increments of the above
	double	U, V;							//Current energy and volume
	double	UT, UP;							//Perturbed energies (for derivatives)
	double	VT, VP;							//Perturbed volumes (for derivatives)
	double	deltaT = 0.001;					//Change in temperature for evalulating derivatives
	double	deltaP = 0.1;					//Change in pressure for evalulating derivatives
	
	double	J[2][2], Jinv[2][2];			//Matrices for the Jacobian, and inverse of Jacobian
	double	DetInv;							//Inverse determinant of J
	double	R[2];							//Residual at the instant
	double	alpha, alphaP, alphaT;			//Under-relaxation factors
	double	eps = 1E-6;						//Convergence criterion
	double	err	= 1;						//Error at this iteration
	int		i = 0;							//iteration counter
	
	//Set initial values for T and P
	T		= *T_o;
	P		= *P_o;
	
	
	while (err > eps)
	{
		//Calculate volumes and internal energies
		if		( Q_i == 0 )
		{	U	= U_lm(T, P, x_i);
			UT	= U_lm(T+deltaT, P, x_i);
			UP	= U_lm(T, P+deltaP, x_i);
			V	= V_lm(T, P, x_i);
			VT	= V_lm(T+deltaT, P, x_i);
			VP	= V_lm(T, P+deltaP, x_i);	}
		else if	( Q_i == 1 )
		{	U	= U_vm(T, P, y_i);
			UT	= U_vm(T+deltaT, P, y_i);
			UP	= U_vm(T, P+deltaP, y_i);
			V	= V_vm(T, P, y_i);
			VT	= V_vm(T+deltaT, P, y_i);
			VP	= V_vm(T, P+deltaP, y_i);	}
		else
		{	U	= Q_i*U_vm(T, P, y_i) + (1-Q_i)*U_lm(T, P, x_i);
			UT	= Q_i*U_vm(T+deltaT, P, y_i) + (1-Q_i)*U_lm(T+deltaT, P, x_i);
			UP	= Q_i*U_vm(T, P+deltaP, y_i) + (1-Q_i)*U_lm(T, P+deltaP, x_i);
			V	= Q_i*V_vm(T, P, y_i) + (1-Q_i)*V_lm(T, P, x_i);
			VT	= Q_i*V_vm(T+deltaT, P, y_i) + (1-Q_i)*V_lm(T+deltaT, P, x_i);
			VP	= Q_i*V_vm(T, P+deltaP, y_i) + (1-Q_i)*V_lm(T, P+deltaP, x_i);	}
	
		//Calculate residual
		R[0]	= U_i - U;
		R[1]	= V_i - V;
		
		//Calculate Jacobian
		J[0][0]	= -(UT - U)/deltaT;
		J[0][1]	= -(UP - U)/deltaP;
		J[1][0]	= -(VT - V)/deltaT;
		J[1][1]	= -(VP - V)/deltaP;
		
		//Calculate inverse
		DetInv		= 1/ ( J[0][0]*J[1][1] - J[0][1]*J[1][0] );
		Jinv[0][0]	=  DetInv*J[1][1];
		Jinv[0][1]	= -DetInv*J[0][1];
		Jinv[1][0]	= -DetInv*J[1][0];
		Jinv[1][1]	=  DetInv*J[0][0];
		
		//Calculate increments for T, P
		dT	= -( Jinv[0][0]*R[0] + Jinv[0][1]*R[1] );
		dP	= -( Jinv[1][0]*R[0] + Jinv[1][1]*R[1] );
		
		//Under-relaxation, relative to magnitudes of T and P
		alphaT	= (fabs(dT) < 0.2*T) ? 1 : 0.2*T/fabs(dT);
		alphaP	= (fabs(dP) < 0.4*P) ? 1 : 0.4*P/fabs(dP);
		alpha	= min_val(alphaT, alphaP);
		
		//Update T and P
		T	+= alpha*dT;
		P	+= alpha*dP;
		
		//Update error
		err	= alpha*sqrt(dT*dT/(T*T) + dP*dP/(P*P) );
		
		//increment iteration count
		i++;
	}
	
	//Return Temperature and pressure
	*T_o	= T;
	*P_o	= P;
}

//Gets TPQxy from XVU
int TPQxy_XVU (double *T_o, double *P_o, double *Q_o, double *x_o, double *y_o, double X_i, double V_i, double U_i, double Tlim, double Plim)
{
	double	T, P;							//Liquid and vapor ammonia mass concentrations
	double	x, y, Q;						//Concentrations and qualities
	//Increments of temperature and pressure
	double	dT = 0;
	double 	dP = 0;							//Increments of the above
	double	U, V;							//Current energy and volume
	//Perturbed energies (for derivatives)
	double	UT = 0;
	double	UP = 0;
	//Perturbed volumes (for derivatives)
	double	VT = 0;
	double	VP = 0;
	double	QT, QP;							//Perturbed qualities
	double	xT, xP, yT, yP;					//Perturbed concentrations
	
	double	deltaT = 0.1*Tlim;					//Change in temperature for evalulating derivatives
	double	deltaP = 0.1*Plim;					//Change in pressure for evalulating derivatives
	
	double	J[2][2], Jinv[2][2];			//Matrices for the Jacobian, and inverse of Jacobian
	double	DetInv;							//Inverse determinant of J
	double	R[2];							//Residual at the instant
	double	alpha, alphaP, alphaT;			//Under-relaxation factors
	double	eps = 1E-6;						//Convergence criterion
	double	err	= 1;						//Error at this iteration
	int		i = 0;							//iteration counter
	int		itsMax = 120;					//Maximum iteration count
	int		subIts = 0;						//Sub-iterations in the TPX solver
	
	//Set initial guesses
	T		= *T_o;
	P		= *P_o;
	Q		= *Q_o;
	x		= *x_o;
	y		= *y_o;
	
	//Iteration
	while ( ( fabs(err) > eps ) && (i >= 0)  )
	{
		//Get current value of Q,x,y, U, and V
		subIts	= Qxy_TPX (&Q, &x, &y, T, P, X_i);
		if ( subIts < 0 )
		{	i = -1;		break; }
		
		Q	   = median( 0, Q, 1);
		x	   = (Q == 0) ? X_i : x;
		y	   = (Q == 1) ? X_i : y;
		if		( Q == 0 )
		{	U	= U_lm(T, P, x);	V	= V_lm(T, P, x); }
		else if	( Q == 1 )
		{	U	= U_vm(T, P, y);	V	= V_vm(T, P, y); }	
		else
		{	U	= Q*U_vm(T, P, y) + (1-Q)*U_lm(T, P, x);
			V	= Q*V_vm(T, P, y) + (1-Q)*V_lm(T, P, x); }
			
		//For T iterations
		if ( (i == 0) || ((i%2) == 0) )
		{
			//Perturbed values (by T)
			xT = x; yT = y; QT = Q;
			subIts = Qxy_TPX (&QT, &xT, &yT, T+deltaT, P, X_i);
			if ( subIts < 0 )
			{	i = -1;		break; }
			QT	   = median( 0, QT, 1);
			xT	   = (QT == 0) ? X_i : xT;
			yT	   = (QT == 1) ? X_i : yT;
			if		( QT == 0 )
			{	UT	= U_lm(T+deltaT, P, xT);	VT	= V_lm(T+deltaT, P, xT); }
			else if	( QT == 1 )
			{	UT	= U_vm(T+deltaT, P, yT);	VT	= V_vm(T+deltaT, P, yT); }	
			else
			{	UT	= QT*U_vm(T+deltaT, P, yT) + (1-QT)*U_lm(T+deltaT, P, xT);
				VT	= QT*V_vm(T+deltaT, P, yT) + (1-QT)*V_lm(T+deltaT, P, xT); }
		}
		if ( (i == 0) || ((i%2) == 1) )
		{
			//Perturbed values (by P)
			xP = x; yP = y; QP = Q;
			subIts = Qxy_TPX (&QP, &xP, &yP, T, P+deltaP, X_i);
			if ( subIts < 0 )
			{	i = -1;		break; }
			QP	   = median( 0, QP, 1);
			xP	   = (QP == 0) ? X_i : xP;
			yP	   = (QP == 1) ? X_i : yP;
			if		( QP == 0 )
			{	UP	= U_lm(T, P+deltaP, xP);	VP	= V_lm(T, P+deltaP, xP); }
			else if	( QP == 1 )
			{	UP	= U_vm(T, P+deltaP, yP);	VP	= V_vm(T, P+deltaP, yP); }	
			else
			{	UP	= QP*U_vm(T, P+deltaP, yP) + (1-QP)*U_lm(T, P+deltaP, xP);
				VP	= QP*V_vm(T, P+deltaP, yP) + (1-QP)*V_lm(T, P+deltaP, xP); }
		}
		
		//Evaluate Residual
		R[0]	= U_i - U;
		R[1]	= V_i - V;
		
		//Evaluate Jacobian
		J[0][0]	= -(UT - U)/deltaT;
		J[0][1]	= -(UP - U)/deltaP;
		J[1][0]	= -(VT - V)/deltaT;
		J[1][1]	= -(VP - V)/deltaP;
		
		//Calculate inverse
		DetInv		= 1/ ( J[0][0]*J[1][1] - J[0][1]*J[1][0] );
		Jinv[0][0]	=  DetInv*J[1][1];
		Jinv[0][1]	= -DetInv*J[0][1];
		Jinv[1][0]	= -DetInv*J[1][0];
		Jinv[1][1]	=  DetInv*J[0][0];
		
		//Calculate increments for T, P
		dT	= -( Jinv[0][0]*R[0] + Jinv[0][1]*R[1] );
		dP	= -( Jinv[1][0]*R[0] + Jinv[1][1]*R[1] );
		
		//Apply limiter
		//Under-relaxation, relative to magnitudes of T and P
		alphaT	= (fabs(dT) < Tlim) ? 1 : Tlim/fabs(dT);
		alphaP	= (fabs(dP) < 0.40*P) ? 1 : 0.40*P/fabs(dP);
		alpha	= min_val(alphaT, alphaP);
		
		
		//Update T and P
		if 		( (i%2) == 0 )
		{	T	+= alpha*dT; }
		else if	( (i%2) == 1 )
		{	P	+= alpha*dP; }
		//T	+= alpha*dT;
		//P	+= alpha*dP;
		
		//Update error
		err	= alpha*sqrt(dT*dT/(T*T) + dP*dP/(P*P) );
		
		//Update iteration counter
		i++;
		
		//Check iteration count limit
		if (i > itsMax)
			i = -1;

	}
	
	//We split iterations, so:
	if (i >= 0)
		i	/= 2;
	
	//Return final value of these properties
	Qxy_TPX (&Q, &x, &y, T, P, X_i);
	
	//Return results
	*T_o	= T;
	*P_o	= P;
	*Q_o	= Q;
	*x_o	= x;
	*y_o	= y;
	
	return i;
}
//Gets P from TQxy
void P_TQVxy ( double *P_o, double T_i, double Q_i, double V_i, double x_i, double y_i)
{
	//General variables
	double	P = *P_o;			//get input pressure
	double	V;					//current specific volume
	double	deltaP	= 1;		//increment for pressure
	
	//For iterative solution
	double	dP;					//Change in T
	double	Vp, Pp;				//previous volume and pressure
	double	eps = 1E-6;			//Convergence criterion
	double	err	= 1;			//Current error in H
	double	R;					//Current residual
	double	J;					//Jacobian
	int		i = 0;				//Iteration counter
	
	//iterative solve for T
	while (err > eps)
	{
		//Get current volume
		V	= Q_i*V_vm(T_i, P, y_i) + (1-Q_i)*V_lm(T_i, P, x_i);
		
		//Get residual
		R	= V_i	- V;
		
		//Get derivative
		if ( i == 0 )		//Initial run
		{	J	= -( Q_i*V_vm(T_i, P+deltaP, y_i) + (1-Q_i)*V_lm(T_i, P+deltaP, x_i) - V )/deltaP; }
		else
		{	J	= -( V - Vp )/( P - Pp ); }
		
		//Update old values of P, V
		Vp		= V;
		Pp		= P;
		
		//calculate change in T
		dP	= -R/J;
		
		//Apply limiter
		dP	= median( -0.5*P, dP, 0.5*P);
		
		//Update P
		P	+= dP;
			
		//update error
		err = fabs( dP/P );
		
		//update iteration counter
		i++;
	}
	
	//Return result
	*P_o	= P;
}

//Gets PQxy from TXV
int PQxy_TXV ( double *P_o, double *Q_o, double *x_o, double *y_o, double T_i, double X_i, double V_i )
{
	double	P, Pp;							//Concentration and updated concentration
	double	Q, Qp;							//Quality and previous quality
	double	V, Vp;							//Specific volume and previous specific volume
	double	x, y, xp, yp;					//Phasic concentrations
	
	double	dP;								//Concentration increment
	double	mu_lw, mu_la, mu_vw, mu_va;		//Chemical potentials
	double	eps = 1E-6;						//Convergence criterion
	double	err	= 1;						//Error at this iteration
	double	alpha = 1;						//under-relaxation factor
	double	deltaP	= 1;					//step in T for calculating derivatives
	double	errP = err;
	int		i = 0;							//iteration counter
	
	//new iteration scheme parameters
	double	dV_dP;
	double 	R, J, Jinv;
	
	double	Pm = 0;
	double	Vm = 0;							//Previous pressure, volume
	
	//Initialize P, Q, x and y
	P	= *P_o;
	Pp	= P;
	x	= *x_o;
	y	= *y_o;
	Q	= *Q_o;
	Q	= median( 0, Q, 1);
	Qp	= Q;
	
	
	//Get potentials for checking initial error
	mu_lw		= ChemPo_lw (T_i, P, x);
	mu_la		= ChemPo_la (T_i, P, x);
	mu_vw		= ChemPo_vw (T_i, P, y);
	mu_va		= ChemPo_va (T_i, P, y);
	//We use the relative error in the chemical potentials as the uncertainty
	err			= sqrt( pow( mu_lw-mu_vw, 2) + pow( mu_la-mu_va, 2) ) / sqrt( pow( 0.5*(mu_lw+mu_vw), 2) + pow( 0.5*(mu_la+mu_va), 2) );
	//Also get volume for checking initial error
	V			= Q*V_vm(T_i, P, y) + (1-Q)*V_lm(T_i, P, x);
	err			= max_val( err, fabs( (V-V_i)/V ) );
	errP 		= err;
	
	
	//New iteration scheme
	while ( fabs(err) > eps)
	{
		//First get residual
		Qxy_TPX (&Q, &x, &y, T_i, P, X_i);		//Get new Q,x,y
		if ( (Q < 0) || (Q > 1) )
		{
			if ( fabs(X_i-x) < fabs(X_i - y) )
			{	Q = 0; x = X_i; }
			else
			{	Q = 1; y = X_i; }
		}
		
		//x		= (Q <= 0) ? X_i : x;
		//y		= (Q >= 1) ? X_i : y;
		//Q		= median( 0, Q, 1);
		V		= Q*V_vm(T_i, P, y) + (1-Q)*V_lm(T_i, P, x);
		
		//Residual for volume
		R = V_i - V;
		
		//Now get Jacobian (-dV/dP)
		if (i == 0)	//No previous data to use to estimate dH/dT
		{
			Pp		= P+deltaP; 							//New T for estimating dH_dT
			xp	= x; yp = y; Qp = Q;
			Qxy_TPX (&Qp, &xp, &yp, T_i, Pp, X_i);			//Get new Q,x,y
			Qxy_TPX (&Qp, &xp, &yp, T_i, Pp, X_i);			//Get new Q,x,y
			//xp		= (Qp <= 0) ? X_i : xp;
			//yp		= (Qp >= 1) ? X_i : yp;
			//Qp		= median( 0, Qp, 1);
			
			if ( (Qp < 0) || (Qp > 1) )
			{
				if ( fabs(X_i-xp) < fabs(X_i - yp) )
				{	Qp = 0; xp = X_i; }
				else
				{	Qp = 1; yp = X_i; }
			}
			
			Vp		= Qp*V_vm(T_i, Pp, yp) + (1-Qp)*V_lm(T_i, Pp, xp);
			dV_dP  = (Vp - V)/deltaP;
		}
		//Sometimes derivative goes to 0 if changes in Q are too small, in this case, use previous value
		else if (V != Vm)	//Use previous H, T pair to estimate derivative
		{
			dV_dP  = (V - Vm)/(P - Pm);
		}
		//Evaluate Jacobian
		J		= -dV_dP;
		Jinv	= 1/J;
		
		//Set previous values for future estimations of dQ_dX
		Pm		= P;
		Vm		= V;
		
		//Get delta P
		dP		= -Jinv*R;
		
		//Apply limiter
		dP		= median( -0.2*P, dP, 0.2*P );
		
		//Update P
		P	   += alpha*dP;
		
		//Calculate error
		errP   = err;
		err    = dP/P;									//Update error
		
		//Avoid oscillations
		if ( sign(errP) != sign(err) )
		{	alpha	*= 0.8;	}
		else
		{	alpha	= min_val(1, 1.25*alpha); }
		
		//Break if nothing is changing
		if (err == errP)
		{	err = 0; }
		
		//iteration counter
		i++;
	}
	
	
	//Get final results for Q,x,y
	Qxy_TPX (&Q, &x, &y, T_i, P, X_i);
	
	//Return P, Q, x, and y
	*P_o	= P;
	*Q_o	= Q;
	*x_o	= x;
	*y_o	= y;
	
	//return iteration count
	return i;
}

//==================================================================================================
//GENERAL MIXTURE QUANTITIES
//==================================================================================================

//Molar mass of the mixture [kg/kmol]
double MW_m (double X)
{
	return ( X*MW_A + (1-X)*MW_W );
}

//==================================================================================================
//SOME HELPER FUNCTIONS
//==================================================================================================
double min_val(double a, double b)
{	return ( (a < b) ? a : b ); }
	
double max_val(double a, double b)
{	return ( (a > b) ? a : b ); }

double median(double a, double b, double c)
{
	if ( ( (a >= b) && (a <= c) ) || ( (a >= c) && (a <= b) ) )
	{	return a; }
	else if ( ( (b >= a) && (b <= c) ) || ( (b >= c) && (b <= a) ) )
	{	return b; }
	else
	{	return c; }
}

double sign(double a)
{
	return (a > 0) ? 1 : -1;
}

//==================================================================================================
//Matrix inverses
//==================================================================================================
void Inv2x2(double M[2][2], double Minv[2][2])
{
	double	DetInv;				//Inverse of determinant
	
	//Get inverse of determinant
	DetInv		= 1/ ( M[0][0]*M[1][1] - M[0][1]*M[1][0] );
	//Fill in individual entries
	Minv[0][0]	=  DetInv*M[1][1];
	Minv[0][1]	= -DetInv*M[0][1];
	Minv[1][0]	= -DetInv*M[1][0];
	Minv[1][1]	=  DetInv*M[0][0];
}

void Inv3x3(double M[3][3], double Minv[3][3])
{
	double	DetInv;				//Inverse of determinant
	
	//Get inverse of determinant
	DetInv		= 1/ ( M[2][0]*M[0][1]*M[1][2] - M[2][0]*M[0][2]*M[1][1] - M[1][0]*M[0][1]*M[2][2] + M[1][0]*M[0][2]*M[2][1] + M[0][0]*M[1][1]*M[2][2] - M[0][0]*M[1][2]*M[2][1] );
	//Fill in individual entries
	Minv[0][0]	=  DetInv*( M[1][1]*M[2][2] - M[2][1]*M[1][2] );
	Minv[0][1]	=  DetInv*( M[0][2]*M[2][1] - M[0][1]*M[2][2] );
	Minv[0][2]	=  DetInv*( M[0][1]*M[1][2] - M[0][2]*M[1][1] );
	Minv[1][0]	=  DetInv*( M[2][0]*M[1][2] - M[1][0]*M[2][2] );
	Minv[1][1]	=  DetInv*( M[0][0]*M[2][2] - M[2][0]*M[0][2] );
	Minv[1][2]	=  DetInv*( M[1][0]*M[0][2] - M[0][0]*M[1][2] );
	Minv[2][0]	=  DetInv*( M[1][0]*M[2][1] - M[2][0]*M[1][1] );
	Minv[2][1]	=  DetInv*( M[2][0]*M[0][1] - M[0][0]*M[2][1] );
	Minv[2][2]	=  DetInv*( M[0][0]*M[1][1] - M[1][0]*M[0][1] );
}
