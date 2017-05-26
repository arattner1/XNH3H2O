/*******************************************************
XNH3H2O.h									05/23/2011
Alexander Rattner							STSL
MEX library for getting thermophysical properties of
ammonia water mixtures
*******************************************************/

//Include libraries
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string.h>

//Include local code
#include "NH3H2O.h"
#include "kdtree.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //***********************INITIALIZATION SECTION*************************
    //*******VARIABLES ASSOCIATIED WITH MEX STUFF
    mxArray *mode_in_m;					//Function mode input
    mwSize buflen;						//buffer length for the mode
    mxArray *in1_m;						//Input 1
    mxArray *in2_m;						//Input 2
    mxArray *in3_m;						//Input 3
    mxArray *in4_m, *in5_m, *in6_m;     //Inputs 4-6, for when an initial guess is provided
    mxArray *out1_m;					//Output 1
    //Doubles of the above structures
    double *in1, *in2, *in3, *in4, *in5, *in6, *out1;
    //Mode data
    char *mode;
    //Prefix and suffix for the mode
    char *mode_pre, *mode_suf;
    //Filename for loading data
    char *fname;
    
    
    //*******BASIC VARIABLES
    Propset P;							//Holds property data
    Propset Pp;                         //Property value at previous time
    
    //*******STATIC VARIABLES (PERSIST BETWEEN RUNS OF XNH3H2O)
    //Database of property points
    static Proplib		PLib;
    //property tree structures
	static void			*kdTPX;
	static void			*kdPXH;
	static void			*kdTPQ;
	static void			*kdPXQ;
	static void			*kdXVU;
	static void			*kdTXV;
	static Trees		TS;
	static int			Initialized = 0;			//Set to 1 after initialization
    
    //Associate input values
    if	 	(nrhs == 7)	//We received a mode, state triple, and initial guess
    {
    	mode  = mxArrayToString(prhs[0]);
    	in1_m = mxDuplicateArray(prhs[1]);
    	in1   = mxGetPr(in1_m);
    	in2_m = mxDuplicateArray(prhs[2]);
    	in2   = mxGetPr(in2_m);
    	in3_m = mxDuplicateArray(prhs[3]);
    	in3   = mxGetPr(in3_m);
        in4_m = mxDuplicateArray(prhs[4]);
    	in4   = mxGetPr(in4_m);
        in5_m = mxDuplicateArray(prhs[5]);
    	in5   = mxGetPr(in5_m);
        in6_m = mxDuplicateArray(prhs[6]);
    	in6   = mxGetPr(in6_m);
    }
    if	 	(nrhs == 4)	//We received a mode and a state triple
    {
    	mode  = mxArrayToString(prhs[0]);
    	in1_m = mxDuplicateArray(prhs[1]);
    	in1   = mxGetPr(in1_m);
    	in2_m = mxDuplicateArray(prhs[2]);
    	in2   = mxGetPr(in2_m);
    	in3_m = mxDuplicateArray(prhs[3]);
    	in3   = mxGetPr(in3_m);
    }
    else if	(nrhs == 2)	//We received a mode and a string (presumably INIT and filename)
    {
    	mode  = mxArrayToString(prhs[0]);
    	fname = mxArrayToString(prhs[1]);
    }
    else if	(nrhs == 1)	//We only received a mode (presumably END and filename)
    {
    	mode  = mxArrayToString(prhs[0]);
    }
    
    //Prepare output structure, size depends on mode
    if      ( strcmp(mode, "DTPX_DXVU") == 0 )
    {   out1_m = mxCreateDoubleMatrix(1,9,mxREAL); }
    else if ( strcmp(mode, "XVU_TPX") == 0 )
    {   out1_m = mxCreateDoubleMatrix(1,3,mxREAL); }
    else
    {   out1_m = mxCreateDoubleMatrix(1,1,mxREAL); }
    
    plhs[0] = out1_m;
    out1  = mxGetPr(out1_m);
    
    //************INITIALIZATION COMPLETE, ACTUAL CODE BELOW*************
    if		( ( strcmp(mode, "INIT") == 0 ) && (nrhs == 2) )
    {
		if ( Initialized == 0 )
		{
			//Load libraries
			LoadPropData (fname, &PLib);
			BuildPropTree("TPX", &PLib, &kdTPX);
			TS.hasTPX = 1; TS.TPX = kdTPX;
			BuildPropTree("PXH", &PLib, &kdPXH);
			TS.hasPXH = 1; TS.PXH = kdPXH;
			BuildPropTree("TPQ", &PLib, &kdTPQ);
			TS.hasTPQ = 1; TS.TPQ = kdTPQ;
			BuildPropTree("PXQ", &PLib, &kdPXQ);
			TS.hasPXQ = 1; TS.PXQ = kdPXQ;
			BuildPropTree("XVU", &PLib, &kdXVU);
			TS.hasXVU = 1; TS.XVU = kdXVU;
			BuildPropTree("TXV", &PLib, &kdTXV);
			TS.hasTXV = 1; TS.TXV = kdTXV;
            
			//Set initialized
			Initialized = 1;
			*out1 = 1;
		}
		else
		{
			mexPrintf("XNH3H2O already initialized\n");
			*out1 = -1;
		}
    }
    else if	( ( strcmp(mode, "SAVE") == 0 ) && (nrhs == 2) )
    {
		if ( Initialized == 1 )
		{
			//Save Library
			SavePropData (fname, &PLib);
			*out1 = 1;
		}
		else
		{
			mexPrintf("XNH3H2O not yet initialized\n");
			*out1 = -1;
		}
    }
    else if ( strcmp(mode, "END") == 0 )
    {
    	if ( Initialized == 1 )
		{
    		//Clean up
			mxFree(PLib.D);
			kd_free(kdTPX);
			kd_free(kdPXH);
			kd_free(kdTPQ);
			kd_free(kdPXQ);
			kd_free(kdXVU);
			kd_free(kdTXV);
	    	Initialized = 0;
    		*out1 = 1;
    	}
    	else
    	{
    		mexPrintf("XNH3H2O not yet initialized\n");
    		*out1 = -1;
    	}
    }
    else if ( strcmp(mode, "DTPX_DXVU") == 0 )
    {
        if ( Initialized == 1 )
		{
    		//Gets a vector with the contents:
            //dT/dX, dT/dv, dT/du, dP/dX, dP/dv, dP/du, dX/dX, dX/dv, dX/du
            //Initialize propset
            P.T	= *in1;		P.P = *in2;		P.psi = *in3;
            //Get derivatives
            getDTPX_DXVU(&P, &TS, &PLib, out1);
    	}
    	else
    	{
    		mexPrintf("XNH3H2O not yet initialized\n");
    		out1[0] = -1;
    	} 
    }
    //Input contains, X_new, V_new, U_new, T_old, P_old, X_old
    else if ( ( strcmp(mode, "XVU_TPX") == 0 ) && (nrhs == 7) )
    {
        if ( Initialized == 1 )
		{
            //Set XVU in P
            P.psi = *in1;   P.v = *in2;     P.u = *in3;
            //Set TPX in Pp
            Pp.T = *in4;    Pp.P = *in5;    Pp.psi = *in6;
    		//Update TPX from XVU
            getTPX_XVU(&P, &P, &Pp,  &TS, &PLib);
            //Set output values
            out1[0] = P.T; out1[1] = P.P; out1[2] = P.psi;
    	}
    	else
    	{
    		mexPrintf("XNH3H2O not yet initialized\n");
    		out1[0] = -1;
    	} 
    }
    //Triples input, already initialized
    else if ( (nrhs == 4) && (Initialized == 1) )
    {
    	//First break out prefix and suffix
    	mode_pre	= strtok(mode, "_");
    	mode_suf	= strtok(NULL, "_");
    	
    	//Now set the right function inputs
    	if		( strcmp(mode_pre, "TPX") == 0 )
    	{	P.T	= *in1;		P.P = *in2;		P.psi = *in3; }
    	else if	( strcmp(mode_pre, "TPQ") == 0 )
    	{	P.T	= *in1;		P.P = *in2;		P.Q = *in3; }
    	else if	( strcmp(mode_pre, "PXH") == 0 )
    	{	P.P	= *in1;		P.psi = *in2;	P.h = *in3; }
    	else if	( strcmp(mode_pre, "PXQ") == 0 )
    	{	P.P	= *in1;		P.psi = *in2;	P.Q = *in3; }
    	else if	( strcmp(mode_pre, "XVU") == 0 )
    	{	P.psi = *in1;	P.v = *in2;		P.u = *in3; }
    	else if	( strcmp(mode_pre, "TXV") == 0 )
    	{	P.T = *in1;		P.psi = *in2;	P.v = *in3; }
    	
    	//Make the property call
    	SetPropsSearch(mode_pre, &P, &TS, &PLib);
    	
    	//Return the right datum
    	if		( strcmp(mode_suf, "T") == 0 )	{	*out1 = P.T; }
    	else if	( strcmp(mode_suf, "P") == 0 )	{	*out1 = P.P; }
    	else if	( strcmp(mode_suf, "X") == 0 )	{	*out1 = P.psi; }
    	else if	( strcmp(mode_suf, "Q") == 0 )	{	*out1 = P.Q; }
    	else if	( strcmp(mode_suf, "H") == 0 )	{	*out1 = P.h; }
    	else if	( strcmp(mode_suf, "U") == 0 )	{	*out1 = P.u; }
    	else if	( strcmp(mode_suf, "V") == 0 )	{	*out1 = P.v; }
    	else if	( strcmp(mode_suf, "S") == 0 )	{	*out1 = P.s; }
    	else if	( strcmp(mode_suf, "x") == 0 )	{	*out1 = P.x; }
    	else if	( strcmp(mode_suf, "y") == 0 )	{	*out1 = P.y; }
    }
    //Triples input, not yet initialized
    else if ( (nrhs == 4) && (Initialized == 0) )
    {
    	mexPrintf("XNH3H2O not yet initialized\n");
    	*out1 = -1;
    }
    return;
}
