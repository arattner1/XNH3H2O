/*******************************************************
NH3H2O.h									02/28/2011
Alexander Rattner							STSL
Contains constants, structures, and function prototypes
for the thermophysical properties of Ammonia-Water fluid 
mixtures.
*******************************************************/


#ifndef NH3H2O_H
#define NH3H2O_H

//include external libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//include local libraries
#include "kdtree.h"

//Include MATLAB libraries
#include <matrix.h>
#include <mex.h>

//==================================================================================================
//Universal constants
#define T_B		100			//Reference temperature in [K]
#define P_B		1.0E6		//Reference pressure in [Pa]
#define R_U		8314.47		//Universal gas constant [J/kmol-K]
#define MW_W	18.01530	//Molecular weight of water in [kg/kmol]
#define MW_A	17.03056	//Molecular weight of ammonia in [kg/kmol]
//Reference states of water
#define PC_W	22.064E6	//Critical pressure of water in [Pa]
#define PT_W	611.8321	//Triple-point pressure of water in [Pa]
#define TC_W	647.1		//Critical temperature of water in [K]
#define TT_W	273.16		//Triple-point Temperature of water in [K]
#define TB_W	373.15		//Normal boiling temperature of water in [K]
//Reference states of ammonia
#define PC_A	11.339E6	//Critical pressure of ammonia in [Pa]
#define PT_A	6090.9		//Triple-point pressure of ammonia in [Pa]
#define TC_A	405.4		//Critical temperature of ammonia in [K]
#define TT_A	195.5001	//Triple-point Temperature of ammonia in [K]
#define TB_A	239.8225	//Normal boiling temperature of ammonia in [K]
//==================================================================================================

//==================================================================================================
//Constants for Ammonia Water property evaluations
//orignally reported in: Ziegler, B., Trepp, C. Equation of state for Ammonia-Water. International 
//Journal of Refrigeration, Volume 7 (2), March 1984, P 101-106.
//FOR WATER
#define	H_R_L_O_W	 2.1821141E1		//Reference reduced enthalpy for liquid water
#define	H_R_V_O_W	 6.0965058E1		//Reference reduced enthalpy for vapor water
#define	S_R_L_O_W	 5.733498			//Reference reduced entropy for liquid water
#define	S_R_V_O_W	 1.3453430E1		//Reference reduced entropy for vapor water
#define	T_R_O_W		 5.0705				//Reference reduced temperature for water
#define	P_R_O_W		 3.0000				//Reference reduced pressure for water
//Constants for Gibbs equations for water
#define	A_1_W		 2.748796E-2
#define	A_2_W		-1.016665E-5
#define	A_3_W		-4.452025E-3
#define	A_4_W		 8.389246E-4
#define	B_1_W		 1.214557E1
#define	B_2_W		-1.898065
#define	B_3_W		 2.911966E-1
#define	C_1_W		 2.136131E-2
#define	C_2_W		-3.169291E1
#define	C_3_W		-4.634611E4
#define C_4_W		 0.0
#define	D_1_W		 4.019170
#define	D_2_W		-5.175550E-2
#define	D_3_W		 1.951939E-2
//Constants for finding saturation curve (Velasco et al., 2007)
#define F_1_W		 0.17027
#define F_2_W		-0.21182
//FOR AMMONIA
#define	H_R_L_O_A	 4.878573			//Reference reduced enthalpy for liquid ammonia
#define	H_R_V_O_A	 2.6468873E1		//Reference reduced enthalpy for vapor ammonia
#define	S_R_L_O_A	 1.644773			//Reference reduced entropy for liquid ammonia
#define	S_R_V_O_A	 8.339026		//Reference reduced entropy for vapor ammonia
#define	T_R_O_A		 3.2252				//Reference reduced temperature for ammonia
#define	P_R_O_A		 2.0000				//Reference reduced pressure for ammonia
//Constants for Gibbs equations for ammonia
#define	A_1_A		 3.971423E-2
#define	A_2_A		-1.790557E-5
#define	A_3_A		-1.308905E-2
#define	A_4_A		 3.752836E-3
#define	B_1_A		 1.634519E1
#define	B_2_A		-6.508119
#define	B_3_A		 1.448937
#define	C_1_A		-1.049377E-2
#define	C_2_A		-8.288224
#define	C_3_A		-6.647257E2
#define C_4_A		-3.045352E3
#define	D_1_A		 3.673647
#define	D_2_A		 9.989629E-2
#define	D_3_A		 3.617622E-2
//Constants for finding saturation curve (Velasco et al., 2007)
#define F_1_A		 0.12442
#define F_2_A		-0.17763
//Constants for excess Gibbs energy of mixing
#define	E_1			-4.1733398E1
#define	E_2			 2.414E-2
#define	E_3			 6.702285 //!!!
#define	E_4			-1.1475E-2
#define	E_5			 6.3608967E1
#define	E_6			-6.2490768E1
#define	E_7			 1.761064
#define	E_8			 8.626E-3
#define	E_9			 3.87983E-1
#define	E_10		-4.772E-3
#define	E_11		-4.648107
#define	E_12		 8.36376E-1
#define	E_13		-3.553627
#define	E_14		 9.04E-4
#define	E_15		 2.4361723E1
#define	E_16		-2.0736547E1
//==================================================================================================

//==================================================================================================
//Reference values, taken at (T = 273.16C and P = 611.7 Pa)
//For liquid water:
#define H_REF_L_W	 -0.133
#define G_REF_L_W	 -0.2101
#define U_REF_L_W	 -0.1403
#define S_REF_L_W	  0
//For vapor water:
#define G_REF_V_W   334.2034
#define H_REF_V_W  1208.3
#define S_REF_V_W     3.2004
#define U_REF_V_W  1202.5
//For liquid ammonia:
#define H_REF_L_A    -0.002
#define G_REF_L_A    -1.1846
#define U_REF_L_A     0.0017
#define S_REF_L_A     0.0043
//For vapor water:
#define G_REF_V_A -5419.3
#define H_REF_V_A -5347.9
#define S_REF_V_A    -5.5428
#define U_REF_V_A -5547.1

//==================================================================================================

//forward definitions of structs
typedef struct Propset Propset;
typedef struct Proplib Proplib;
typedef struct Trees Trees;

//Contains the set of properties for a mixture
struct Propset
{
	//Thermodynamic properties
	double		T;			//Temperature in [K]
	double		P;			//Pressure in [Pa]
	double		psi;		//Mass fraction of ammonia [-]
	double		x, y;		//Mass fraction of ammonia in liquid and vapor phases, respectively
	double		Q;			//Quality of the mixture [-]
	double		MW;			//Molar mass in [kg/kmol]
	double		rho;		//Density in [kg/m^3]
	double		v;			//Specifc volume in [m^3/kg]
	double		cp;			//Specific heat in [J/kg-K]
	double		h;			//Enthalpy in [J/kg]
	double		u;			//Internal energy in [J/kg]
	double		s;			//Entropy in [J/kg-K]
	
	//Transport properties
	double		mu;			//Dynamic viscosity in [kg/m-s]
	double		nu;			//Kinematic viscosity in [m^2/s]
	double		k;			//Thermal conductivity in [W/m-K]
	double		alpha;		//Thermal diffusivity in [m^2/s]
	double		D_ab;		//Binary diffusion coefficient in [m^2/s]
	double		sigma;		//Surface tension in [N/m]
	
};

//Contains a library of property points
struct Proplib
{
	Propset		*D;				//Array of property data
	int			n;				//Current number of data points
	int			N;				//Maximum number of data points that have been malloced
	int			NBlock;			//By default, malloc for additional data points in units of NBlock
};

//structure of kdtrees
struct Trees
{
	int				hasTPX;			//If this struct contains a TPX tree
	void			*TPX;			//TPX tree
	int				hasPXH;			//If this struct contains a PXH tree
	void			*PXH;			//PXH tree
	int				hasTPQ;			//If this struct contains a TPQ tree
	void			*TPQ;			//TPQ tree
	int				hasPXQ;			//If this struct contains a PXQ tree
	void			*PXQ;			//PXQ tree
	int				hasXVU;			//If this struct contains a XVU tree
	void			*XVU;			//XVU tree
	int				hasTXV;			//If this struct contains a TXV tree
	void			*TXV;			//TXV tree
};

//function prototypes
//Top level functions
void SetProps(char *Mode, Propset *P);
void SetPropsSearch(char *Mode, Propset *P, Trees *TS, Proplib *PL);
void PrintProps(Propset *P);
void PropTreeSearch(double a, double b, double c, Propset *P, void *Tree);
void LoadPropData (char *fName, Proplib *L);
void SavePropData (char *fname, Proplib *L);
//Individual query modes
int getTPX(Propset *P, void *Tree);
int getTPXLite(Propset *P, void *Tree); //To be called from XVU routine
int getPXH(Propset *P, void *Tree);
int getTPQ(Propset *P, void *Tree);
int getPXQ(Propset *P, void *Tree);
int getXVU(Propset *P, void *Tree);
int getXVU2(Propset *Pin,  Trees *TS, Proplib *PL);
int getTXV(Propset *P, void *Tree);
//Special modes
void getDTPX_DXVU(Propset *P, Trees *TS, Proplib *PL, double *Jvec);
void getTPX_XVU(Propset *PTPX_n, Propset *PXVU_n, Propset *PTPX_p,  Trees *TS, Proplib *PL);
//Manipulate property database
void BuildPropTree(char *Mode, Proplib *L, void **Tree);
void AddPropPoint(char *Mode, Propset *P, Proplib *L, Trees *TS);
//General water properties
double Psat_w (double T);
double Tsat_w (double P);
//Liquid water properties
double G_lw (double T, double P);
double H_lw (double T, double P);
double V_lw (double T, double P);
double U_lw (double T, double P);
double Rho_lw (double T, double P);
double S_lw (double T, double P);
double Cp_lw (double T, double P);
double ChemPo_lw (double T, double P, double x);
//Vapor water properties
double G_vw 	(double T, double P);
double H_vw 	(double T, double P);
double V_vw 	(double T, double P);
double U_vw 	(double T, double P);
double Rho_vw 	(double T, double P);
double S_vw 	(double T, double P);
double Cp_vw 	(double T, double P);
double ChemPo_vw (double T, double P, double y);
//General water properties
double Psat_a (double T);
double Tsat_a (double P);
//Liquid ammonia properties
double G_la 	(double T, double P);
double H_la 	(double T, double P);
double V_la 	(double T, double P);
double U_la 	(double T, double P);
double Rho_la 	(double T, double P);
double S_la 	(double T, double P);
double Cp_la 	(double T, double P);
double ChemPo_la (double T, double P, double x);
//Vapor ammonia properties
double G_va 	(double T, double P);
double H_va 	(double T, double P);
double V_va 	(double T, double P);
double U_va 	(double T, double P);
double Rho_va 	(double T, double P);
double S_va 	(double T, double P);
double Cp_va 	(double T, double P);
double ChemPo_va (double T, double P, double y);
//Vapor mixtures
double H_vm 	(double T, double P, double y);
double V_vm 	(double T, double P, double y);
double Rho_vm 	(double T, double P, double y);
double S_vm 	(double T, double P, double y);
double U_vm 	(double T, double P, double y);
double G_vm 	(double T, double P, double y);
double Cp_vm 	(double T, double P, double y);
double Tdew		(double P, double y);
double ydew		(double T, double P);
double Pdew		(double y, double v);
double Pdew_TY	(double T_i, double Y_i);
//Liquid mixtures
double G_lm 	(double T, double P, double x);
double H_lm 	(double T, double P, double x);
double V_lm 	(double T, double P, double x);
double Rho_lm 	(double T, double P, double x);
double S_lm 	(double T, double P, double x);
double U_lm 	(double T, double P, double x);
double Cp_lm 	(double T, double P, double x);
double Tbub		(double P, double x);
double xbub		(double T, double P);
double Pbub		(double x, double v);
double Pbub_TX	(double T_i, double X_i);
//General mixtures
double MW_m		(double X);

//Derivatives for applying Newton-Rhapson for iterative solutions
double d_mu_lw_d_x (double T, double P, double x);
double d_mu_la_d_x (double T, double P, double x);
double d_mu_vw_d_y (double T, double P, double y);
double d_mu_va_d_y (double T, double P, double y);

//Full iterative solution (return iteration count)
int Qxy_TPX (double *Q_o, double *x_o, double *y_o, double T_i, double P_i, double X_i);	//Gets Quality, x, and y from TPX
int TQxy_PXH (double *T_o, double *Q_o, double *x_o, double *y_o, double P_i, double X_i, double H_i, double Tlim);	//Gets Temperature, Quality, x and y from PXH
int Txy_PXQ (double *T_o, double *x_o, double *y_o, double P_i, double X_i, double Q_i, double Tlim);
void T_PQHxy (double *T_o, double P_i, double Q_i, double H_i, double x_i, double y_i);		//Gets Temperature from PQHxy
void T_PHx_l (double *T_o, double P_i, double H_i, double x_i);		//Find temperature for the liquid mixture
void T_PHy_v (double *T_o, double P_i, double H_i, double y_i);		//Find temperature for the vapor mixture
void TP_QVUxy (double *T_o, double *P_o, double Q_i, double V_i, double U_i, double x_i, double y_i);		//Gets T,P from QVUxy
int TPQxy_XVU (double *T_o, double *P_o, double *Q_o, double *x_o, double *y_o, double X_i, double V_i, double U_i, double Tlim, double Plim);		//Gets T,P from QVUxy
void P_TQVxy ( double *P_o, double T_i, double Q_i, double V_i, double x_i, double y_i);
int PQxy_TXV ( double *P_o, double *Q_o, double *x_o, double *y_o, double T_i, double X_i, double V_i );
//A couple helper functions
double min_val(double a, double b);
double max_val(double a, double b);
double median(double a, double b, double c);
double sign(double a);
//Matrix inverses
void Inv2x2(double M[2][2], double Minv[2][2]);
void Inv3x3(double M[3][3], double Minv[3][3]);

#endif
