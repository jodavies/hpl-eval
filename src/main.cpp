/*
 * Arbitrary precision evaluation of expressions written in terms of
 * Harmonic Polylogarithms, natural logarithms, simple polynomials and
 * denominators of x and Riemann zeta functions.
 * 
 * Uses the GiNaC framework for evaluation.
 * 
 * Expressions produced and formatted appropriately by FORM.
 * 
 * Mellin moment evaluation routine expects expression to be split into
 * regular and singular (the plus-distribution) pieces.
 * 
 * Convolution (not yet implemented, will be unusably slow?) requires
 * regular, singular and local pieces to be passed to the function.
 */

#include <iostream>
#include <iomanip>
using namespace std;
#include <ginac/ginac.h>
using namespace GiNaC;


// Include expressions.
#include "EC3NP3.h"
#include "ECLNP3.h"
#include "EC2NP3.h"
#include "smallx.h"



// Function pointer for passing expressions
typedef ex (*ExpressionPtr)(symbol x, numeric nf);

// Need a null expression of the same signature as ExpressionPtr
ex NullExpr(symbol x, numeric nf) {
	return 0;
}


// Evaluate at a logarithmically ~evenly spaced range of x values, for
// plotting and fitting. Must add together regular and singular part of
// the expression.
// Print results for nfuncs expressions, passed in as arrays of ExpressionPtr.
void EvaluateAtX(int nfuncs, ExpressionPtr * regularExpressions,
                 ExpressionPtr * singularExpressions, numeric nf,
                 double scaleFactor);


// Compute specified mellin moment. This is VERY SLOW compared to the
// fortran code it is intended to replace. This is just a check anyway,
// it is far better to evaluate moments in Mellin space...
void EvaluateMellinMoment(ExpressionPtr regular, ExpressionPtr singular,
                          ExpressionPtr local, numeric nval, numeric nf);


int main(void)
{
	Digits = 15;
	symbol x("x");

	// number of flavours
	numeric nf(4);

	ExpressionPtr regularExpressions[] = {EC2NP3A,
	                                      c2ns3OSXLL,
	                                      c2ns3OSXNLL,
	                                      c2ns3OSXNNLL,
	                                      c2ns3OSXNNNLL,
	                                      c2ns3OSXNNNNLL,
	                                      ECLNP3A,
	                                      cLns3OSXLL,
	                                      cLns3OSXNLL,
	                                      cLns3OSXNNLL,
	                                      EC3NP3A,
	                                      c3ns3ESXLL,
	                                      c3ns3ESXNLL,
	                                      c3ns3ESXNNLL,
	                                      c3ns3ESXNNNLL,
	                                      c3ns3ESXNNNNLL,
	                                     };

	ExpressionPtr singularExpressions[] ={EC2NP3B,
	                                      NullExpr,
	                                      NullExpr,
	                                      NullExpr,
	                                      NullExpr,
	                                      NullExpr,
	                                      ECLNP3B,
	                                      NullExpr,
	                                      NullExpr,
	                                      NullExpr,
	                                      EC3NP3B,
	                                      NullExpr,
	                                      NullExpr,
	                                      NullExpr,
	                                      NullExpr,
	                                      NullExpr,
	                                     };

	EvaluateAtX(16, regularExpressions, singularExpressions, nf, 1.0/2000.0);




	// Investigate accuracy of evaluations near x=1, for integration issues.

/*	ex integrand = (pow(x, nval-numeric(1))*EC3NP3A(x,nf));
//	ex integrand = (pow(x, nval-numeric(1)-numeric(1))*EC3NP3B(x,nf));

	cout << integrand << endl << endl;
//	cout << integrand.evalf() << endl << endl;



	cout << numeric(1.0E-12) << " " << integrand.subs(x==numeric(1.0E-12)).evalf() << endl;
	cout << numeric(0.5) << " " << integrand.subs(x==numeric(0.5)).evalf() << endl;
//	cout << numeric(0.99999999) << " " << integrand.subs(x==numeric(0.99999999)).evalf() << endl;
	cout << numeric(0.99999999) << " " << integrand.evalf().subs(x==numeric(0.99999999)).evalf() << endl;
	cout << numeric(0.25) << " " << integrand.subs(x==numeric(0.25)).evalf() << endl;
	cout << numeric(0.75) << " " << integrand.subs(x==numeric(0.75)).evalf() << endl;
	cout << numeric(0.125) << " " << integrand.subs(x==numeric(0.125)).evalf() << endl;
	cout << numeric(0.375) << " " << integrand.subs(x==numeric(0.375)).evalf() << endl << endl;

	//cout << integral(x, numeric(1,1000000000000), numeric(99999999,100000000),  integrand).evalf() << endl;
	
//	adaptivesimpson(x, numeric(1,1000000000000), numeric(99999999,100000000),  integrand);*/


/*	EvaluateMellinMoment(&EC3NP3A, &EC3NP3B, &EC3NP3C, numeric(2), nf);
	EvaluateMellinMoment(&EC3NP3A, &EC3NP3B, &EC3NP3C, numeric(4), nf);
	EvaluateMellinMoment(&EC3NP3A, &EC3NP3B, &EC3NP3C, numeric(6), nf);
	EvaluateMellinMoment(&EC3NP3A, &EC3NP3B, &EC3NP3C, numeric(8), nf);
	EvaluateMellinMoment(&EC3NP3A, &EC3NP3B, &EC3NP3C, numeric(10), nf);
	EvaluateMellinMoment(&EC3NP3A, &EC3NP3B, &EC3NP3C, numeric(12), nf);*/

	return 0;
}



// Evaluate at a logarithmically ~evenly spaced range of x values, for
// plotting and fitting. Must add together regular and singular part of
// the expression.
void EvaluateAtX(int nfuncs, ExpressionPtr * regularExpressions,
                 ExpressionPtr * singularExpressions, numeric nf, double scaleFactor)
{
	// large-x x values, not computed in loop.
	numeric XL[10] = {numeric(995,1000),
	                  numeric(999,1000),
	                  numeric(9995,10000),
	                  numeric(9999,10000),
	                  numeric(99995,100000),
	                  numeric(99999,100000),
	                  numeric(999995,1000000),
	                  numeric(999999,1000000),
	                  numeric(9999995,10000000),
	                  numeric(9999999,10000000)};
	symbol x("x");
	
	// Loop over x values and evaluate. Each 50 K1 values below 150 extends
	// range down by a power of 10. -200: 10^-8,
	//                              -150: 10^-7,
	//                              -100: 10^-6, ...
	for (int K1 = -200; K1 <= 249; K1++) {
		numeric xval,ax;
		if (K1 <= 150) {
			ax = - numeric(92103404,10000000)
			     + numeric(K1*1000000,21714724);
			xval = exp(ax);
		}
		else if (K1 <= 239) {
			xval = numeric(K1-140, 100);
		}
		else {
			xval = XL[K1-240];
		}

		cout << evalf(xval);
		for (int i = 0; i < nfuncs; i++) {
			cout << " "
			     << scaleFactor*real_part(evalf(+regularExpressions[i](x,nf).subs(x == xval)
			                        +singularExpressions[i](x,nf).subs(x == xval) ));
		}
		cout << endl;
	}

}



// Compute specified mellin moment. This is VERY SLOW compared to the
// fortran code it is intended to replace. This is just a check anyway,
// it is far better to evaluate moments in Mellin space...
void EvaluateMellinMoment(ExpressionPtr regular, ExpressionPtr singular,
                          ExpressionPtr local, numeric nval, numeric nf)
{
	symbol x("x");

	numeric low(1,1000000000000);
	numeric high(999999999,1000000000);

	ex integrand = (+ pow(x, nval-numeric(1)) * regular(x,nf)
	                +(pow(x, nval-numeric(1))-numeric(1)) * singular(x,nf)
	               );

	integral::max_integration_level = 150; //default == 15
	integral::relative_integration_error = 1.0E-8; //default = 1.0E-8

//	THIS IS NOT ACCURATE
//	ex mellinIntegral = integral(x, low, high, integrand);
//	cout << real_part(mellinIntegral.evalf() + evalf(local(x,nf).subs(x == low))) << endl;

	// Need to call adaptivesimpson directly, to avoid an evalf before
	// symbol substitution, which causes numerics to blow up.
	ex mellinIntegral = adaptivesimpson(x, low, high, integrand);
	cout << real_part(mellinIntegral) + evalf(local(x,nf).subs(x==low)) << endl;
}
