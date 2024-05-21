#ifndef POSTPROCESS_INLINE_H
#define POSTPROCESS_INLINE_H

#include <math.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <vector>
#include "Common-inline.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using Eigen::Matrix3d;
using Eigen::Vector3d;

using namespace std::complex_literals;

typedef std::vector<std::complex<double>> COMP_VEC;
typedef std::tuple< std::vector<double>, std::vector<double>, std::vector<double> > VSH_MODES_REAL;
typedef std::tuple< std::vector<double>, std::vector<double>, std::vector<double>, COMP_VEC, COMP_VEC, COMP_VEC > VSH_MODES_REAL_COMPLEX;


inline double LmMdivbyLpM(int l, int m) // (l-m)!/(l+m)!
{
	double result=1.0;
	if(m==0) return result;
	else if(m>0) for(int x=0; x<2*m; x++) result /= (double)(l+m-x);
	else for(int x=0; x<2*m; x++) result *= (double)(l+m-x);
	return result;
}

inline double dPlm_dTheta(int l, int m, double cosTheta)
{
	double Plmp1 = boost::math::legendre_p(l,m+1,cosTheta);
	double Plmm1 = boost::math::legendre_p(l,m-1,cosTheta);
	double result = 0.5*Plmp1 - 0.5*(l+m)*(l-m+1.0)*Plmm1;
	return result;
}

inline double realYlm(int l,int m,double cosTheta, double phi)
{
  double result =0.;
  if(m==0 || l==0) result=std::sqrt((2.0*l+1.0)/(4.0*M_PI)) * boost::math::legendre_p(l,m,cosTheta);

  if(m>0){
   // double lMinusMFact = boost::math::factorial<double>(l-m);
   // double lPlusMFact = boost::math::factorial<double>(l+m);
    double LmMdivbyLpMFac = LmMdivbyLpM(l,m);
    //result = sqrt(2.0) * sqrt((2.0*l+1.0)/(4.0*PI) * lMinusMFact/lPlusMFact) * boost::math::legendre_p(l,m,cosTheta) * cos(m*phi);
    result = std::sqrt(2.0) * std::sqrt((2.0*l+1.0)/(4.0*M_PI) * LmMdivbyLpMFac) * boost::math::legendre_p(l,m,cosTheta) * std::cos(m*phi);
  }else if(m<0){
    double absm = std::abs(m);
    //double lMinusMFact = boost::math::factorial<double>(l-absm);
    //double lPlusMFact = boost::math::factorial<double>(l+absm);
    double LmMdivbyLpMFac = LmMdivbyLpM(l,absm);
    result = std::sqrt(2.0) * std::sqrt(2.0*l/(4.0*M_PI) * LmMdivbyLpMFac) * boost::math::legendre_p(l,absm,cosTheta) * std::sin(absm*phi);
  }
  return result;
}

inline std::complex<double> real_to_complex_prefactor(int m, double phi)
{
	std::complex<double> result(1,0);
	if(m==0) return result;
	else if(m<0)
	{
		double sm = std::sin(std::abs(m)*phi);
		double cm = std::cos(std::abs(m)*phi);
		result = std::sqrt(2.)*sm/(cm-1i*sm);
	}
	else
	{
		double sm = std::sin(m*phi);
		double cm = std::cos(m*phi);
		result = std::pow(-1.,m)*std::sqrt(2.)*cm/(cm+1i*sm);
	}
	return result;
}


inline double dYlm_dTheta(int l, int m, double cosTheta, double phi)
{
  double result =0;
  if(l==0) return result;

  if(m==0) result= std::sqrt((2.0*l+1.0)/(4.0*M_PI)) * dPlm_dTheta(l, m, cosTheta);

  if(m>0){
    //double lMinusMFact = boost::math::factorial<double>(l-m);
    //double lPlusMFact = boost::math::factorial<double>(l+m);
    double LmMdivbyLpMFac = LmMdivbyLpM(l, m);
    result = std::sqrt(2.0) * std::sqrt((2.0*l+1.0)/(4.0*M_PI) * LmMdivbyLpMFac) * dPlm_dTheta(l, m, cosTheta) * std::cos(m*phi);
  }else if(m<0){
    double absm = std::abs(m);
    //double lMinusMFact = boost::math::factorial<double>(l-absm);
    //double lPlusMFact = boost::math::factorial<double>(l+absm);
    double LmMdivbyLpMFac = LmMdivbyLpM(l, absm);
    result = std::sqrt(2.0) * std::sqrt((2.0*l+1.0)/(4.0*M_PI) * LmMdivbyLpMFac) * dPlm_dTheta(l, absm, cosTheta) * std::sin(absm*phi);
  }
  return result;
}

inline double dYlm_dPhi(int l, int m, double cosTheta, double phi)
{
  if(m==0 || l==0) return 0.0;

  double result =0;
  if(m>0){
    //double lMinusMFact = boost::math::factorial<double>(l-m);
    //double lPlusMFact = boost::math::factorial<double>(l+m);
    double LmMdivbyLpMFac = LmMdivbyLpM(l,m);
    result = std::sqrt(2.0) * std::sqrt((2.0*l+1.0)/(4.0*M_PI) * LmMdivbyLpMFac) * boost::math::legendre_p(l,m,cosTheta) * (-1.0*m) * std::sin(m*phi);
  }else if(m<0){
    double absm = std::abs(m);
    //double lMinusMFact = boost::math::factorial<double>(l-absm);
    //double lPlusMFact = boost::math::factorial<double>(l+absm);
    double LmMdivbyLpMFac = LmMdivbyLpM(l,absm);
    result = std::sqrt(2.0) * std::sqrt((2.0*l+1.0)/(4.0*M_PI) * LmMdivbyLpMFac) * boost::math::legendre_p(l,absm,cosTheta) * (1.0*absm) * std::cos(absm*phi);
  }
  return result;
}

inline double dYlm_dPhiDivBySin(int l, int m, double cosTheta, double phi)
{
  double result =0;

  if(m==0 || l==0) return result;
  else if(m>0){
    //double lMinusMFact = boost::math::factorial<double>(l-m);
    //double lPlusMFact = boost::math::factorial<double>(l+m);
    double LmMdivbyLpMFac = LmMdivbyLpM(l,m);
    double fp = 0.5 * ( boost::math::legendre_p(l+1,m+1,cosTheta) + (l-m+1.0)*(l-m+2.0)*boost::math::legendre_p(l+1,m-1,cosTheta) );
    result = std::sqrt(2.0) * std::sqrt((2.0*l+1.0)/(4.0*M_PI) * LmMdivbyLpMFac) * fp * std::sin(m*phi);
  }else if(m<0){
    double absm = std::abs(m);
    //double lMinusMFact = boost::math::factorial<double>(l-absm);
    //double lPlusMFact = boost::math::factorial<double>(l+absm);
    double LmMdivbyLpMFac = LmMdivbyLpM(l,absm);
    double fp = -0.5 * ( boost::math::legendre_p(l+1,absm+1,cosTheta) + (l-absm+1.0)*(l-absm+2.0)*boost::math::legendre_p(l+1,absm-1,cosTheta) );
    result = std::sqrt(2.0) * std::sqrt((2.0*l+1.0)/(4.0*M_PI) * LmMdivbyLpMFac) * fp * std::cos(absm*phi);
  }
  return result;
}

inline std::pair<Vector3d,Vector3d> PSIlmPHIlm(int l, int m, double cosTheta, double phi, Vector3d theta_hat, Vector3d phi_hat)
{
  double t = dYlm_dTheta(l, m, cosTheta, phi);
  double p = dYlm_dPhiDivBySin(l, m, cosTheta, phi);
  Vector3d PSIlm = t*theta_hat + p * phi_hat;
  Vector3d PHIlm = -p*theta_hat + t*phi_hat;
  return std::make_pair(PSIlm, PHIlm);
}

inline Vector3d build_vector_from_one_mode(int l, int m, double er, double e1, double e2, double cosTheta, double phi, Vector3d r_hat, Vector3d theta_hat, Vector3d phi_hat)
{
  Vector3d PSIlm, PHIlm;
  double Ylm = realYlm(l, m, cosTheta, phi);
  std::tie(PSIlm, PHIlm) = PSIlmPHIlm(l, m, cosTheta, phi, theta_hat, phi_hat);
  return ( r_hat * (Ylm*er) + PSIlm*e1 + PHIlm*e2 );
}

// inline Vector3d build_vector_from_one_complex_mode(int l, int m, double cr, double c1, double c2, double cosTheta, double phi, Vector3d r_hat, Vector3d theta_hat, Vector3d phi_hat)
// {
//   Vector3d PSIlm, PHIlm;
//   double Ylm = realYlm(l, m, cosTheta, phi);
//   std::tie(PSIlm, PHIlm) = PSIlmPHIlm(l, m, cosTheta, phi, theta_hat, phi_hat);
//   return ( r_hat * (Ylm*er) + PSIlm*e1 + PHIlm*e2 );
// }

inline Vector3d get_purely_rotational_part_from_E2(std::vector<double> E2, double cosTheta, double phi, Vector3d theta_hat, Vector3d phi_hat)
{
  Vector3d PSIlm, PHIlm;
  std::tie(PSIlm, PHIlm) = PSIlmPHIlm(1, -1, cosTheta, phi, theta_hat, phi_hat);
  Vector3d rot_part = PHIlm*E2[1];
  std::tie(PSIlm, PHIlm) = PSIlmPHIlm(1, 0, cosTheta, phi, theta_hat, phi_hat);
  rot_part += PHIlm*E2[2];
  std::tie(PSIlm, PHIlm) = PSIlmPHIlm(1, 1, cosTheta, phi, theta_hat, phi_hat);
  rot_part += PHIlm*E2[3];

  return rot_part;
}

#endif //ifndef
