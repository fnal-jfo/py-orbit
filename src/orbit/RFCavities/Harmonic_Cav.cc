#include "Harmonic_Cav.hh"
#include "ParticleMacroSize.hh"
#define FMT_HEADER_ONLY
#include <fmt/format.h>

#include <iostream>
#include <cmath>
#include <cstdio>

#include "Bunch.hh"
#include "OrbitConst.hh"

using namespace OrbitUtils;

// Constructor
Harmonic_Cav::Harmonic_Cav(double ZtoPhi,
                           double dESync,
                           double RFHNum,
                           double RFVoltage,
                           double RFPhase)
  : CppPyWrapper(NULL),
       ZtoPhi_(ZtoPhi),
       dESync_(dESync),
       RFHNum_(RFHNum),
    RFVoltage_(RFVoltage),
      RFPhase_(RFPhase)
{

  correction = no_correction;
  
  if ( pTmp =getenv("HCAV_CORRECTION")) {

     HCAV_CORRECTION = std::string(pTmp);

     if (HCAV_CORRECTION == "POSITION") {
           correction = position_correction; 
     }
     else if (HCAV_CORRECTION == "ENERGY") {
           correction = energy_correction; 
     }
  }

}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

// Destructor
Harmonic_Cav::~Harmonic_Cav()
{}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

void Harmonic_Cav::setZtoPhi(double ZtoPhi)
{
  ZtoPhi_ = ZtoPhi;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

double Harmonic_Cav::getZtoPhi()
{
  return ZtoPhi_;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

void Harmonic_Cav::setdESync(double dESync)
{
  dESync_ = dESync;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

double Harmonic_Cav::getdESync()
{
  return dESync_;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

void Harmonic_Cav::setRFHNum(double RFHNum)
{
  RFHNum = RFHNum_;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

double Harmonic_Cav::getRFHNum()
{
  return RFHNum_;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

void Harmonic_Cav::setRFVoltage(double RFVoltage)
{
  RFVoltage_ = RFVoltage;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

double Harmonic_Cav::getRFVoltage()
{
  return RFVoltage_;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
void Harmonic_Cav::setRFPhase(double RFPhase)
{
  RFPhase_ = RFPhase;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

double Harmonic_Cav::getRFPhase()
{
  return RFPhase_;
}

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

void Harmonic_Cav::trackBunch(Bunch* bunch)
{

  static std::string DEBUG_BCOR; 
  static bool FIRST_TIME     = true;


  double pi =  OrbitConst::PI;
  
  double ZtoPhi    =  ZtoPhi_;
  double RFHNum    =  RFHNum_;
  double RFVoltage =  RFVoltage_;
  double RFPhase   =  (pi * RFPhase_ / 180.0); // this should not be specified in degrees !
 
  bunch->compress();

  SyncPart* syncPart = bunch->getSyncPart();
  double** arr = bunch->coordArr();

  double gma   =  syncPart->getGamma(); //  gamma
  double bta2  =  (1.0-1.0/(gma*gma));  //  beta**2
  double bg2   =  gma*gma-1.0;          //  gma**2-1 = (beta*gamma)**2
  double t     =  syncPart->getTime();
  double mass  =  syncPart->getMass(); 

  double dESync =  bunch->getCharge()*RFVoltage*sin(RFPhase);
  
  double bratio = dESync/(gma*mass);      // dgamma/gamma
         bratio = 1 + bratio/bg2;         // 1 + (dgamma/gamma)*1/(beta**2*gamma***2) = 1 + dbeta/beta = beta2/beta1
	 double pratio = (dESync/(gma*mass));
         pratio = 1 + pratio/bta2;   // 1 + dp/p
         pratio = 1.0/pratio;        // 1 /(1+dpp/p)

  double pos_cor = 1.0;
  double ek_cor  = bratio;
  
  //  ------------------------------------------
  char* pTmp = 0;

  switch (correction) {
     case: no_correction:
           ek_cor  = 1.0;
	   pos_cor = 1.0;
	   std::cout << "DEBUG:  Harmonic_Cav: NO CORRECTION "<< std::endl;
           break;
      case: position_correction:
	   ek_cor  = 1.0;
	   pos_cor = bratio;
	   std::cout << "DEBUG: Harmonic_Cav: APPLYING POSITION_CORRECTION"<< std::endl;
	   break;
      case: energy_correction:
	   ek_cor  = bratio;
	   pos_cor = 1.0;
	   std::cout << "DEBUG:  Harmonic_Cav: APPLYING ENERGY_CORRECTION "<< std::endl;
	   break;
      
    };
  

 for(int i=0; i< bunch->getSize(); ++i) {

    double phase  = (-ZtoPhi * arr[i][4]* RFHNum);


    arr[i][4]     *= pos_cor;
   
    double dERF   = bunch->getCharge()*RFVoltage*sin(RFPhase +phase);
    arr[i][5]      = ek_cor*arr[i][5] + (dERF - dESync); // velocity correction

    //arr[i][5]      = ek_cor*arr[i][5] + bunch->getCharge()*RFVoltage*cos(RFPhase)*phase; // rf linearization       

    arr[i][1]    *= pratio; // adiabatic damping 
    arr[i][3]    *= pratio; // adiabatic damping 
  } 

}
