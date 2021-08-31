//2016-05-02, uploaded by nakai
//2016-04-01, uploaded by nakai
//2015-11-02, uploaded by yokkaich
//2015-11-02, uploaded by komatsu
//2015-08-10, uploaded by nakai
//2015-03-27, uploaded by yokkaich
//2015-01-05, modified by yokkaich
//2015-01-05, modified by yokkaich
//2014-08-27, modified by kawama
//2014-05-07, modified by kawama
//2014-04-30, modified by kawama
//2014-04-24, modified by kawama
//2013-11-14, modified by kawama
//2013-05-15, modified by kawama
//2013-05-13, modified by kawama
#include "E16ANA_Transform.hh"
#include "E16ANA_Geometry.hh"
#include "E16ANA_GeometryV2.hh"
//#include <Geant4/CLHEP/Units/PhysicalConstants.h> 
#include <TVector3.h>
#include <TRotation.h>
#include <iostream> 
using namespace std;
using namespace E16Geo2D;
const double PI=3.14159265358979;
/*** GTR & GEM ID 
// From the beam upstream
                             +y
  //- 0-||- 3-||- 6-||- 9-||-12-||-15-||-18-||-21-||-24-\\
+x||- 1-||- 4-||- 7-||-10-||-13-||-16-||-19-||-22-||-25-||-x
  \\- 2-||- 5-||- 8-||-11-||-14-||-17-||-20-||-23-||-26-//
                             -y
//From the upside
                             +z 
  //- 2-||- 2-||- 2-||- 2-||- 2-||- 2-||- 2-||- 2-||- 2-\\
+x||- 1-||- 1-||- 1-||- 1-||- 1-||- 1-||- 1-||- 1-||- 1-||-x
  \\- 0-||- 0-||- 0-||- 0-||- 0-||- 0-||- 0-||- 0-||- 0-//
                             -z
//consecutive ID number of each GEM is defined as 
// idTr*3+idGEM
***/

 E16ANA_Transform::E16ANA_Transform(){
}

 E16ANA_Transform::~E16ANA_Transform(){
}

TVector3 E16ANA_Transform::GetLPos(
      const TVector3 &gPos,
      const E16ANA_GeometryV1 *geom, int layer_id, int module_id){
   return GetLPos(gPos,
         geom->GTRtheta[module_id], geom->GTRy[layer_id][module_id], geom->GTRz[layer_id][module_id],
         geom->GTRdx[layer_id][module_id], geom->GTRdy[layer_id][module_id], geom->GTRdz[layer_id][module_id],
         geom->GTRrotx[layer_id][module_id], geom->GTRroty[layer_id][module_id], geom->GTRrotz[layer_id][module_id],
         geom->GTRdx_frame[module_id], geom->GTRdy_frame[module_id], geom->GTRdz_frame[module_id],
         geom->GTRrotx_frame[module_id], geom->GTRroty_frame[module_id], geom->GTRrotz_frame[module_id],
         geom->GTRcx_frame[module_id], geom->GTRcy_frame[module_id], geom->GTRcz_frame[module_id]
         );
}

TVector3 E16ANA_Transform::GetGPos(
      const TVector3 &lPos,
      const E16ANA_GeometryV1 *geom, int layer_id, int module_id){
   return GetGPos(lPos,
         geom->GTRtheta[module_id], geom->GTRy[layer_id][module_id], geom->GTRz[layer_id][module_id],
         geom->GTRdx[layer_id][module_id], geom->GTRdy[layer_id][module_id], geom->GTRdz[layer_id][module_id],
         geom->GTRrotx[layer_id][module_id], geom->GTRroty[layer_id][module_id], geom->GTRrotz[layer_id][module_id],
         geom->GTRdx_frame[module_id], geom->GTRdy_frame[module_id], geom->GTRdz_frame[module_id],
         geom->GTRrotx_frame[module_id], geom->GTRroty_frame[module_id], geom->GTRrotz_frame[module_id],
         geom->GTRcx_frame[module_id], geom->GTRcy_frame[module_id], geom->GTRcz_frame[module_id]
         );
}

TVector3 E16ANA_Transform::GetLPos(
      const TVector3 &gPos,
      E16ANA_GeometryV2 *geom, int layer_id, int module_id){
   return geom->GetLPos(gPos, layer_id, module_id);
}

TVector3 E16ANA_Transform::GetGPos(
      const TVector3 &lPos,
      E16ANA_GeometryV2 *geom, int layer_id, int module_id){
   return geom->GetGPos(lPos, layer_id, module_id);
}

TVector3 E16ANA_Transform::
GetLPos(const TVector3 &gPos,  double thetaC, double YC, double ZC, // theta,y,z coordinate of detector center
      double dXL, double dYL, double dZL, // local delta xyz
      double rotXL, double rotYL, double rotZL, // local rotation
      double dXF, double dYF, double dZF, // delta xyz of frame
      double rotXF, double rotYF, double rotZF, // rotation of frame
      double rotCXF, double rotCYF, double rotCZF) // rotation center of frame
{
   thetaC*=PI/180.;
   //phiC*=PI/180.;
   rotXL*=PI/180.;
   rotYL*=PI/180.;
   rotZL*=PI/180.;
   rotXF*=PI/180.;
   rotYF*=PI/180.;
   rotZF*=PI/180.;
   TRotation rotC;
   //rotC.RotateX(phiC);
   rotC.RotateY(thetaC);
   TRotation rotinvC=rotC.Inverse();
   TRotation rotR;
   rotR.RotateX(rotXL);
   rotR.RotateY(rotYL);
   rotR.RotateZ(rotZL);
   TRotation rotinvR=rotR.Inverse();
   //TVector3 cent= rotC*TVector3(0,YC,ZC);
   TRotation rotRF;
   rotRF.RotateX(rotXF);
   rotRF.RotateY(rotYF);
   rotRF.RotateZ(rotZF);
   TRotation rotinvRF=rotRF.Inverse();
   TVector3 cent= TVector3(0,YC,ZC);
   TVector3 dxyz= TVector3(dXL,dYL,dZL);
   TVector3 fxyz = TVector3(dXF,dYF,dZF);
   TVector3 fcent = TVector3(rotCXF,rotCYF,rotCZF);
   //cout << cent.X() <<" "<<cent.Y()<<" "<<cent.Z()<<endl;
   //return rotinvR*(rotinvC*(gPos-cent)+dxyz);
   //return rotinvR*(rotinvC*gPos-cent)-dxyz;
   //return rotinvC*rotinvR*(gPos-cent);
   TVector3 temp = rotinvC*gPos-fxyz;
   temp = rotinvRF*(temp-fcent)+fcent;
   return rotinvR*(temp-(cent+dxyz));  
}

TVector3 E16ANA_Transform::
GetGPos(const TVector3 &lPos,  double thetaC, double YC, double ZC, // theta,y,z coordinate of detector center
      double dXL, double dYL, double dZL, // local delta xyz
      double rotXL, double rotYL, double rotZL, // local rotation
      double dXF, double dYF, double dZF, // delta xyz of frame
      double rotXF, double rotYF, double rotZF, // rotation of frame
      double rotCXF, double rotCYF, double rotCZF) // rotation center of frame
{
   thetaC*=PI/180.;
   //phiC*=PI/180.;
   rotXL*=PI/180.;
   rotYL*=PI/180.;
   rotZL*=PI/180.;
   rotXF*=PI/180.;
   rotYF*=PI/180.;
   rotZF*=PI/180.;
   TRotation rotC;
   //rotC.RotateX(phiC);
   rotC.RotateY(thetaC);
   //TRotation rotinvC=rotC.Inverse();
   TRotation rotR;
   rotR.RotateX(rotXL);
   rotR.RotateY(rotYL);
   rotR.RotateZ(rotZL);
   TRotation rotRF;
   rotRF.RotateX(rotXF);
   rotRF.RotateY(rotYF);
   rotRF.RotateZ(rotZF);
   //TRotation rotinvR=rotR.Inverse();
   //TVector3 cent= rotC*TVector3(0,YC,ZC);
   TVector3 cent= TVector3(0,YC,ZC);
   TVector3 dxyz= TVector3(dXL,dYL,dZL);
   TVector3 fxyz = TVector3(dXF,dYF,dZF);
   TVector3 fcent = TVector3(rotCXF,rotCYF,rotCZF);
   //cout << cent.X() <<" "<<cent.Y()<<" "<<cent.Z()<<endl;
   //return rotC*(cent+rotR*(lPos+dxyz));
   //return rotC*(rotR*lPos-dxyz)+cent;
   TVector3 temp = rotR*lPos+(cent+dxyz);
   temp = rotRF*(temp-fcent)+fcent;
   return rotC*(temp+fxyz);
}

TVector3 E16ANA_Transform::
GetLPos(const TVector3 &gPos,  double thetaC, double YC, double ZC, 
      double dXL, double dYL, double dZL,
      double rotXL, double rotYL, double rotZL)
{
   thetaC*=PI/180.;
   //phiC*=PI/180.;
   rotXL*=PI/180.;
   rotYL*=PI/180.;
   rotZL*=PI/180.;
   TRotation rotC;
   //rotC.RotateX(phiC);
   rotC.RotateY(thetaC);
   TRotation rotinvC=rotC.Inverse();
   TRotation rotR;
   rotR.RotateX(rotXL);
   rotR.RotateY(rotYL);
   rotR.RotateZ(rotZL);
   TRotation rotinvR=rotR.Inverse();
   //TVector3 cent= rotC*TVector3(0,YC,ZC);
   TVector3 cent= TVector3(0,YC,ZC);
   TVector3 dxyz= TVector3(dXL,dYL,dZL);
   //cout << cent.X() <<" "<<cent.Y()<<" "<<cent.Z()<<endl;
   //return rotinvR*(rotinvC*(gPos-cent)+dxyz);
   //return rotinvR*(rotinvC*gPos-cent)-dxyz;
   //return rotinvC*rotinvR*(gPos-cent);
   return rotinvR*(rotinvC*gPos-(cent+dxyz));
}

TVector3 E16ANA_Transform::
GetGPos(const TVector3 &lPos,  double thetaC, double YC, double ZC, 
      double dXL, double dYL, double dZL,
      double rotXL, double rotYL, double rotZL)
{
   thetaC*=PI/180.;
   //phiC*=PI/180.;
   rotXL*=PI/180.;
   rotYL*=PI/180.;
   rotZL*=PI/180.;
   TRotation rotC;
   //rotC.RotateX(phiC);
   rotC.RotateY(thetaC);
   TRotation rotinvC=rotC.Inverse();
   TRotation rotR;
   rotR.RotateX(rotXL);
   rotR.RotateY(rotYL);
   rotR.RotateZ(rotZL);
   TRotation rotinvR=rotR.Inverse();
   //TVector3 cent= rotC*TVector3(0,YC,ZC);
   TVector3 cent= TVector3(0,YC,ZC);
   TVector3 dxyz= TVector3(dXL,dYL,dZL);
   //cout << cent.X() <<" "<<cent.Y()<<" "<<cent.Z()<<endl;
   //return rotC*(cent+rotR*(lPos+dxyz));
   //return rotC*(rotR*lPos-dxyz)+cent;
   return rotC*(rotR*lPos+(cent+dxyz));
}

TVector3 E16ANA_Transform::
GetLPos(const TVector3 &gPos,  double thetaC, double dy, double dz, 
      double thetaR, double phiR)
{
   thetaC*=PI/180.;
   //phiC*=PI/180.;
   thetaR*=PI/180.;
   phiR*=PI/180.;
   TRotation rotC;
   //rotC.RotateX(phiC);
   rotC.RotateY(thetaC);
   TRotation rotinvC=rotC.Inverse();
   TRotation rotR;
   rotR.RotateX(phiR);
   rotR.RotateY(thetaR);
   TRotation rotinvR=rotR.Inverse();
   TVector3 cent= rotC*TVector3(0,dy,dz);
   //cout << cent.X() <<" "<<cent.Y()<<" "<<cent.Z()<<endl;
   return rotinvR*rotinvC*(gPos-cent);
   //return rotinvC*rotinvR*(gPos-cent);
}


TVector3 E16ANA_Transform::
GetGPos(const TVector3 &lPos, double thetaC, double dy, double dz, 
      double thetaR, double phiR)
{
   thetaC*=PI/180.;
   //phiC*=PI/180.;
   thetaR*=PI/180.;
   phiR*=PI/180.;
   TRotation rotC;
   //rotC.RotateX(phiC);
   rotC.RotateY(thetaC);
   TRotation rotinvC=rotC.Inverse();
   TRotation rotR;
   rotR.RotateX(phiR);
   rotR.RotateY(thetaR);
   TRotation rotinvR=rotR.Inverse();
   TVector3 cent= TVector3(0,dy,dz);
   //cout << cent.X() <<" "<<cent.Y()<<" "<<cent.Z()<<endl;
   return rotC*(cent+rotR*lPos);
}

TVector3 E16ANA_Transform::
GetLMom(const TVector3 &gMom,  double thetaC, double phiC, 
      double rotXL, double rotYL, double rotZL)
{
   thetaC*=PI/180.;
   //phiC*=PI/180.;
   rotXL*=PI/180.;
   rotYL*=PI/180.;
   rotZL*=PI/180.;
   TRotation rotC;
   //rotC.RotateX(phiC);
   rotC.RotateY(thetaC);
   TRotation rotinvC=rotC.Inverse();
   TRotation rotR;
   rotR.RotateX(rotXL);
   rotR.RotateY(rotYL);
   rotR.RotateZ(rotZL);
   TRotation rotinvR=rotR.Inverse();
   //return rotinvC*rotinvR*(gMom);
   TVector3 lMom=rotinvR*rotinvC*(gMom);
   //cout << "thetaC="<<theta<<"phiC="<<phi<<endl;
   //cout << "gMom="<<gMom.X()<<" "<<gMom.Y()<<" "<<gMom.Z()<<endl;
   //cout << "lMom="<<lMom.X()<<" "<<lMom.Y()<<" "<<lMom.Z()<<endl;
   return lMom;
   
}

TVector3 E16ANA_Transform::
GetLMom(const TVector3 &gMom,  double thetaC, double phiC, double thetaR, double phiR)
{
   thetaC*=PI/180.;
   //phiC*=PI/180.;
   thetaR*=PI/180.;
   phiR*=PI/180.;
   TRotation rotC;
   //rotC.RotateX(phiC);
   rotC.RotateY(thetaC);
   TRotation rotinvC=rotC.Inverse();
   TRotation rotR;
   rotR.RotateX(phiR);
   rotR.RotateY(thetaR);
   TRotation rotinvR=rotR.Inverse();
   //return rotinvC*rotinvR*(gMom);
   TVector3 lMom=rotinvR*rotinvC*(gMom);
   //cout << "thetaC="<<theta<<"phiC="<<phi<<endl;
   //cout << "gMom="<<gMom.X()<<" "<<gMom.Y()<<" "<<gMom.Z()<<endl;
   //cout << "lMom="<<lMom.X()<<" "<<lMom.Y()<<" "<<lMom.Z()<<endl;
   return lMom;
   
}

TVector3 E16ANA_Transform::
GetGMom(const TVector3 &lMom,  double thetaC, double phiC, 
      double rotXL, double rotYL, double rotZL)
{
   thetaC*=PI/180.;
   //phiC*=PI/180.;
   rotXL*=PI/180.;
   rotYL*=PI/180.;
   rotZL*=PI/180.;
   TRotation rotC;
   //rotC.RotateX(phiC);
   rotC.RotateY(thetaC);
   TRotation rotinvC=rotC.Inverse();
   TRotation rotR;
   rotR.RotateX(rotXL);
   rotR.RotateY(rotYL);
   rotR.RotateZ(rotZL);
   TRotation rotinvR=rotR.Inverse();
   return rotC*rotR*(lMom);
}

TVector3 E16ANA_Transform::
GetGMom(const TVector3 &lMom,  double thetaC, double phiC, double thetaR, double phiR)
{
   thetaC*=PI/180.;
   //phiC*=PI/180.;
   thetaR*=PI/180.;
   phiR*=PI/180.;
   TRotation rotC;
   //rotC.RotateX(phiC);
   rotC.RotateY(thetaC);
   TRotation rotinvC=rotC.Inverse();
   TRotation rotR;
   rotR.RotateX(phiR);
   rotR.RotateY(thetaR);
   TRotation rotinvR=rotR.Inverse();
   return rotC*rotR*(lMom);
}
