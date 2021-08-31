//2016-05-02, uploaded by nakai
//2016-04-01, uploaded by nakai
//2015-12-31, uploaded by yokkaich
//2015-11-02, uploaded by yokkaich
//2015-11-02, uploaded by komatsu
//2015-08-10, uploaded by nakai
//2015-03-27, uploaded by yokkaich
//2015-01-05, uploaded by yokkaich
//2015-01-05, uploaded by yokkaich
//2014-08-27, uploaded by kawama
//2014-05-07, uploaded by kawama
//2014-04-30, uploaded by kawama
//2014-04-24, uploaded by kawama
//2013-11-14, uploaded by kawama
//2013-05-13, modified by kawama
//2013-05-13, modified by kawama
#ifndef E16Geo2D_E16ANA_Transform_HH
#define E16Geo2D_E16ANA_Transform_HH

#include "TVector3.h"

namespace E16Geo2D {

class E16ANA_GeometryV1;
class E16ANA_GeometryV2;

class E16ANA_Transform 
{
   public:
   static TVector3 GetLPos(const TVector3 &gPos, const E16ANA_GeometryV1 *geom, int layer_id, int module_id);
   static TVector3 GetGPos(const TVector3 &lPos, const E16ANA_GeometryV1 *geom, int layer_id, int module_id);
   // for backward compatibility (should be eliminated)
   static TVector3 GetLPos(const TVector3 &gPos, E16ANA_GeometryV2 *geom, int layer_id, int module_id);
   static TVector3 GetGPos(const TVector3 &lPos, E16ANA_GeometryV2 *geom, int layer_id, int module_id);

static TVector3 GetLPos(const TVector3 &gPos,  double thetaC, double YC, double ZC, // theta,y,z coordinate of detector center
      double dXL, double dYL, double dZL, // local delta xyz
      double rotXL, double rotYL, double rotZL, // local rotation
      double dXF, double dYF, double dZF, // delta xyz of frame
      double rotXF, double rotYF, double rotZF, // rotation of frame
      double rotCXF, double rotCYF, double rotCZF); // rotation center of frame

static TVector3 GetGPos(const TVector3 &lPos,  double thetaC, double YC, double ZC, // theta,y,z coordinate of detector center
      double dXL, double dYL, double dZL, // local delta xyz
      double rotXL, double rotYL, double rotZL, // local rotation
      double dXF, double dYF, double dZF, // delta xyz of frame
      double rotXF, double rotYF, double rotZF, // rotation of frame
      double rotCXF, double rotCYF, double rotCZF); // rotation center of frame

   static TVector3 GetLPos(const TVector3 &gPos,  double thetaC, 
         double YC, double ZC, 
         double dXL, double dYL, double dZL,
         double rotXL, double rotYL, double rotZL);
   static TVector3 GetLPos(const TVector3 &gPos, double thetaC, 
         double dy, double dz, double thetaR, double phiR);
   static TVector3 GetGPos(const TVector3 &gPos,  double thetaC, 
         double YC, double ZC, 
         double dXL, double dYL, double dZL,
         double rotXL, double rotYL, double rotZL);
   static TVector3 GetGPos(const TVector3 &lPos, double thetaC, 
         double dy, double dz, double thetaR, double phiR);
   static TVector3 GetLMom(const TVector3 &gMom, double thetaC, double phiC, 
         double rotXL, double rotYL, double rotZL);
   static TVector3 GetLMom(const TVector3 &gMom, double thetaC, double phiC, 
         double thetaR, double phiR);
   static TVector3 GetGMom(const TVector3 &lMom, double thetaC, double phiC, 
         double rotXL, double rotYL, double rotZL);
   static TVector3 GetGMom(const TVector3 &lMom, double thetaC, double phiC, 
         double thetaR, double phiR);

   private:
   //int idTr_, idGEM_;
   E16ANA_Transform();
   ~E16ANA_Transform();
};

};
#endif
