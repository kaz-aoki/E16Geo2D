//2016-05-02, uploaded by nakai
//2016-04-01, uploaded by nakai
//2015-11-14, uploaded by nakai
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
#ifndef E16Geo2D_E16ANA_Geometry_hh
#define E16Geo2D_E16ANA_Geometry_hh

#include <TVector3.h>
#include <Math/Plane3D.h>

#include "E16ANA_Transform.hh"

namespace E16Geo2D{

class E16ANA_Geometry {
public:
   virtual ~E16ANA_Geometry(){};
   virtual TVector3 GetLPos(const TVector3 &gPos, int layer_id, int module_id) = 0;
   virtual TVector3 GetGPos(const TVector3 &lPos, int layer_id, int module_id) = 0;
};

class E16ANA_GeometryV1 : public E16ANA_Geometry
{
public:
   E16ANA_GeometryV1(char* filename);
   ~E16ANA_GeometryV1();
   void ReadParam();
   /*TVector3 GetLPos(const TVector3 &gPos, double theta, double phi, double d);
   TVector3 GetGPos(const TVector3 &lPos, double theta, double phi, double d) const;
   TVector3 GetLMom(const TVector3 &gMom, double theta, double phi);
   TVector3 GetGMom(const TVector3 &lMom, double theta, double phi) const;*/
   double GTRy[3][33];
   double GTRz[3][33];
   double GTRtheta[33];
   double HBDy[27];
   double HBDz[27];
   double GTRdx[3][33];
   double GTRdy[3][33];
   double GTRdz[3][33];
   double GTRrotx[3][33];
   double GTRroty[3][33];
   double GTRrotz[3][33];
   double GTRdx_frame[33];
   double GTRdy_frame[33];
   double GTRdz_frame[33];
   double GTRcx_frame[33];
   double GTRcy_frame[33];
   double GTRcz_frame[33];
   double GTRrotx_frame[33];
   double GTRroty_frame[33];
   double GTRrotz_frame[33];
   TVector3 GetDetectorCenter(int layer_id, int module_id);
   TVector3 GetLPos(const TVector3 &gPos, int layer_id, int module_id){
      return E16ANA_Transform::GetLPos(gPos, this, layer_id, module_id);
   };
   TVector3 GetGPos(const TVector3 &lPos, int layer_id, int module_id){
      return E16ANA_Transform::GetGPos(lPos, this, layer_id, module_id);
   };
   ROOT::Math::Plane3D GetPlane(int layer_id, int module_id, double local_z = 0.0);

   void WriteGeometryFile(char *filename);
   void SetRandomGeometryError(int seed);

   void SetGTRdxSmear(double val){
      for(int i=0; i<3; i++){
         for(int j=0; j<3; j++){
            GTRdx_smear[i][j] = val;
         }
      }
   };
   void SetGTRdxSmear(double val, int layer, int stage){
      GTRdx_smear[layer][stage] = val;
   };
   void SetGTRdySmear(double val){
      for(int i=0; i<3; i++){
         for(int j=0; j<3; j++){
            GTRdy_smear[i][j] = val;
         }
      }
   };
   void SetGTRdySmear(double val, int layer, int stage){
      GTRdy_smear[layer][stage] = val;
   };
   void SetGTRdzSmear(double val){
      for(int i=0; i<3; i++){
         for(int j=0; j<3; j++){
            GTRdz_smear[i][j] = val;
         }
      }
   };
   void SetGTRdzSmear(double val, int layer, int stage){
      GTRdz_smear[layer][stage] = val;
   };
   void SetGTRrotxSmear(double val){
      for(int i=0; i<3; i++){
         for(int j=0; j<3; j++){
            GTRrotx_smear[i][j] = val;
         }
      }
   };
   void SetGTRrotySmear(double val){
      for(int i=0; i<3; i++){
         for(int j=0; j<3; j++){
            GTRroty_smear[i][j] = val;
         }
      }
   };
   void SetGTRrotzSmear(double val){
      for(int i=0; i<3; i++){
         for(int j=0; j<3; j++){
            GTRrotz_smear[i][j] = val;
         }
      }
   };
   void SetGTRdxFrameSmear(double val){GTRdx_frame_smear = val;};
   void SetGTRdyFrameSmear(double val){GTRdy_frame_smear = val;};
   void SetGTRdzFrameSmear(double val){GTRdz_frame_smear = val;};
   void SetGTRrotxFrameSmear(double val){GTRrotx_frame_smear = val;};
   void SetGTRrotyFrameSmear(double val){GTRroty_frame_smear = val;};
   void SetGTRrotzFrameSmear(double val){GTRrotz_frame_smear = val;};

private:
   int idTr_, idGEM_;
   char* filename_;
   double GTRdx_smear[3][3];
   double GTRdy_smear[3][3];
   double GTRdz_smear[3][3];
   double GTRrotx_smear[3][3];
   double GTRroty_smear[3][3];
   double GTRrotz_smear[3][3];
   double GTRdx_frame_smear;
   double GTRdy_frame_smear;
   double GTRdz_frame_smear;
   double GTRrotx_frame_smear;
   double GTRroty_frame_smear;
   double GTRrotz_frame_smear;

};

};
#endif
