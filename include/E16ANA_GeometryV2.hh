//2021-02-27, uploaded by nakai
//2020-10-27, uploaded by yokkaich
//2016-11-22, uploaded by nakai
//2016-05-02, uploaded by nakai
//2016-04-01, uploaded by nakai
//E16ANA_GeometryV2.hh 201023 by W.Nakai
//    Last modified at <2020-10-24 20:13:46 >

#ifndef E16Geo2D_E16ANA_GeometryV2_hh
#define E16Geo2D_E16ANA_GeometryV2_hh

#include <string>
#include <sstream>
#include <map>

#include <TMath.h>
#include <TVector3.h>
#include <TRotation.h>
#include <Math/Plane3D.h>

#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>

#include <G4SystemOfUnits.hh>

#include "E16ANA_Geometry.hh"

namespace E16Geo2D{

class E16ANA_DetectorGeometry {
   // Abstract class
   // A specific geometry (e.g. planar, cylindrical, ...) should be derived from this class
public:
   E16ANA_DetectorGeometry(const std::string &type, int _module_id, int _layer_id) :
      detector_type(type), module_id(_module_id), layer_id(_layer_id)
   {};
   virtual ~E16ANA_DetectorGeometry(){};
   virtual TVector3 GetGPos(const TVector3 &lpos) const = 0; // local  -> global (mm)
   virtual TVector3 GetLPos(const TVector3 &gpos) const = 0; // global -> local  (mm)
   virtual TVector3 GetGMom(const TVector3 &lmom) const = 0;
   virtual TVector3 GetLMom(const TVector3 &gmom) const = 0;

   virtual G4ThreeVector GetGPos(const G4ThreeVector &lpos) const {
      TVector3 gpos = GetGPos(TVector3(lpos.x()/mm, lpos.y()/mm, lpos.z()/mm));
      return G4ThreeVector(gpos.X()*mm, gpos.Y()*mm, gpos.Z()*mm);
   }; // local  -> global
   virtual G4ThreeVector GetLPos(const G4ThreeVector &gpos) const {
      TVector3 lpos = GetLPos(TVector3(gpos.x()/mm, gpos.y()/mm, gpos.z()/mm));
      return G4ThreeVector(lpos.X()*mm, lpos.Y()*mm, lpos.Z()*mm);
   }; // global -> local 
   virtual G4ThreeVector GetGMom(const G4ThreeVector &lmom) const {
      TVector3 gmom = GetGMom(TVector3(lmom.x()/GeV, lmom.y()/GeV, lmom.z()/GeV));
      return G4ThreeVector(gmom.X()*GeV, gmom.Y()*GeV, gmom.Z()*GeV);
   };
   virtual G4ThreeVector GetLMom(const G4ThreeVector &gmom) const {
      TVector3 lmom = GetLMom(TVector3(gmom.x()/GeV, gmom.y()/GeV, gmom.z()/GeV));
      return G4ThreeVector(lmom.X()*GeV, lmom.Y()*GeV, lmom.Z()*GeV);
   };

   virtual TVector3 GetDetectorCenter() const = 0;
   virtual G4ThreeVector GetDetectorCenterG4() const = 0;
   virtual TRotation GetRotation() const = 0;
   virtual G4RotationMatrix GetRotationG4() const = 0;

   // translation along global axis
   virtual void Translate(const TVector3 &delta_r) = 0;
   // rotation around global axis
   // alpha           : rotation angle around global-X
   // beta            : rotation angle around global-Y
   // gamma           : rotation angle around global-Z
   // rotation_center : global coordinates of rotatetion center
   virtual void Rotate(double alpha, double beta, double gamma, const TVector3 &rotation_center) = 0;
   virtual void Rotate(double alpha, double beta, double gamma) = 0;
   // translation along local axis
   virtual void LocalTranslate(const TVector3 &delta_r) = 0;

   // cross determination and calculation of cross point
   virtual bool IsCrossed(const G4ThreeVector &gpos0, const G4ThreeVector &gpos1, G4ThreeVector &cross_lpos, double local_z = 0.0) const {
      TVector3 cp_temp;
      bool flag = IsCrossed(TVector3(gpos0.x()/mm, gpos0.y()/mm, gpos0.z()/mm),
                            TVector3(gpos1.x()/mm, gpos1.y()/mm, gpos1.z()/mm),
                            cp_temp, local_z);
      cross_lpos.set(cp_temp.X()*mm, cp_temp.Y()*mm, cp_temp.Z()*mm);
      return flag;

   };
   virtual bool IsCrossed(const TVector3 &gpos0, const TVector3 &gpos1, TVector3 &cross_lpos, double local_z = 0.0) const = 0;
   virtual bool IsCrossed(const double gpos0[3], const double gpos1[3], TVector3 &cross_lpos, double local_z = 0.0) const = 0;

   virtual void Reset(){};
   std::string GetDetectorType() const {return detector_type;};
   int GetModuleID() const {return module_id;};
   int GetLayerID() const {return layer_id;};

protected:
   std::string detector_type; // "GTR", "HBD", or "SSD"
   int module_id;
   int layer_id;

};

class E16ANA_PlanarGeometry : public E16ANA_DetectorGeometry {
   // planar type detector
public:
   E16ANA_PlanarGeometry(const std::string &type, int module_id, int layer_id);
   //E16ANA_PlanarGeometry(const double _xyz[3], const double _rot[3]);
   //E16ANA_PlanarGeometry(const double _xyz[3], const double _rot[3], const double _dxyz[3], const double _drot[3]);
   ~E16ANA_PlanarGeometry();
   TVector3 GetGPos(const TVector3 &lpos) const;
   TVector3 GetLPos(const TVector3 &gpos) const;
   TVector3 GetGMom(const TVector3 &lmom) const;
   TVector3 GetLMom(const TVector3 &gmom) const;
   TVector3 GetDetectorCenter() const {return detector_center;};
   G4ThreeVector GetDetectorCenterG4() const {
      return G4ThreeVector(
            detector_center.X()*mm,
            detector_center.Y()*mm,
            detector_center.Z()*mm
            );
   };
   TRotation GetRotation() const {return rotation;};
   G4RotationMatrix GetRotationG4() const {
      return G4RotationMatrix(
            CLHEP::HepRep3x3(
               rotation.XX(),rotation.XY(),rotation.XZ(),
               rotation.YX(),rotation.YY(),rotation.YZ(),
               rotation.ZX(),rotation.ZY(),rotation.ZZ()
               ));
   };

   // translation along global axis
   void Translate(const TVector3 &delta_r);
   // rotation around global axis
   // alpha           : rotation angle around global-X
   // beta            : rotation angle around global-Y
   // gamma           : rotation angle around global-Z
   // rotation_center : global coordinates of rotateion center
   void Rotate(double alpha, double beta, double gamma, const TVector3 &rotation_center);
   void Rotate(double alpha, double beta, double gamma){
      Rotate(alpha, beta, gamma, TVector3(0,0,0));
   };
   void LocalTranslate(const TVector3 &delta_r);

   ROOT::Math::Plane3D GetPlane(double local_z = 0.0) const;
   bool IsCrossed(const TVector3 &gpos0, const TVector3 &gpos1, TVector3 &cross_lpos, double local_z = 0.0) const {
      double r0[3] = {gpos0.X(), gpos0.Y(), gpos0.Z()};
      double r1[3] = {gpos1.X(), gpos1.Y(), gpos1.Z()};
      return IsCrossed(r0, r1, cross_lpos, local_z);
   };
   bool IsCrossed(const double gpos0[3], const double gpos1[3], TVector3 &cross_lpos, double local_z = 0.0) const;
   void Reset(){
      detector_center.SetXYZ(0,0,0);
      rotation = TRotation();
   };

private:
   // global coordinates of local origin (0,0,0)
   TVector3 detector_center;
   TRotation rotation;

};

class E16ANA_GeometryV2 : public E16ANA_Geometry {
public:
   E16ANA_GeometryV2(); // for testing purposes only
   E16ANA_GeometryV2(const std::string &filename);
   ~E16ANA_GeometryV2();

   void SetDesignValues(); // for testing purposes only

   // for backward compatibility (cm -> cm)
   TVector3 GetGPos(const TVector3 &lpos, int layer_id, int module_id){
      return (gtr_geometry[module_id][layer_id]->GetGPos(lpos*10.0))*0.1;
   };
   TVector3 GetLPos(const TVector3 &gpos, int layer_id, int module_id){
      return (gtr_geometry[module_id][layer_id]->GetLPos(gpos*10.0))*0.1;
   };
   TVector3 GetGMom(const TVector3 &lmom, int layer_id, int module_id){
      return gtr_geometry[module_id][layer_id]->GetGMom(lmom);
   }
   TVector3 GetLMom(const TVector3 &gmom, int layer_id, int module_id){
      return gtr_geometry[module_id][layer_id]->GetLMom(gmom);
   }

   // Ex. geometry->GTR(module,layer)->GetLPos(gpos); (mm -> mm)
   const E16ANA_DetectorGeometry* GTR(int module, int layer) const {
      return gtr_geometry[module][layer];
   };
   const E16ANA_DetectorGeometry* GTR1(int module) const {
      return gtr_geometry[module][0];
   };
   const E16ANA_DetectorGeometry* GTR2(int module) const {
      return gtr_geometry[module][1];
   };
   const E16ANA_DetectorGeometry* GTR3(int module) const {
      return gtr_geometry[module][2];
   };
   const E16ANA_DetectorGeometry* HBD(int module) const {
      return hbd_geometry[module];
   }
  const E16ANA_DetectorGeometry* SSD(int module,int layer) const {
      return ssd_geometry[module][layer];
   };
   const E16ANA_DetectorGeometry* SSD(int module) const {
      return ssd_geometry[module][0];
   };
   const E16ANA_DetectorGeometry* LGVD(int module) const {
      return lgvd_geometry[module];
   };
   const E16ANA_DetectorGeometry* LG(int module, int block) const {
      return lg_geometry[module][block];
   };

   // conversion from V1 to V2
   void SetValuesFromV1(E16ANA_GeometryV1 *geom_v1);

   void PrintParams(std::ostream &ros);

   // write and read txt file
   void WriteGeometryFile(const std::string &filename);
   void ReadGeometryFile(const std::string &filename){
      ReadGeometryFile(filename, params);
   };

   struct Params_t {
      double xyz[3];
      double rot[3];
      double rot_center[3];
      Params_t(){
         for(int i=0; i<3; i++){
            xyz[i] = 0.0;
            rot[i] = 0.0;
            rot_center[i] = 0.0;
         }
      };
      Params_t(const double pars[9]){
         for(int i=0; i<3; i++){
            xyz[i] = pars[i];
            rot[i] = pars[i+3];
            rot_center[i] = pars[i+6];
         }
      };
      Params_t(const Params_t &rhs){
         for(int i=0; i<3; i++){
            xyz[i] = rhs.xyz[i];
            rot[i] = rhs.rot[i];
            rot_center[i] = rhs.rot_center[i];
         }
      };
      void SetXYZ(double x, double y, double z){
         xyz[0] = x;
         xyz[1] = y;
         xyz[2] = z;
      };
      void SetRotation(double x, double y, double z){
         rot[0] = x;
         rot[1] = y;
         rot[2] = z;
      };
      void SetRotationCenter(double x, double y, double z){
         rot_center[0] = x;
         rot_center[1] = y;
         rot_center[2] = z;
      };
   };

   //int GetNumGTRModules(){return n_gtr_modules;};
   //int GetNumHBDModules(){return n_hbd_modules;};
   //int GetNumSSDModules(){return n_ssd_modules;};
   //void SetNumSSDModules(int n){n_ssd_modules = n;};

private:

   typedef std::map<std::string, Params_t> ParamsArray_t;

   // nodes of geometry tree
   class Node_t {
      public:
         Node_t(const std::string &_name);
         virtual ~Node_t();
         void AddChild(Node_t *node);
         virtual void Translate(const double r[3]){
            Translate(TVector3(r[0], r[1], r[2]));
         };
         virtual void Translate(const TVector3 &delta_r);
         virtual void Rotate(const double rot[3], const double cent[3]){
            Rotate(rot[0], rot[1], rot[2], TVector3(cent[0], cent[1], cent[2]));
         };
         virtual void Rotate(double alpha, double beta, double gamma, const TVector3 &rotation_center);
         virtual void Rotate(double alpha, double beta, double gamma);
         virtual void Reset();
         virtual Node_t* FindChild(const std::string &key);
         std::string GetName(){return name;};
         void TranslateAndRotate(const Params_t &_par);
         void TranslateAndRotate(const std::string &tree_str, ParamsArray_t &_params);
         void Print(int depth = 0);
      private:
         std::map<std::string, Node_t*> children;
      protected:
         std::string name;
   };

   // leafs of geometry tree
   class Leaf_t : public Node_t {
      public:
         //Leaf_t(const std::string &_name, double _xyz[3], double _rot[3], double _dxyz[3], double _drot[3]);
         Leaf_t(const std::string &_name, E16ANA_DetectorGeometry *_detector);
         ~Leaf_t();
         void Translate(const TVector3 &delta_r){
            detector->Translate(delta_r);
         };
         void Rotate(double alpha, double beta, double gamma, const TVector3 &rotation_center){
            //std::cout << __func__ << "1 : name = " << name << std::endl;
            detector->Rotate(alpha, beta, gamma, rotation_center);
         };
         void Rotate(double alpha, double beta, double gamma){
            //std::cout << __func__ << "2 : name = " << name << std::endl;
            detector->Rotate(alpha, beta, gamma);
         };
         void Reset(){
            detector->Reset();
         };
         E16ANA_DetectorGeometry* GetDetector(){return detector;};
      private:
         E16ANA_DetectorGeometry *detector;
   };

   int n_gtr_modules;   int n_gtr_layers;
   int n_ssd_modules;   int n_ssd_layers;
   int n_hbd_modules;
   int n_lg_modules;    int n_lg_blocks;

   // character strings for geometry tree
   // Ex. "E16:GTR:Frame05:Layer02:Module16"
   //     "E16:HBD:Module01"
   //     "E16:SSD:Module01"
   const static std::string delimiter;
   const static std::string design_str;
   const static std::string error_str;
   const static std::string e16_str;
   const static std::string gtr_str;
   const static std::string hbd_str;
   const static std::string ssd_str;
   const static std::string lgvd_str;
   const static std::string lg_str;
   const static std::string frame_str;
   const static std::string layer_str;
   const static std::string chamber_str;
   const static std::string module_str;
   const static std::string block_str;

   Node_t *geometry_tree;
   void ConstructGeometryTree();
   /*
   E16ANA_DetectorGeometry *gtr_geometry[n_gtr_modules][n_gtr_layers];
   E16ANA_DetectorGeometry *hbd_geometry[n_hbd_modules];
   //E16ANA_DetectorGeometry *lg_geometry[n_gtr_modules];
   E16ANA_DetectorGeometry *ssd_geometry[n_ssd_modules];
   */
   E16ANA_DetectorGeometry ***gtr_geometry;
   E16ANA_DetectorGeometry **hbd_geometry;
   //E16ANA_DetectorGeometry **lg_geometry;
   E16ANA_DetectorGeometry **lgvd_geometry;
   E16ANA_DetectorGeometry ***lg_geometry;
   E16ANA_DetectorGeometry ***ssd_geometry;

   ParamsArray_t params;

   // design values (for testing purposes only)
   //double design_gtr_theta[n_gtr_modules/3];
   double design_gtr_theta[11];
   // A-frame
   double design_gtr_ya[3][3];
   double design_gtr_za[3][3];
   // B-frame
   double design_gtr_yb[3][3];
   double design_gtr_zb[3][3];

   void Initialize();
   void SetValues(Node_t *tree_root, ParamsArray_t &_params);
   void ReadGeometryFile(const std::string &filename, ParamsArray_t &_params);

   // string utilities for geometry tree
   std::string GetNumString(const std::string &str, int n){
      std::ostringstream sout;
      sout << str << std::setw(2) << std::setfill('0') << n;
      return sout.str();
   };
   std::string GetTreeStringGTR(int frame);
   std::string GetTreeStringGTR(int frame, int layer);
   std::string GetTreeStringGTR(int frame, int layer, int module);
   std::string GetTreeStringHBD(int frame);
   std::string GetTreeStringHBD(int frame ,int module);
   std::string GetTreeStringSSD(int frame);
   std::string GetTreeStringSSD(int frame, int layer);
   std::string GetTreeStringSSD(int frame, int layer, int module);
   std::string GetTreeStringLGVD(int frame);
   std::string GetTreeStringLGVD(int frame ,int module);
   std::string GetTreeStringLG(int frame);
   std::string GetTreeStringLG(int frame, int module);
   std::string GetTreeStringLG(int frame, int module, int block);
};
}
#endif // E16ANA_GeometryV2_hh
