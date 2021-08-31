//2021-02-27, uploaded by nakai
//2020-10-27, uploaded by yokkaich
//2016-11-22, uploaded by nakai
//2016-04-01, uploaded by nakai
#include "E16ANA_GeometryV2.hh"

#include <fstream>
#include <bitset>

#include <TMath.h>

using namespace std;
using namespace E16Geo2D;

/* -------- E16ANA_PlanarGeometry -------- */

E16ANA_PlanarGeometry::E16ANA_PlanarGeometry(const std::string &type, int module_id, int layer_id) :
   E16ANA_DetectorGeometry(type, module_id, layer_id),
   detector_center(0.,0.,0.), rotation()
{
}

#if 0
E16ANA_PlanarGeometry::E16ANA_PlanarGeometry(const double _xyz[3], const double _rot[3]) :
   detector_center(_xyz[0], _xyz[1], _xyz[2]), rotation()
{
   rotation.RotateX(_rot[0]);
   rotation.RotateY(_rot[1]);
   rotation.RotateZ(_rot[2]);
}

E16ANA_PlanarGeometry::E16ANA_PlanarGeometry(const double _xyz[3], const double _rot[3], const double _dxyz[3], const double _drot[3]) :
   detector_center(_xyz[0]+_dxyz[0], _xyz[1]+_dxyz[1], _xyz[2]+_dxyz[2]), rotation()
{
   rotation.RotateX(_drot[0]);
   rotation.RotateY(_drot[1]);
   rotation.RotateZ(_drot[2]);

   rotation.RotateX(_rot[0]);
   rotation.RotateY(_rot[1]);
   rotation.RotateZ(_rot[2]);
/*
   for(int i=0; i<3; i++){
      xyz[i] = _xyz[i];
      rot[i] = _rot[i];
      dxyz[i] = _dxyz[i];
      drot[i] = _drot[i];
   }
*/
}
#endif

E16ANA_PlanarGeometry::~E16ANA_PlanarGeometry(){
   //std::cout << "E16ANA_PlanarGeometry delete" << std::endl;
}

TVector3 E16ANA_PlanarGeometry::GetGPos(const TVector3 &lpos) const {
  //  cerr<<"rotation "<<GetRotationG4()<<endl;
   return rotation*lpos+detector_center;
}

TVector3 E16ANA_PlanarGeometry::GetLPos(const TVector3 &gpos) const {
   return rotation.Inverse()*(gpos-detector_center);
}

TVector3 E16ANA_PlanarGeometry::GetGMom(const TVector3 &lmom) const {
   return rotation*lmom;
}

TVector3 E16ANA_PlanarGeometry::GetLMom(const TVector3 &gmom) const {
   return rotation.Inverse()*gmom;
}

ROOT::Math::Plane3D E16ANA_PlanarGeometry::GetPlane(double local_z) const {
   ROOT::Math::Plane3D::Point origin;
   origin = GetGPos(TVector3(0.,0.,local_z));
   TVector3 normal;
   normal.SetXYZ(0.,0.,1.);
   normal = rotation*normal;
   ROOT::Math::Plane3D::Vector normal_temp(normal.X(), normal.Y(), normal.Z());
   ROOT::Math::Plane3D plane(normal_temp, origin);
   return plane;
}

bool E16ANA_PlanarGeometry::IsCrossed(const double gpos0[], const double gpos1[], TVector3 &cross_lpos, double local_z) const {
   ROOT::Math::Plane3D plane = GetPlane(local_z);
   ROOT::Math::Plane3D::Point p0(gpos0[0], gpos0[1], gpos0[2]);
   ROOT::Math::Plane3D::Point p1(gpos1[0], gpos1[1], gpos1[2]);
   double dist0 = plane.Distance(p0);
   double dist1 = plane.Distance(p1);
   if((dist0*dist1) > 0.0){
      return false;
   }
   ROOT::Math::Plane3D::Point cp = p0+(p1-p0)*fabs(dist0)/(fabs(dist0)+fabs(dist1));
   cross_lpos = GetLPos(TVector3(cp.X(), cp.Y(), cp.Z()));
   return true;
}

void E16ANA_PlanarGeometry::Translate(const TVector3 &delta_r){
   detector_center += delta_r;
}

void E16ANA_PlanarGeometry::Rotate(double alpha, double beta, double gamma, const TVector3 &rotation_center){
   TRotation temp_rot;
   temp_rot.RotateX(alpha);
   temp_rot.RotateY(beta);
   temp_rot.RotateZ(gamma);
   detector_center = temp_rot*(detector_center-rotation_center)+rotation_center;
   rotation = temp_rot*rotation;
}

void E16ANA_PlanarGeometry::LocalTranslate(const TVector3 &delta_r){
   TVector3 g_dir = this->GetGMom(delta_r);
   detector_center += g_dir;
}

/* -------- E16ANA_GeometryV2 -------- */

std::ostream& operator<<(std::ostream &os, const E16ANA_GeometryV2::Params_t &rhs){
   for(int i=0; i<3; i++){
      os << " " << std::setw(12) << rhs.xyz[i];
   }
   for(int i=0; i<3; i++){
      os << " " << std::setw(12) << rhs.rot[i]/TMath::Pi()*180.0;
   }
   for(int i=0; i<3; i++){
      os << " " << std::setw(12) << rhs.rot_center[i];
   }
   return os;
}

std::istream& operator>>(std::istream &is, E16ANA_GeometryV2::Params_t &rhs){
   double temp;
   for(int i=0; i<3; i++){
      is >> temp;
      rhs.xyz[i] = temp;
   }
   for(int i=0; i<3; i++){
      is >> temp;
      rhs.rot[i] = temp*(TMath::Pi()/180.0);
   }
   for(int i=0; i<3; i++){
      is >> temp;
      rhs.rot_center[i] = temp;
   }
   return is;
}

E16ANA_GeometryV2::Node_t::Node_t(const std::string &_name) : 
   name(_name)
{
}

E16ANA_GeometryV2::Node_t::~Node_t(){
   //std::cout << "E16ANA_GeometryV2::Node_t delete, name = " << name << std::endl;
   std::map<std::string, Node_t*>::iterator it;
   for(it = children.begin(); it != children.end(); it++){
      Node_t *node = it->second;
      if(node != NULL){
         delete node;
      }
   }
}

//void E16ANA_GeometryV2::Node_t::AddChild(const std::string &key, Node_t *node){
void E16ANA_GeometryV2::Node_t::AddChild(Node_t *node){
   std::string key = node->GetName();
   children[key] = node;
}

void E16ANA_GeometryV2::Node_t::Translate(const TVector3 &delta_r){
   std::map<std::string, Node_t*>::iterator it;
   for(it = children.begin(); it != children.end(); it++){
      Node_t *node = it->second;
      node->Translate(delta_r);
   }
}

void E16ANA_GeometryV2::Node_t::Rotate(double alpha, double beta, double gamma, const TVector3 &rotation_center){
   std::map<std::string, Node_t*>::iterator it;
   for(it = children.begin(); it != children.end(); it++){
      Node_t *node = it->second;
      node->Rotate(alpha, beta, gamma, rotation_center);
   }
}

void E16ANA_GeometryV2::Node_t::Rotate(double alpha, double beta, double gamma){
   std::map<std::string, Node_t*>::iterator it;
   for(it = children.begin(); it != children.end(); it++){
      Node_t *node = it->second;
      node->Rotate(alpha, beta, gamma);
   }
}

void E16ANA_GeometryV2::Node_t::Reset(){
   std::map<std::string, Node_t*>::iterator it;
   for(it = children.begin(); it != children.end(); it++){
      Node_t *node = it->second;
      node->Reset();
   }
}

E16ANA_GeometryV2::Node_t* E16ANA_GeometryV2::Node_t::FindChild(const std::string &key){
   std::map<std::string, Node_t*>::iterator it = children.find(key);
   if(it == children.end()){
      return NULL;
   }else{
      return it->second;
   }
}

void E16ANA_GeometryV2::Node_t::TranslateAndRotate(const Params_t &_par){
   this->Rotate(_par.rot, _par.rot_center);
   this->Translate(_par.xyz);
}

void E16ANA_GeometryV2::Node_t::TranslateAndRotate(const std::string &tree_str, ParamsArray_t &_params){
   std::map<std::string, Node_t*>::iterator it;
   for(it = children.begin(); it != children.end(); it++){
      Node_t *node = it->second;
      node->TranslateAndRotate(tree_str+name+delimiter, _params);
   }
   std::string key = design_str+delimiter+tree_str+name;
   std::string key_error = error_str+delimiter+tree_str+name;

   // operation order is error -> design
   this->TranslateAndRotate(_params[key_error]);
   this->TranslateAndRotate(_params[key]);
}

void E16ANA_GeometryV2::Node_t::Print(int depth){
   std::cout << "E16ANA_GeometryV2::Node_t::Print() : " << std::setw((depth)*4) << "" << name <<std::endl;
   std::map<std::string, Node_t*>::iterator it;
   for(it = children.begin(); it != children.end(); it++){
      Node_t *node = it->second;
      node->Print(depth+1);
   }
}

E16ANA_GeometryV2::Leaf_t::Leaf_t(const std::string &_name, E16ANA_DetectorGeometry *_detector) : 
   Node_t(_name), detector(_detector)
{
}

E16ANA_GeometryV2::Leaf_t::~Leaf_t(){
   //std::cout << "E16ANA_GeometryV2::Leaf_t delete" << std::endl;
   if(detector != NULL){
      delete detector;
   }
}

const std::string E16ANA_GeometryV2::delimiter   = ":";
const std::string E16ANA_GeometryV2::design_str  = "D"; // Design
const std::string E16ANA_GeometryV2::error_str   = "A"; // Adjust
const std::string E16ANA_GeometryV2::e16_str     = "E16";
const std::string E16ANA_GeometryV2::gtr_str     = "GTR";
const std::string E16ANA_GeometryV2::hbd_str     = "HBD";
const std::string E16ANA_GeometryV2::ssd_str     = "SSD";
const std::string E16ANA_GeometryV2::lgvd_str    = "LGVD";
const std::string E16ANA_GeometryV2::lg_str      = "LG";
const std::string E16ANA_GeometryV2::frame_str   = "Frame";
const std::string E16ANA_GeometryV2::layer_str   = "Layer";
const std::string E16ANA_GeometryV2::chamber_str = "Chamber";
const std::string E16ANA_GeometryV2::module_str  = "Module";
const std::string E16ANA_GeometryV2::block_str   = "Block";

E16ANA_GeometryV2::E16ANA_GeometryV2() : 
  n_gtr_modules(33), n_gtr_layers(3), n_ssd_modules(33), n_ssd_layers(1),n_hbd_modules(27), n_lg_modules(27), n_lg_blocks(60)
{
   Initialize();
   SetDesignValues();
   std::cerr << "E16ANA_GeometryV2 initialized by design values." << std::endl;
}

E16ANA_GeometryV2::E16ANA_GeometryV2(const std::string &filename) :
   n_gtr_modules(33), n_gtr_layers(3), n_ssd_modules(33), n_ssd_layers(1), n_hbd_modules(27), n_lg_modules(27), n_lg_blocks(60)
{
   std::ifstream ifs(filename);
   if(ifs.fail()){
      std::cerr << __func__ << " (" << __FILE__ << "/" << __LINE__ << ") exit: geom file cannot open...: " << filename << std::endl;
      exit(1);
   }
   std::cerr << "E16ANA_GeometryV2 initialized by " << filename << std::endl;

   ReadGeometryFile(filename, params);
   Initialize();
   SetValues(geometry_tree, params);
}

void E16ANA_GeometryV2::Initialize(){

   // design values are used for testing purposes only
   design_gtr_theta[ 0] =  121.53/180.0*TMath::Pi();
   design_gtr_theta[ 1] =   97.87/180.0*TMath::Pi();
   design_gtr_theta[ 2] =   74.21/180.0*TMath::Pi();
   design_gtr_theta[ 3] =   50.55/180.0*TMath::Pi();
   design_gtr_theta[ 4] =   26.89/180.0*TMath::Pi();
   design_gtr_theta[ 5] =    0.00/180.0*TMath::Pi();
   design_gtr_theta[ 6] =  -26.89/180.0*TMath::Pi();
   design_gtr_theta[ 7] =  -50.55/180.0*TMath::Pi();
   design_gtr_theta[ 8] =  -74.21/180.0*TMath::Pi();
   design_gtr_theta[ 9] =  -97.87/180.0*TMath::Pi();
   design_gtr_theta[10] = -121.53/180.0*TMath::Pi();

   // design_gtr_XX[stage][layer]
   //                   Layer-1                       Layer-2                       Layer-3
   design_gtr_ya[0][0] = -137.6; design_gtr_ya[0][1] = -262.0; design_gtr_ya[0][2] = -384.6; // Lower
   design_gtr_ya[1][0] =    0.0; design_gtr_ya[1][1] =    0.0; design_gtr_ya[1][2] =    0.0; // Middle
   design_gtr_ya[2][0] =  137.6; design_gtr_ya[2][1] =  262.0; design_gtr_ya[2][2] =  384.6; // Upper

   design_gtr_yb[0][0] = -125.0; design_gtr_yb[0][1] = -242.6; design_gtr_yb[0][2] = -363.1;
   design_gtr_yb[1][0] =    0.0; design_gtr_yb[1][1] =    0.0; design_gtr_yb[1][2] =    0.0;
   design_gtr_yb[2][0] =  125.0; design_gtr_yb[2][1] =  242.6; design_gtr_yb[2][2] =  363.1;

   design_gtr_za[0][0] =  237.5; design_gtr_za[0][1] =  450.3; design_gtr_za[0][2] =  636.2;
   design_gtr_za[1][0] =  237.5; design_gtr_za[1][1] =  450.3; design_gtr_za[1][2] =  636.2;
   design_gtr_za[2][0] =  237.5; design_gtr_za[2][1] =  450.3; design_gtr_za[2][2] =  636.2;

   design_gtr_zb[0][0] =  197.7; design_gtr_zb[0][1] =  404.9; design_gtr_zb[0][2] =  582.7;
   design_gtr_zb[1][0] =  197.7; design_gtr_zb[1][1] =  404.9; design_gtr_zb[1][2] =  582.7;
   design_gtr_zb[2][0] =  197.7; design_gtr_zb[2][1] =  404.9; design_gtr_zb[2][2] =  582.7;

   ConstructGeometryTree();
   //SetValues(geometry_tree, params);
}

E16ANA_GeometryV2::~E16ANA_GeometryV2(){
   if(geometry_tree != NULL){
      delete geometry_tree;
   }
}

void E16ANA_GeometryV2::ConstructGeometryTree(){
   gtr_geometry = new E16ANA_DetectorGeometry**[n_gtr_modules];
   for(int i=0; i<n_gtr_modules; i++){
      gtr_geometry[i] = new E16ANA_DetectorGeometry*[n_gtr_layers];
   }
   ssd_geometry = new E16ANA_DetectorGeometry**[n_ssd_modules];
   for(int i=0; i<n_ssd_modules; i++){
      ssd_geometry[i] = new E16ANA_DetectorGeometry*[n_ssd_layers];
   }
   hbd_geometry = new E16ANA_DetectorGeometry*[n_hbd_modules];
   lgvd_geometry = new E16ANA_DetectorGeometry*[n_lg_modules];
   lg_geometry = new E16ANA_DetectorGeometry**[n_lg_modules];
   for(int i=0; i<n_lg_modules; i++){
      lg_geometry[i] = new E16ANA_DetectorGeometry*[n_lg_blocks];
   }
   geometry_tree = new Node_t(e16_str);
   Node_t *gtr = new Node_t(gtr_str);
   Node_t *hbd = new Node_t(hbd_str);
   Node_t *ssd = new Node_t(ssd_str);
   Node_t *lgvd = new Node_t(lgvd_str);
   Node_t *lg = new Node_t(lg_str);
   Leaf_t *gtr_leaf[n_gtr_modules][n_gtr_layers];
   Leaf_t *hbd_leaf[n_hbd_modules];
   Leaf_t *ssd_leaf[n_ssd_modules][n_ssd_layers];
   Leaf_t *lgvd_leaf[n_lg_modules];
   Leaf_t *lg_leaf[n_lg_modules][n_lg_blocks];
   geometry_tree->AddChild(gtr);
   geometry_tree->AddChild(hbd);
   geometry_tree->AddChild(ssd);
   geometry_tree->AddChild(lgvd);
   geometry_tree->AddChild(lg);

   std::string node_name;
   // construct detector geometry (each chamber)
   for(int i=0; i<n_gtr_modules; i++){
      node_name = GetNumString(chamber_str, i);
      for(int j=0; j<n_gtr_layers; j++){
         gtr_geometry[i][j] = new E16ANA_PlanarGeometry(gtr_str,i,j);
         gtr_leaf[i][j] = new Leaf_t(node_name, gtr_geometry[i][j]);
      }
   }
   for(int i=0; i<n_hbd_modules; i++){
      node_name = GetNumString(chamber_str, i);
      hbd_geometry[i] = new E16ANA_PlanarGeometry(hbd_str,i,0);
      hbd_leaf[i] = new Leaf_t(node_name, hbd_geometry[i]);
   }
#if 0
   for(int i=0; i<n_ssd_modules; i++){
      node_name = GetNumString(chamber_str, i);
      ssd_geometry[i] = new E16ANA_PlanarGeometry(ssd_str,i,0);
      ssd_leaf[i] = new Leaf_t(node_name, ssd_geometry[i]);
   }
#else
   for(int i=0; i<n_ssd_modules; i++){
      node_name = GetNumString(chamber_str, i);
      for(int j=0; j<n_ssd_layers; j++){
         ssd_geometry[i][j] = new E16ANA_PlanarGeometry(ssd_str,i,j);
         ssd_leaf[i][j] = new Leaf_t(node_name, ssd_geometry[i][j]);
      }
   }
#endif
   for(int i=0; i<n_lg_modules; i++){
      node_name = GetNumString(chamber_str, i);
      lgvd_geometry[i] = new E16ANA_PlanarGeometry(lgvd_str,i,0);
      lgvd_leaf[i] = new Leaf_t(node_name, lgvd_geometry[i]);
   }
   for(int i=0; i<n_lg_modules; i++){
      for(int j=0; j<n_lg_blocks; j++){
         node_name = GetNumString(block_str, j);
         lg_geometry[i][j] = new E16ANA_PlanarGeometry(lg_str,i,j);
         lg_leaf[i][j] = new Leaf_t(node_name, lg_geometry[i][j]);
      }
   }

   // construct GTR tree
   for(int i=0; i<n_gtr_modules/3; i++){
      node_name = GetNumString(frame_str, i);
      Node_t *gtr_frame = new Node_t(node_name);
      for(int j=0; j<n_gtr_layers; j++){
         node_name = GetNumString(layer_str, j);
         Node_t *gtr_ladder = new Node_t(node_name);
         gtr_ladder->AddChild(gtr_leaf[i*3  ][j]);
         gtr_ladder->AddChild(gtr_leaf[i*3+1][j]);
         gtr_ladder->AddChild(gtr_leaf[i*3+2][j]);
         gtr_frame->AddChild(gtr_ladder);
      }
      gtr->AddChild(gtr_frame);
   }

   // construct HBD tree
   for(int i=0; i<n_hbd_modules/3; i++){
      node_name = GetNumString(frame_str, i);
      Node_t *hbd_frame = new Node_t(node_name);
      hbd_frame->AddChild(hbd_leaf[i*3  ]);
      hbd_frame->AddChild(hbd_leaf[i*3+1]);
      hbd_frame->AddChild(hbd_leaf[i*3+2]);
      hbd->AddChild(hbd_frame);
   }

   // construct SSD tree
#if 0
   for(int i=0; i<n_ssd_modules/3; i++){
      node_name = GetNumString(frame_str, i);
      Node_t *ssd_frame = new Node_t(node_name);
      ssd_frame->AddChild(ssd_leaf[i*3  ]);
      ssd_frame->AddChild(ssd_leaf[i*3+1]);
      ssd_frame->AddChild(ssd_leaf[i*3+2]);
      ssd->AddChild(ssd_frame);
   }
#else
   for(int i=0; i<n_ssd_modules/3; i++){
      node_name = GetNumString(frame_str, i);
      Node_t *ssd_frame = new Node_t(node_name);
      for(int j=0; j<n_ssd_layers; j++){
         node_name = GetNumString(layer_str, j);
         Node_t *ssd_ladder = new Node_t(node_name);
         ssd_ladder->AddChild(ssd_leaf[i*3  ][j]);
         ssd_ladder->AddChild(ssd_leaf[i*3+1][j]);
         ssd_ladder->AddChild(ssd_leaf[i*3+2][j]);
         ssd_frame->AddChild(ssd_ladder);
      }
      ssd->AddChild(ssd_frame);
   }
#endif

   // construct LGVD tree
   for(int i=0; i<n_lg_modules/3; i++){
      node_name = GetNumString(frame_str, i);
      Node_t *lgvd_frame = new Node_t(node_name);
      lgvd_frame->AddChild(lgvd_leaf[i*3  ]);
      lgvd_frame->AddChild(lgvd_leaf[i*3+1]);
      lgvd_frame->AddChild(lgvd_leaf[i*3+2]);
      lgvd->AddChild(lgvd_frame);
   }

   // construct LG tree
   std::bitset<64> valid_lg_block; // TODO : more sophisticated impl
   for(int i=0; i<6; i++){
      for(int j=0; j<6; j++){
         valid_lg_block.set(i*10 + j);
      }
   }
   valid_lg_block.set(26);   valid_lg_block.set(36);
   for(int i=0; i<n_lg_modules/3; i++){
      node_name = GetNumString(frame_str, i);
      Node_t *lg_frame = new Node_t(node_name);
      for(int j=0; j<3; j++){
         int module_id = i*3+j;
         node_name = GetNumString(module_str, module_id);
         Node_t *lg_module = new Node_t(node_name);
         for(int k=0; k<n_lg_blocks; k++){
            if (valid_lg_block[k]) {
               lg_module->AddChild(lg_leaf[module_id][k]);
            }
         }
         lg_frame->AddChild(lg_module);
      }
      lg->AddChild(lg_frame);
   }

   //geometry_tree->Print();
}

void E16ANA_GeometryV2::SetDesignValues(){
   params.clear();
   geometry_tree->Reset();
   // GTR local translate
   for(int i=0; i<n_gtr_modules; i++){
      int module = i;
      int frame = module/3;
      int is_a = (frame+1)%2;
      int stage = module%3;
      for(int j=0; j<n_gtr_layers; j++){
         int layer = j;
         double y, z;
         if(is_a){
            y = design_gtr_ya[stage][layer];
            z = design_gtr_za[stage][layer];
         }else{
            y = design_gtr_yb[stage][layer];
            z = design_gtr_zb[stage][layer];
         }
         //std::cout << GetTreeStringGTR(frame, layer, module) << " : " << y << ", " << z << std::endl;
         params[GetTreeStringGTR(frame, layer, module)].SetXYZ(0.,y,z);
      }
   }
   // GTR frame rotation
   for(int i=0; i<n_gtr_modules/3; i++){
      //std::cout << GetTreeStringGTR(i) << " : " << design_gtr_theta[i] << std::endl;
      params[GetTreeStringGTR(i)].SetRotation(0.,design_gtr_theta[i],0.);
   }
   SetValues(geometry_tree, params);
}

void E16ANA_GeometryV2::SetValuesFromV1(E16ANA_GeometryV1 *geom_v1){
   params.clear();
   geometry_tree->Reset();

   // GTR local
   for(int i=0; i<n_gtr_modules; i++){
      int module = i;
      int frame = module/3;
      for(int j=0; j<n_gtr_layers; j++){
         int layer = j;
         double y = geom_v1->GTRy[layer][module]*10.0;
         double z = geom_v1->GTRz[layer][module]*10.0;
         double dx = geom_v1->GTRdx[layer][module]*10.0;
         double dy = geom_v1->GTRdy[layer][module]*10.0;
         double dz = geom_v1->GTRdz[layer][module]*10.0;
         double drotx = geom_v1->GTRrotx[layer][module]/180.0*TMath::Pi();
         double droty = geom_v1->GTRroty[layer][module]/180.0*TMath::Pi();
         double drotz = geom_v1->GTRrotz[layer][module]/180.0*TMath::Pi();
         //params[GetTreeStringGTR(frame, layer, module)].SetXYZ(0.,y,z);
         //params[GetTreeStringGTR(frame, layer, module)+delimiter+error_str].SetXYZ(dx,dy,dz);
         //params[GetTreeStringGTR(frame, layer, module)+delimiter+error_str].SetRotation(drotx,droty,drotz);
         params[design_str+delimiter+GetTreeStringGTR(frame, layer, module)].SetXYZ(0.,y,z);
         params[error_str+delimiter+GetTreeStringGTR(frame, layer, module)].SetXYZ(dx,dy,dz);
         params[error_str+delimiter+GetTreeStringGTR(frame, layer, module)].SetRotation(drotx,droty,drotz);

      }
   }
   // GTR frame
   for(int i=0; i<n_gtr_modules/3; i++){
      double roty = geom_v1->GTRtheta[i*3]/180.0*TMath::Pi();
      double dx = geom_v1->GTRdx_frame[i*3]*10.0;
      double dy = geom_v1->GTRdy_frame[i*3]*10.0;
      double dz = geom_v1->GTRdz_frame[i*3]*10.0;
      double drotx = geom_v1->GTRrotx_frame[i*3]/180.0*TMath::Pi();
      double droty = geom_v1->GTRroty_frame[i*3]/180.0*TMath::Pi();
      double drotz = geom_v1->GTRrotz_frame[i*3]/180.0*TMath::Pi();
      double cx = geom_v1->GTRcx_frame[i*3]*10.0;
      double cy = geom_v1->GTRcy_frame[i*3]*10.0;
      double cz = geom_v1->GTRcz_frame[i*3]*10.0;
      params[design_str+delimiter+GetTreeStringGTR(i)].SetRotation(0.,roty,0.);
      params[error_str+delimiter+GetTreeStringGTR(i)].SetXYZ(dx,dy,dz);
      params[error_str+delimiter+GetTreeStringGTR(i)].SetRotation(drotx,droty,drotz);
      params[error_str+delimiter+GetTreeStringGTR(i)].SetRotationCenter(cx,cy,cz);
   }
   // HBD
   for(int i=0; i<n_hbd_modules; i++){
      double y = geom_v1->HBDy[i]*10.0;
      double z = geom_v1->HBDz[i]*10.0;
      double roty = 30.76*(4-i/3)/180.0*TMath::Pi();
      params[design_str+delimiter+GetTreeStringHBD(i/3,i)].SetXYZ(0.,y,z);
      params[design_str+delimiter+GetTreeStringHBD(i/3)].SetRotation(0.,roty,0.);
   }
   // LGVD
   for(int i=0; i<n_lg_modules; i++){
      double y = 905.5*(i%3-1);
      double z = 1400.0;
      double roty = 30.76*(4-i/3)/180.0*TMath::Pi();
      params[design_str+delimiter+GetTreeStringLGVD(i/3,i)].SetXYZ(0.,y,z);
      params[design_str+delimiter+GetTreeStringLGVD(i/3)].SetRotation(0.,roty,0.);
   }
   // SSD
   for(int i=0; i<n_ssd_modules; i++){
      int id_offset = (n_gtr_modules-n_ssd_modules)/2;
      double y = geom_v1->GTRy[0][i+id_offset]*10.0/2.0;
      double z = geom_v1->GTRz[0][i+id_offset]*10.0/2.0;
      double roty = geom_v1->GTRtheta[i+id_offset]/180.0*TMath::Pi();
      params[design_str+delimiter+GetTreeStringSSD(i/3,i)].SetXYZ(0.,y,z);
      params[design_str+delimiter+GetTreeStringSSD(i/3)].SetRotation(0.,roty,0.);
   }

   SetValues(geometry_tree, params);
}

void E16ANA_GeometryV2::SetValues(Node_t *tree_root, ParamsArray_t &_params){
   tree_root->TranslateAndRotate("", _params);
}

void E16ANA_GeometryV2::PrintParams(std::ostream &ros){
   ParamsArray_t::iterator it;
   int str_len = (GetTreeStringGTR(0,0,0)+delimiter+error_str).length();
   for(it = params.begin(); it != params.end(); it++){
      ros << std::setw(str_len) << std::left << it->first << std::right << it->second << std::endl;
   }
}

void E16ANA_GeometryV2::WriteGeometryFile(const std::string &filename){
   std::ofstream ofs(filename);
   //std::ofstream ofs(filename.c_str());
   PrintParams(ofs);
}

void E16ANA_GeometryV2::ReadGeometryFile(const std::string &filename, ParamsArray_t &_params){
   _params.clear();
   std::ifstream ifs(filename);
   while(ifs.good()){
      std::string buf;
      std::getline(ifs, buf);
      if(buf[0] == '#' || buf.empty()){
         continue;
      }
      std::istringstream is(buf);
      std::string key;
      is >> key;
      is >> _params[key];
   }
   //SetValues(geometry_tree, _params);
}

std::string E16ANA_GeometryV2::GetTreeStringGTR(int frame){
   return e16_str+delimiter
      +gtr_str+delimiter
      +GetNumString(frame_str, frame);
}

std::string E16ANA_GeometryV2::GetTreeStringGTR(int frame, int layer){
   return GetTreeStringGTR(frame)+delimiter+GetNumString(layer_str, layer);
}

std::string E16ANA_GeometryV2::GetTreeStringGTR(int frame, int layer, int module){
   return GetTreeStringGTR(frame, layer)+delimiter+GetNumString(chamber_str, module);
}

std::string E16ANA_GeometryV2::GetTreeStringHBD(int frame){
   return e16_str+delimiter
      +hbd_str+delimiter
      +GetNumString(frame_str, frame);
}

std::string E16ANA_GeometryV2::GetTreeStringHBD(int frame, int module){
   return GetTreeStringHBD(frame)+delimiter+GetNumString(chamber_str, module);
}

std::string E16ANA_GeometryV2::GetTreeStringLGVD(int frame){
   return e16_str+delimiter
      +lgvd_str+delimiter
      +GetNumString(frame_str, frame);
}

std::string E16ANA_GeometryV2::GetTreeStringLGVD(int frame, int module){
   return GetTreeStringLGVD(frame)+delimiter+GetNumString(chamber_str, module);
}

std::string E16ANA_GeometryV2::GetTreeStringLG(int frame){
   return e16_str+delimiter
      +lg_str+delimiter
      +GetNumString(frame_str, frame);
}

std::string E16ANA_GeometryV2::GetTreeStringLG(int frame, int module){
   return GetTreeStringLG(frame)+delimiter+GetNumString(module_str, module);
}

std::string E16ANA_GeometryV2::GetTreeStringLG(int frame, int module, int block){
   return GetTreeStringLG(frame, module)+delimiter+GetNumString(block_str, block);
}

std::string E16ANA_GeometryV2::GetTreeStringSSD(int frame){
   return e16_str+delimiter
      +ssd_str+delimiter
      +GetNumString(frame_str, frame);
}

std::string E16ANA_GeometryV2::GetTreeStringSSD(int frame, int module){
   return GetTreeStringSSD(frame)+delimiter+GetNumString(chamber_str, module);
}

