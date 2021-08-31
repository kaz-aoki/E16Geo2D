//#include <TGeoMatrix.h>
//#include <TGeoManager.h>
#include "E16Geo2DMan.h"
#include <TMath.h>
#include <iostream>
#include <bitset>
#include "TH2D.h"
#include "TLine.h"

#include "E16ANA_GeometryV2.hh"
#include "E16ANA_GeoUtil.hh"

using namespace std;
//using namespace E16Geo;

E16Geo2DMan::E16Geo2DMan():geom_v2(nullptr),frame(nullptr)
{
  //mgr = new TGeoManager("E16Geo","E16 Geometry");
  string geomfilename="/e16/w/r01/aoki/E16/database/geometry/v2/geometry_Run0b_210226_design.dat";
  geom_v2 = new E16Geo2D::E16ANA_GeometryV2(geomfilename);
  //frame = new TH2D("frame","E16",100,-1500,1500,100,-1500,1500);
}

E16Geo2DMan::~E16Geo2DMan(){
  delete geom_v2;
}

void E16Geo2DMan::draw(const char* option){
  if ( !frame ) {
    frame = new TH2D("frame","E16",100,-1500,1500,100,-1500,1500);
  }
  frame->Draw();
  draw_gtr();
  draw_hbd();
}

void E16Geo2DMan::draw_gtr(){
  TLine line;
  const double halfsize[3] = { 50., 100., 150. };
  for( int i = 1;i<=32;i+=3){
    if ( i == 16 ) continue;
    if ( i == 1 ) continue;
    if ( i == 31 ) continue;
    /// GTR1
    TVector3 lpos(-halfsize[0],0.,0.);
    TVector3 pos = geom_v2->GTR1(i)->GetGPos(lpos);
    lpos.SetXYZ(halfsize[0],0.,0.);
    TVector3 pos2 = geom_v2->GTR1(i)->GetGPos(lpos);
    line.DrawLine(pos.X(),pos.Z(),pos2.X(),pos2.Z());

    /// GTR2
    lpos.SetXYZ(-halfsize[1],0.,0.);
    pos = geom_v2->GTR2(i)->GetGPos(lpos);
    lpos.SetXYZ(halfsize[1],0.,0.);
    pos2 = geom_v2->GTR2(i)->GetGPos(lpos);
    line.DrawLine(pos.X(),pos.Z(),pos2.X(),pos2.Z());
    
    // GTR3
    lpos.SetXYZ(-halfsize[2],0.,0.);
    pos = geom_v2->GTR3(i)->GetGPos(lpos);
    lpos.SetXYZ(halfsize[2],0.,0.);
    pos2 = geom_v2->GTR3(i)->GetGPos(lpos);
    line.DrawLine(pos.X(),pos.Z(),pos2.X(),pos2.Z());
  }
}


void E16Geo2DMan::draw_hbd(){
   TLine line;
   const double halfsize = 300.;
   for( int i = 1;i<=26;i+=3){
     if ( i == 13 ) continue;
     TVector3 lpos(-halfsize,0.,0.);
     TVector3 pos = geom_v2->HBD(i)->GetGPos(lpos);
     lpos.SetXYZ(halfsize,0.,0.);
     TVector3 pos2 = geom_v2->HBD(i)->GetGPos(lpos);
     line.DrawLine(pos.X(),pos.Z(),pos2.X(),pos2.Z());
   }
} 
/*

void E16Geo2DMan::build_using_geom()
{
  if ( geom_v2 != nullptr ) {
    cout << "build_using_geom() aloready called. will not do again." << endl;
  }
  string geomfilename="/e16/w/r01/aoki/E16/database/geometry/v2/geometry_Run0b_210226_design.dat";
  geom_v2 = new E16ANA_GeometryV2(geomfilename);
  
  build_world();
  build_gtr(0);
  build_gtr(1);
  build_gtr(2);
  build_hbd();
  build_lg_using_geom();
 
  mgr->CloseGeometry();
}

void E16Geo2DMan::build()
{
  build_world();
  build_gtr(0);
  build_gtr(1);
  build_gtr(2);
  build_hbd();
  build_lg();
 
  mgr->CloseGeometry();
}

void E16Geo2DMan::build_world()
{
  E16MedMan* med_man = E16MedMan::get_instance();
  TGeoMedium* med = med_man->get("Vacuum");
  top = mgr->MakeBox("TOP",med,500.*E16Geo::cm,500.*E16Geo::cm,500.*E16Geo::cm);
  mgr->SetTopVolume(top);
}


void E16Geo2DMan::build_gtr(int igtr)
{
  E16MedMan* med_man = E16MedMan::get_instance();
  TGeoMedium* med = med_man->get("GTR");
  med->Print();

  double size = 10.*(igtr+1)*E16Geo::cm;
  double thick = 0.4*E16Geo::cm;

  TString str; str.Form("gtr%d",igtr+1);

  TGeoVolume* gtr = mgr->MakeBox(str,med,size/2.,size/2.,thick/2.);
  //TGeoVolume* gtr2 = mgr->MakeBox("gtr2",med,size/2.,size/2.,thick/2.);
  if (igtr == 1 ) gtr->SetLineColor(kCyan);
  gtr->SetLineColor(kBlue);
  gtr->SetTransparency(50);
  double GTR_R1[] = { 23.7*E16Geo::cm, 45.0*E16Geo::cm, 63.6*E16Geo::cm };
  double GTR_R2[] = { 19.5*E16Geo::cm, 40.5*E16Geo::cm, 58.3*E16Geo::cm };

  double R1 = GTR_R1[igtr];
  double R2 = GTR_R2[igtr];
  constexpr double pi = TMath::Pi();
  double theta[]= {26.3/180.*pi,
		   50.6/180.*pi,
		   74.2/180.*pi,
		   97.9/180.*pi,
		   120.9/180.*pi };

  //auto *test = new TGeoTranslation(0.,0.,0.);
  //top->AddNode(gtr,200,test);
  for(int j = 0;j<2;j++){
    double s = (j==1)? 1.:-1.;
    for( int i = 0;i<4;i++){
      double R = (i%2==0)? R1 : R2;
      auto* t1 = new TGeoTranslation(R*sin(s*theta[i]),0.,R*cos(s*theta[i]));
      auto* r1 = new TGeoRotation();
      r1->RotateY(s*theta[i]*180./pi);
      auto* tt = new TGeoCombiTrans(*t1,*r1);
      top->AddNodeOverlap(gtr,i+100+j*10,tt);
      //cout << i << " " << theta[i] << endl;
   }
  }
  //TGeoTranslation test2(0.,0.,100.);
  //top->AddNodeOverlap(gtr2,201,&test2);
}

void E16Geo2DMan::build_hbd()
{
  E16MedMan* med_man = E16MedMan::get_instance();
  TGeoMedium* med = med_man->get("HBD");
  med->Print();

  double size = 60.*E16Geo::cm;
  double thick = 1.*E16Geo::cm;

  TGeoVolume* hbd = mgr->MakeBox("hbd",med,size/2.,size/2.,thick/2.);
  hbd->SetLineColor(kGreen);
  double R = 120.6*E16Geo::cm;
  constexpr double pi = TMath::Pi();
  double theta[]= {31./180.*pi,
		   62./180.*pi,
		   93./180.*pi,
		   124./180.*pi
  };
  //hbd->RandomPoints(10000);
  //auto *test = new TGeoTranslation(0.,0.,0.);
  //top->AddNode(gtr,200,test);
  for(int j = 0;j<2;j++){
    double s = (j==1)? 1.:-1.;
    for( int i = 0;i<4;i++){
      auto* t1 = new TGeoTranslation(R*sin(s*theta[i]),0.,R*cos(s*theta[i]));
      auto* r1 = new TGeoRotation();
      r1->RotateY(s*theta[i]*180./pi);
      auto* tt = new TGeoCombiTrans(*t1,*r1);
      top->AddNodeOverlap(hbd,i+200+j*10,tt);
      //cout << i << " " << theta[i] << endl;
   }
  }
  //TGeoTranslation test2(0.,0.,100.);
  //top->AddNodeOverlap(gtr2,201,&test2);
}

void E16Geo2DMan::build_lg()
{
  E16MedMan* med_man = E16MedMan::get_instance();
  TGeoMedium* med = med_man->get("LG");
  med->Print();
  
  double size = 70.*E16Geo::cm;
  double thick = 10.*E16Geo::cm;
  
  double h = 135*mm/2.;
  double ss = 142.3*mm/2.;
  double ls = 154.1*mm/2.;
  //double lg_block[]={0,-h,ls,-h,ss,h,0.,h};
  double lg_block[]={0,-h,ls,-h,ss,h,0.,h,
		     0,-h,ls,-h,ss,h,0.,h};
  TGeoArb8 *arb = new TGeoArb8(122.*mm/2.,lg_block);
  TGeoVolume* vol_lgblock = new TGeoVolume("lg_block",arb,med);
  TGeoVolume* lg = mgr->MakeBox("lg",med,size/2.,size/2.,thick/2.);
  lg->SetLineColor(kMagenta);
  lg->SetTransparency(50);
  double R = 140.*E16Geo::cm;
  constexpr double pi = TMath::Pi();
  double theta[]= {30./180.*pi,
		   60./180.*pi,
		   90./180.*pi,
		   120./180.*pi
  };

  //auto *test = new TGeoTranslation(0.,0.,0.);
  //top->AddNode(gtr,200,test);
  for(int j = 0;j<2;j++){
    double s = (j==1)? 1.:-1.;
    for( int i = 0;i<4;i++){
      auto* t1 = new TGeoTranslation(R*sin(s*theta[i]),0.,R*cos(s*theta[i]));
      auto* r1 = new TGeoRotation();
      r1->RotateY(s*theta[i]*180./pi);
      auto* tt = new TGeoCombiTrans(*t1,*r1);
      top->AddNodeOverlap(lg,i+300+j*10,tt);
      top->AddNodeOverlap(vol_lgblock,i+300+j*10,tt);
      //cout << i << " " << theta[i] << endl;
    }
  }
  //TGeoTranslation test2(0.,0.,100.);
  //top->AddNodeOverlap(gtr2,201,&test2);
}

*/
/*

void E16Geo2DMan::build_lg_using_geom()
{
  E16MedMan* med_man = E16MedMan::get_instance();
  TGeoMedium* med = med_man->get("LG");
  med->Print();
  
  int mod_id = 109;
  int block = 0;
  TVector3 lpos; lpos.SetXYZ(0,0,0);

  std::bitset<64> valid_lg_block;
  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      valid_lg_block.set(i*10+j);
    }
  }
  valid_lg_block.set(26);   valid_lg_block.set(36);

  const double height = 16.*E16Geo::cm;
  const double width = 13.*E16Geo::cm;
  const double depth = 12.*E16Geo::cm;
  constexpr double pi = TMath::Pi();

  for( int imodule = 101; imodule <= 109; imodule++){
    if (imodule == 105 ) continue;
    for(int block = 0;block<56;block++){
      if ( valid_lg_block[block] ) {
	int mod = E16Geo::ModuleID_2020to2013_27(imodule);
	TVector3 gpos = geom_v2->LG(mod,block)->GetGPos(lpos);
	TGeoVolume* lg = mgr->MakeBox("lg",med,width/2.,height/2.,depth/2.);
	lg->SetLineColor(kBlue);
	auto* t1 = new TGeoTranslation(gpos.X(),gpos.Y(),gpos.Z());
	auto* r1 = new TGeoRotation;
	double theta = gpos.Theta();
	if ( gpos.X() < 0 ) theta = -theta;
	r1->RotateY(theta*180./pi);
	auto* tt = new TGeoCombiTrans(*t1,*r1);
	top->AddNodeOverlap(lg,300+block+imodule*10,tt);
      }
    }
  }
}

*/

/*
void E16Geo2DMan::draw_tracks() {
  mgr->DrawTracks();
}

*/
