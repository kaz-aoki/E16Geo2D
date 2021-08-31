//2020-11-26, uploaded by yokkaich
//2016-05-02, uploaded by nakai
//2016-04-01, uploaded by nakai
//2015-11-14, uploaded by nakai
//2015-05-29, uploaded by yokkaich
//2015-03-27, uploaded by yokkaich
//2015-01-20, uploaded by yokkaich
//2015-01-20, uploaded by yokkaich
//2015-01-05, modified by yokkaich
//2015-01-05, modified by yokkaich
//2014-08-27, modified by kawama
//2014-05-07, modified by kawama
//2014-04-30, modified by kawama
//2014-04-24, modified by kawama
//2013-11-14, modified by kawama
//2013-05-13, modified by kawama

//E16ANA_Geometry.cc 150118 by S. Yokkaichi
//    Last modified at <>
//

#include <fstream> 
#include <iostream> 
#include <stdlib.h> 


#include <TVector3.h>
#include <TRotation.h>
#include <TRandom3.h>
#include <Math/Rotation3D.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>

#include "E16ANA_Geometry.hh"
#include "E16ANA_Transform.hh"

using namespace std;
using namespace E16Geo2D;
const double PI=3.14159265358979;
const int MAX_LENGTH=128;

/*** GTR & GEM ID 
// From the beam upstream
                             +y
  //- 0-||- 3-||- 6-||- 9-||-12-||-15-||-18-||-21-||-24-\\
+x||- 1-||- 4-||- 7-||-10-||-13-||-16-||-19-||-22-||-25-||-x
  \\- 2-||- 5-||- 8-||-11-||-14-||-17-||-20-||-23-||-26-//
                             -y


                                    +y
   //- 0-||- 3-||- 6-||- 9-||-12-||-15-||-18-||-21-||-24-||-27-||-30\\
+x ||- 1-||- 4-||- 7-||-10-||-13-||-16-||-19-||-22-||-25-||-28-||-31||  -x
   \\- 2-||- 5-||- 8-||-11-||-14-||-17-||-20-||-23-||-26-||-29-||-32//
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

 E16ANA_GeometryV1::E16ANA_GeometryV1(char* filename)
:filename_(filename),
   //GTRdx_smear(0.02), GTRdy_smear(0.02), GTRdz_smear(0.05),
   //GTRrotx_smear(0.004*180.0/PI), GTRroty_smear(0.004*180.0/PI), GTRrotz_smear(0.004*180.0/PI),
   GTRdx_frame_smear(0.02), GTRdy_frame_smear(0.02), GTRdz_frame_smear(0.05),
   GTRrotx_frame_smear(0.004*180.0/PI), GTRroty_frame_smear(0.004*180.0/PI), GTRrotz_frame_smear(0.004*180.0/PI)
{
  GTRdx_smear[0][0] = 0.005; GTRdx_smear[1][0] = 0.005; GTRdx_smear[2][0] = 0.005;
  GTRdx_smear[0][1] = 0.005; GTRdx_smear[1][1] = 0.005; GTRdx_smear[2][1] = 0.005;
  GTRdx_smear[0][2] = 0.005; GTRdx_smear[1][2] = 0.005; GTRdx_smear[2][2] = 0.005;

  GTRdy_smear[0][0] = 0.005; GTRdy_smear[1][0] = 0.005; GTRdy_smear[2][0] = 0.005;
  GTRdy_smear[0][1] = 0.005; GTRdy_smear[1][1] = 0.005; GTRdy_smear[2][1] = 0.005;
  GTRdy_smear[0][2] = 0.005; GTRdy_smear[1][2] = 0.005; GTRdy_smear[2][2] = 0.005;

  GTRdz_smear[0][0] = 0.01; GTRdz_smear[1][0] = 0.01; GTRdz_smear[2][0] = 0.01;
  GTRdz_smear[0][1] = 0.01; GTRdz_smear[1][1] = 0.01; GTRdz_smear[2][1] = 0.01;
  GTRdz_smear[0][2] = 0.01; GTRdz_smear[1][2] = 0.01; GTRdz_smear[2][2] = 0.01;

  GTRrotx_smear[0][0] = 0.00025*180.0/PI; GTRrotx_smear[1][0] = 0.00025*180.0/PI; GTRrotx_smear[2][0] = 0.00025*180.0/PI;
  GTRrotx_smear[0][1] = 0.00025*180.0/PI; GTRrotx_smear[1][1] = 0.00025*180.0/PI; GTRrotx_smear[2][1] = 0.00025*180.0/PI;
  GTRrotx_smear[0][2] = 0.00025*180.0/PI; GTRrotx_smear[1][2] = 0.00025*180.0/PI; GTRrotx_smear[2][2] = 0.00025*180.0/PI;

  GTRroty_smear[0][0] = 0.001*180.0/PI; GTRroty_smear[1][0] = 0.001*180.0/PI; GTRroty_smear[2][0] = 0.001*180.0/PI;
  GTRroty_smear[0][1] = 0.001*180.0/PI; GTRroty_smear[1][1] = 0.001*180.0/PI; GTRroty_smear[2][1] = 0.001*180.0/PI;
  GTRroty_smear[0][2] = 0.001*180.0/PI; GTRroty_smear[1][2] = 0.001*180.0/PI; GTRroty_smear[2][2] = 0.001*180.0/PI;

  GTRrotz_smear[0][0] = 0.0007*180.0/PI; GTRrotz_smear[1][0] = 0.0007*180.0/PI; GTRrotz_smear[2][0] = 0.0007*180.0/PI;
  GTRrotz_smear[0][1] = 0.0007*180.0/PI; GTRrotz_smear[1][1] = 0.0007*180.0/PI; GTRrotz_smear[2][1] = 0.0007*180.0/PI;
  GTRrotz_smear[0][2] = 0.0007*180.0/PI; GTRrotz_smear[1][2] = 0.0007*180.0/PI; GTRrotz_smear[2][2] = 0.0007*180.0/PI;

  ifstream ifs(filename);
  //  ifs.open(filename);
  //  if( ! ifs.is_open() ) {
  if( ifs.fail() ) {
    cerr<<__func__<<" ("<<__FILE__<<"/"<<__LINE__<<") exit: geom file cannot open...: "<<filename_<<endl;
        exit(1);
  }
  ReadParam();
  cerr<<"E16ANA_GeometryV1 initialized by "<<filename_<<endl;
}

 E16ANA_GeometryV1::~E16ANA_GeometryV1(){
}

ROOT::Math::Plane3D E16ANA_GeometryV1::GetPlane(int layer_id, int module_id, double local_z){
   ROOT::Math::Plane3D plane;
   ROOT::Math::Plane3D::Vector normal;
   ROOT::Math::Plane3D::Point origin;
   origin = this->GetGPos(TVector3(0.0,0.0,local_z),layer_id,module_id);
   normal.SetXYZ(0.0,0.0,1.0);
   normal = ROOT::Math::RotationX(this->GTRrotx[layer_id][module_id]/180.0*M_PI)*normal;
   normal = ROOT::Math::RotationY(this->GTRroty[layer_id][module_id]/180.0*M_PI)*normal;
   normal = ROOT::Math::RotationZ(this->GTRrotz[layer_id][module_id]/180.0*M_PI)*normal;
   normal = ROOT::Math::RotationX(this->GTRrotx_frame[module_id]/180.0*M_PI)*normal;
   normal = ROOT::Math::RotationY(this->GTRroty_frame[module_id]/180.0*M_PI)*normal;
   normal = ROOT::Math::RotationZ(this->GTRrotz_frame[module_id]/180.0*M_PI)*normal;
   normal = ROOT::Math::RotationY(this->GTRtheta[module_id]/180.0*M_PI)*normal;

#ifdef E16_ROOT6
   plane=ROOT::Math::Plane3D(normal,origin);
#else
   plane=ROOT::Math::Plane3D::Plane3D(normal,origin);
#endif


   return plane;
}

void E16ANA_GeometryV1::ReadParam(){
   ifstream ifs(filename_);
   char line[MAX_LENGTH];
   double val1, val2, val3;
   int igy=0,jgy=0;
   int igz=0,jgz=0;
   int iga=0;
   int ihy=0;
   int ihz=0;
   int igdx=0,jgdx=0;
   int igdy=0,jgdy=0;
   int igdz=0,jgdz=0;
   int igrotx=0,jgrotx=0;
   int igroty=0,jgroty=0;
   int igrotz=0,jgrotz=0;
   int igdx_frame=0,igdy_frame=0,igdz_frame=0;
   int igcx_frame=0,igcy_frame=0,igcz_frame=0;
   int igrotx_frame=0,igroty_frame=0,igrotz_frame=0;
   int iline=0;
   int iline2=0;
   while(ifs.getline(line,sizeof(line))){
      if(sscanf(line,"%lf %lf %lf", &val1, &val2, &val3)==3){
         if(iline<33){
            GTRy[jgy][igy]=val1;igy++;
            GTRy[jgy][igy]=val2;igy++;
            GTRy[jgy][igy]=val3;igy++;
            if(igy==33){
               jgy++;
               igy=0;
            }
         }
         else if (iline<33*2){
            GTRz[jgz][igz]=val1;igz++;
            GTRz[jgz][igz]=val2;igz++;
            GTRz[jgz][igz]=val3;igz++;
            if(igz==33){
               jgz++;
               igz=0;
            }
         }
         else if (iline<33*2+9){
            HBDy[ihy]=val1;ihy++;
            HBDy[ihy]=val2;ihy++;
            HBDy[ihy]=val3;ihy++;
         }
         else if (iline<33*2+18){
            HBDz[ihz]=val1;ihz++;
            HBDz[ihz]=val2;ihz++;
            HBDz[ihz]=val3;ihz++;
         }
         else if (iline<33*2+29){
            GTRtheta[iga]=val1;iga++;
            GTRtheta[iga]=val2;iga++;
            GTRtheta[iga]=val3;iga++;
         }
         //alignment data for GTR
         else if (iline<33*2+29+33){
            //cout << "GTRdx setting... : iline = " << iline << endl;
            //cout << "    (val1,val2,val3) = (" << val1 << ", " << val2 << ", " << val3 << ")" <<endl;
            GTRdx[jgdx][igdx]=val1;igdx++;
            GTRdx[jgdx][igdx]=val2;igdx++;
            GTRdx[jgdx][igdx]=val3;igdx++;
            if(igdx==33){
               jgdx++;
               igdx=0;
            }
         }
         else if (iline<33*2+29+33*2){
            //cout << "GTRdy setting... : iline = " << iline << endl;
            //cout << "    (val1,val2,val3) = (" << val1 << ", " << val2 << ", " << val3 << ")" <<endl;
            GTRdy[jgdy][igdy]=val1;igdy++;
            GTRdy[jgdy][igdy]=val2;igdy++;
            GTRdy[jgdy][igdy]=val3;igdy++;
            if(igdy==33){
               jgdy++;
               igdy=0;
            }
         }
         else if (iline<33*2+29+33*3){
            GTRdz[jgdz][igdz]=val1;igdz++;
            GTRdz[jgdz][igdz]=val2;igdz++;
            GTRdz[jgdz][igdz]=val3;igdz++;
            if(igdz==33){
               jgdz++;
               igdz=0;
            }
         }
         else if (iline<33*2+29+33*4){
            GTRrotx[jgrotx][igrotx]=val1;igrotx++;
            GTRrotx[jgrotx][igrotx]=val2;igrotx++;
            GTRrotx[jgrotx][igrotx]=val3;igrotx++;
            if(igrotx==33){
               jgrotx++;
               igrotx=0;
            }
         }
         else if (iline<33*2+29+33*5){
            GTRroty[jgroty][igroty]=val1;igroty++;
            GTRroty[jgroty][igroty]=val2;igroty++;
            GTRroty[jgroty][igroty]=val3;igroty++;
            if(igroty==33){
               jgroty++;
               igroty=0;
            }
         }
         else if (iline<33*2+29+33*6){
         //else {
            GTRrotz[jgrotz][igrotz]=val1;igrotz++;
            GTRrotz[jgrotz][igrotz]=val2;igrotz++;
            GTRrotz[jgrotz][igrotz]=val3;igrotz++;
            if(igrotz==33){
               jgrotz++;
               igrotz=0;
            }
         }
         iline++;
      }else if(sscanf(line,"%lf", &val1)==1){
         if (iline2<11){
            GTRdx_frame[igdx_frame]=val1;igdx_frame++;
            GTRdx_frame[igdx_frame]=val1;igdx_frame++;
            GTRdx_frame[igdx_frame]=val1;igdx_frame++;
         }
         else if (iline2<11*2){
            GTRdy_frame[igdy_frame]=val1;igdy_frame++;
            GTRdy_frame[igdy_frame]=val1;igdy_frame++;
            GTRdy_frame[igdy_frame]=val1;igdy_frame++;
         }
         else if (iline2<11*3){
            GTRdz_frame[igdz_frame]=val1;igdz_frame++;
            GTRdz_frame[igdz_frame]=val1;igdz_frame++;
            GTRdz_frame[igdz_frame]=val1;igdz_frame++;
         }
         else if (iline2<11*4){
            GTRcx_frame[igcx_frame]=val1;igcx_frame++;
            GTRcx_frame[igcx_frame]=val1;igcx_frame++;
            GTRcx_frame[igcx_frame]=val1;igcx_frame++;
         }
         else if (iline2<11*5){
            GTRcy_frame[igcy_frame]=val1;igcy_frame++;
            GTRcy_frame[igcy_frame]=val1;igcy_frame++;
            GTRcy_frame[igcy_frame]=val1;igcy_frame++;
         }
         else if (iline2<11*6){
            GTRcz_frame[igcz_frame]=val1;igcz_frame++;
            GTRcz_frame[igcz_frame]=val1;igcz_frame++;
            GTRcz_frame[igcz_frame]=val1;igcz_frame++;
         }
         else if (iline2<11*7){
            GTRrotx_frame[igrotx_frame]=val1;igrotx_frame++;
            GTRrotx_frame[igrotx_frame]=val1;igrotx_frame++;
            GTRrotx_frame[igrotx_frame]=val1;igrotx_frame++;
         }
         else if (iline2<11*8){
            GTRroty_frame[igroty_frame]=val1;igroty_frame++;
            GTRroty_frame[igroty_frame]=val1;igroty_frame++;
            GTRroty_frame[igroty_frame]=val1;igroty_frame++;
         }
         else if (iline2<11*9){
            GTRrotz_frame[igrotz_frame]=val1;igrotz_frame++;
            GTRrotz_frame[igrotz_frame]=val1;igrotz_frame++;
            GTRrotz_frame[igrotz_frame]=val1;igrotz_frame++;
         }
         iline2++;
      }
   }
}

TVector3 E16ANA_GeometryV1::GetDetectorCenter(int layer_id, int module_id){
   TVector3 cent(0.0,0.0,0.0); // local (0,0,0) corresponds global coordinate of detector center.
   return E16ANA_Transform::GetGPos(cent, this, layer_id, module_id);
}

void E16ANA_GeometryV1::WriteGeometryFile(char *filename){
   ofstream ofs(filename);
   ofs << "#GTR 10cm center-y" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRy[0][module_id]   << " "
          << GTRy[0][module_id+1] << " "
          << GTRy[0][module_id+2] << endl;
   }
   ofs << "#GTR 20cm center-y" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRy[1][module_id]   << " "
          << GTRy[1][module_id+1] << " "
          << GTRy[1][module_id+2] << endl;
   }
   ofs << "#GTR 30cm center-y" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRy[2][module_id]   << " "
          << GTRy[2][module_id+1] << " "
          << GTRy[2][module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#GTR 10cm center-z" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRz[0][module_id]   << " "
          << GTRz[0][module_id+1] << " "
          << GTRz[0][module_id+2] << endl;
   }
   ofs << "#GTR 20cm center-z" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRz[1][module_id]   << " "
          << GTRz[1][module_id+1] << " "
          << GTRz[1][module_id+2] << endl;
   }
   ofs << "#GTR 30cm center-z" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRz[2][module_id]   << " "
          << GTRz[2][module_id+1] << " "
          << GTRz[2][module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#HBD center-y" << endl;
   for(int module_id=0; module_id<27; module_id+=3){
      ofs << HBDy[module_id]   << " "
          << HBDy[module_id+1] << " "
          << HBDy[module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#HBD center-z" << endl;
   for(int module_id=0; module_id<27; module_id+=3){
      ofs << HBDz[module_id]   << " "
          << HBDz[module_id+1] << " "
          << HBDz[module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#GTR angle" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRtheta[module_id]   << " "
          << GTRtheta[module_id+1] << " "
          << GTRtheta[module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#GTR1 Local dx" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdx[0][module_id]   << " "
          << GTRdx[0][module_id+1] << " "
          << GTRdx[0][module_id+2] << endl;
   }
   ofs << "#GTR2 Local dx" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdx[1][module_id]   << " "
          << GTRdx[1][module_id+1] << " "
          << GTRdx[1][module_id+2] << endl;
   }
   ofs << "#GTR3 Local dx" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdx[2][module_id]   << " "
          << GTRdx[2][module_id+1] << " "
          << GTRdx[2][module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#GTR1 Local dy" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdy[0][module_id]   << " "
          << GTRdy[0][module_id+1] << " "
          << GTRdy[0][module_id+2] << endl;
   }
   ofs << "#GTR2 Local dy" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdy[1][module_id]   << " "
          << GTRdy[1][module_id+1] << " "
          << GTRdy[1][module_id+2] << endl;
   }
   ofs << "#GTR3 Local dy" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdy[2][module_id]   << " "
          << GTRdy[2][module_id+1] << " "
          << GTRdy[2][module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#GTR1 Local dz" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdz[0][module_id]   << " "
          << GTRdz[0][module_id+1] << " "
          << GTRdz[0][module_id+2] << endl;
   }
   ofs << "#GTR2 Local dz" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdz[1][module_id]   << " "
          << GTRdz[1][module_id+1] << " "
          << GTRdz[1][module_id+2] << endl;
   }
   ofs << "#GTR3 Local dz" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdz[2][module_id]   << " "
          << GTRdz[2][module_id+1] << " "
          << GTRdz[2][module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#GTR1 Local rotx" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRrotx[0][module_id]   << " "
          << GTRrotx[0][module_id+1] << " "
          << GTRrotx[0][module_id+2] << endl;
   }
   ofs << "#GTR2 Local rotx" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRrotx[1][module_id]   << " "
          << GTRrotx[1][module_id+1] << " "
          << GTRrotx[1][module_id+2] << endl;
   }
   ofs << "#GTR3 Local rotx" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRrotx[2][module_id]   << " "
          << GTRrotx[2][module_id+1] << " "
          << GTRrotx[2][module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#GTR1 Local roty" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRroty[0][module_id]   << " "
          << GTRroty[0][module_id+1] << " "
          << GTRroty[0][module_id+2] << endl;
   }
   ofs << "#GTR2 Local roty" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRroty[1][module_id]   << " "
          << GTRroty[1][module_id+1] << " "
          << GTRroty[1][module_id+2] << endl;
   }
   ofs << "#GTR3 Local roty" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRroty[2][module_id]   << " "
          << GTRroty[2][module_id+1] << " "
          << GTRroty[2][module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#GTR1 Local rotz" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRrotz[0][module_id]   << " "
          << GTRrotz[0][module_id+1] << " "
          << GTRrotz[0][module_id+2] << endl;
   }
   ofs << "#GTR2 Local rotz" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRrotz[1][module_id]   << " "
          << GTRrotz[1][module_id+1] << " "
          << GTRrotz[1][module_id+2] << endl;
   }
   ofs << "#GTR3 Local rotz" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRrotz[2][module_id]   << " "
          << GTRrotz[2][module_id+1] << " "
          << GTRrotz[2][module_id+2] << endl;
   }
   ofs << endl;

   ofs << "#GTR frame dx" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdx_frame[module_id]   << endl;
   }
   ofs << endl;

   ofs << "#GTR frame dy" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdy_frame[module_id]   << endl;
   }
   ofs << endl;

   ofs << "#GTR frame dz" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRdz_frame[module_id]   << endl;
   }
   ofs << endl;

   ofs << "#GTR frame rotation center x" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRcx_frame[module_id]   << endl;
   }
   ofs << endl;

   ofs << "#GTR frame rotation center y" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRcy_frame[module_id]   << endl;
   }
   ofs << endl;

   ofs << "#GTR frame rotation center z" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRcz_frame[module_id]   << endl;
   }
   ofs << endl;

   ofs << "#GTR frame rotx" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRrotx_frame[module_id]   << endl;
   }
   ofs << endl;

   ofs << "#GTR frame roty" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRroty_frame[module_id]   << endl;
   }
   ofs << endl;

   ofs << "#GTR frame rotz" << endl;
   for(int module_id=0; module_id<33; module_id+=3){
      ofs << GTRrotz_frame[module_id]   << endl;
   }
   ofs << endl;
}

void E16ANA_GeometryV1::SetRandomGeometryError(int seed){
   TRandom3 rand(seed);

   // Local 6 parameters
   for(int layer_id=0; layer_id<3; layer_id++){
      for(int module_id=0; module_id<33; module_id++){
         GTRdx[layer_id][module_id] = rand.Gaus(0.0,GTRdx_smear[layer_id][module_id%3]);
         GTRdy[layer_id][module_id] = rand.Gaus(0.0,GTRdy_smear[layer_id][module_id%3]);
         GTRdz[layer_id][module_id] = rand.Gaus(0.0,GTRdz_smear[layer_id][module_id%3]);
         GTRrotx[layer_id][module_id] = rand.Gaus(0.0,GTRrotx_smear[layer_id][module_id%3]);
         GTRroty[layer_id][module_id] = rand.Gaus(0.0,GTRroty_smear[layer_id][module_id%3]);
         GTRrotz[layer_id][module_id] = rand.Gaus(0.0,GTRrotz_smear[layer_id][module_id%3]);
      }
   }

   // Frame 6 parameters
   for(int module_id=0; module_id<33; module_id+=3){
      GTRdx_frame[module_id] = rand.Gaus(0.0,GTRdx_frame_smear);
      GTRdy_frame[module_id] = rand.Gaus(0.0,GTRdy_frame_smear);
      GTRdz_frame[module_id] = rand.Gaus(0.0,GTRdz_frame_smear);
      GTRrotx_frame[module_id] = rand.Gaus(0.0,GTRrotx_frame_smear);
      GTRroty_frame[module_id] = rand.Gaus(0.0,GTRroty_frame_smear);
      GTRrotz_frame[module_id] = rand.Gaus(0.0,GTRrotz_frame_smear);
      GTRdx_frame[module_id+1] = GTRdx_frame[module_id]; GTRdx_frame[module_id+2] = GTRdx_frame[module_id];
      GTRdy_frame[module_id+1] = GTRdy_frame[module_id]; GTRdy_frame[module_id+2] = GTRdy_frame[module_id];
      GTRdz_frame[module_id+1] = GTRdz_frame[module_id]; GTRdz_frame[module_id+2] = GTRdz_frame[module_id];
      GTRrotx_frame[module_id+1] = GTRrotx_frame[module_id]; GTRrotx_frame[module_id+2] = GTRrotx_frame[module_id];
      GTRroty_frame[module_id+1] = GTRroty_frame[module_id]; GTRroty_frame[module_id+2] = GTRroty_frame[module_id];
      GTRrotz_frame[module_id+1] = GTRrotz_frame[module_id]; GTRrotz_frame[module_id+2] = GTRrotz_frame[module_id];
   }
}

