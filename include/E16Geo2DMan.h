#ifndef __E16Geo2DMan_h_
#define __E16Geo2DMan_h_

class TH2D;

namespace E16Geo2D {
  class E16ANA_GeometryV2;
  const double mm = 1.;
  const double cm = 10.;
};

class E16Geo2DMan
{
 public:
  E16Geo2DMan();
  void draw(const char* option=""); // ogl= OpenGL

  void draw_gtr();
  void draw_hbd();
  
  /*
  void testmethod() { ; }
  
  void build_world();
  void build_gtr(int igtr);
  void build_hbd();

  void build_lg();
  void build_lg_using_geom();

  void draw_tracks();

  E16ANA_GeometryV2* get_geom();
  */
  ~E16Geo2DMan();
 public:
  E16Geo2D::E16ANA_GeometryV2* geom_v2;
  TH2D* frame;
};

#endif

