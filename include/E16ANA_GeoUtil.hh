#ifndef __E16ANA_GeoUtil_hh__
#define __E16ANA_GeoUtil_hh__

namespace E16Geo2D{
  
  inline int ModuleID_2013to2020_33(int id )
  {
    const int moduleID_2013to2020_33[33] = {
      10, 110, 210, 9, 109, 209, 8, 108, 208, 7, 107, 207, 6, 106, 206, 
      5, 105, 205,
      4,  104, 204, 3, 103, 203, 2, 102, 202, 1, 101, 201, 0, 100, 200 };
    return  moduleID_2013to2020_33[ id ];
  }
  inline int ModuleID_2020to2013_33(int id ){//GTR,SSD
    const int module_id_kawama33[3][11] = {
      {30, 27, 24, 21, 18, 15, 12, 9, 6, 3, 0},
    {31, 28, 25, 22, 19, 16, 13, 10, 7, 4, 1},
      {32, 29, 26, 23, 20, 17, 14, 11, 8, 5, 2}};
  return module_id_kawama33[ id / 100][id % 100];
  }
  inline int ModuleID_2013to2020_27( int id ){// //LG,HBD
    const int moduleID_2013to2020_27[27] = {
      9, 109, 209, 8, 108, 208, 7, 107, 207,  6, 106, 206, 
      5, 105, 205, 
      4, 104, 204, 3, 103, 203, 2, 102, 202,  1, 101, 201};
    return  moduleID_2013to2020_27[ id ];
  }
  inline int ModuleID_2020to2013_27(int id ){//LG,HBD
    const int module_id_kawama27[3][9] = {
    {24, 21, 18, 15, 12, 9, 6, 3, 0 },
    {25, 22, 19, 16, 13, 10, 7, 4, 1},
    {26, 23, 20, 17, 14, 11, 8, 5, 2}};
    
    int id2=(id %100) -1;
    
    return module_id_kawama27[ id / 100][ id2 ];
  }
  
};

#endif

