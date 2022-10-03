#ifdef __CLING__
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ namespace supera+;
#pragma link C++ namespace supera::Point3D+;
#pragma link C++ namespace supera::EDep+;
#pragma link C++ namespace supera::Vertex+;
#pragma link C++ namespace std::vector<supera::Point3D>+;
#pragma link C++ namespace std::vector<supera::EDep>+;
#pragma link C++ namespace std::vector<supera::Vertex>+;
#pragma link C++ function supera::SamplePointsFromLine(const Point3D&, const Point3D&, const double&);
#pragma link C++ namespace edep2supera+;
#pragma link C++ class edep2supera::SuperaDriver+;
#endif