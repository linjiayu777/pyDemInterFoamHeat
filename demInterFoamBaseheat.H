#pragma once
//#include "fvCFD.H"
#include "pisoControl.H"
#include "meshSearch.H"
//#include "pimpleControl.H"// is there pimpleControl.H?


#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"


#include <iostream>
#include <fstream>
#include <stdexcept>

#define VEC_GET(T) double T(int i, int j) { return (*T##_)[i][j]; }
#define VEC_SET(T) void set_##T(int i, int j, double v) { (*T##_)[i][j]=v; }
#define VEC_SETGET(T) VEC_SET(T) VEC_GET(T)

#define GET(T) double T(int i) { return (*T##_)[i]; }
#define SET(T) void set_##T(int i, double v) { (*T##_)[i]=v; }
#define SETGET(T) SET(T) GET(T)

class demInterFoamBaseheat {
public:
  demInterFoamBaseheat();
  virtual ~demInterFoamBaseheat();
  int nCells() { return mesh_->nCells(); }
  double nu1() { return nu1_->value(); }
  double nu2() { return nu2_->value(); }
//----------------------------------------------------------
  double k1() { return k1_->value(); }
  double k2() { return k2_->value(); }
  double Cv1() { return Cv1_->value(); }
  double Cv2() { return Cv2_->value(); }
//----------------------------------------------------------
  int nNodes() { return mesh_->nPoints(); }
  int nFaces() { return mesh_->nFaces(); }
  int face_node(int face, int node) { return mesh_->faces()[face][node]; }
  int cell_face(int cell, int face) { return mesh_->cells()[cell][face]; }
  double cell_flux(int cell, int face);
  int cell_near(double x, double y, double z);
  double flux_on_patch(char *patch_name);
  double cell_center(int cell, int j) { return mesh_->cellCentres()[cell][j]; }
  double cell_volume(int cell) { return mesh_->cellVolumes()[cell]; }
  double face_center(int face, int j) { return mesh_->faceCentres()[face][j]; }
  double node_pos(int i, int j) { return mesh_->points()[i][j]; }
  int element(int i, int j) { return mesh_->cellShapes()[i][j]; }
  void set_dt(double t) { runTime_->setDeltaT(t); }
  double dt() { return runTime_->deltaT().value(); }
  void set_endtime(double t) { runTime_->setEndTime(t); }
  virtual void run(double v)=0;
  //void run(double t) { runTime_->setDeltaT(t); }
   
  VEC_SETGET(f)
  SETGET(n)
  VEC_SETGET(U)
  SETGET(p)
  VEC_GET(gradp)
  GET(phi) 
  SETGET(p_rgh)
  SETGET(rho)
  GET(nu)
  GET(Deff)
  //--------------
  SETGET(T)
  //--------------
  SETGET(rhoPhi)
  GET(alphaPhi)


   
  Foam::argList *args_;
  Foam::Time  *runTime_;
  Foam::fvMesh *mesh_;
 // pimpleControl *pimple_;
  
   
  
  dimensionedScalar *ghRef_;
  volScalarField *gh_;
  surfaceScalarField *ghf_; 
  uniformDimensionedVectorField *g_;
  uniformDimensionedScalarField *hRef_;

 //---------------------
  dimensionedScalar *nu1_;
  dimensionedScalar *nu2_;
 //-------------------------------
  dimensionedScalar *k1_;
  dimensionedScalar *k2_;
  dimensionedScalar *Cv1_;
  dimensionedScalar *Cv2_;
 //-------------------------------
  volScalarField *p_rgh_;
  volVectorField *U_;
  volVectorField *f_;
  surfaceScalarField *phi_;
  volScalarField *nu_;
 //-------------------------------
  volScalarField *Deff_;
  volScalarField *T_;
 //-------------------------------
  volScalarField *rho_;
  surfaceScalarField *rhoPhi_;
  volScalarField *n_;
  volScalarField *p_;
  label pRefCell_;
  scalar pRefValue_;
  surfaceScalarField *alphaPhi_; 
  volVectorField *gradp_;
  scalar cumulativeContErr_;
};
