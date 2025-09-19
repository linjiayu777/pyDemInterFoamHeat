#include "demInterFoamBaseheat.H"


demInterFoamBaseheat::demInterFoamBaseheat() {
  int argc=1;
  const char *argv[]={"demInterFoamheat", NULL};
  char ** argv2 = const_cast<char **>(argv);
  args_ = new Foam::argList(argc, argv2);
  if (!args_->checkRootCase()) Foam::FatalError.exit();
  
  runTime_ = new Foam::Time(Foam::Time::controlDictName, *args_);

  mesh_ = new Foam::fvMesh(Foam::IOobject (Foam::fvMesh::defaultRegion,
                                           runTime_->timeName(),
                                           *runTime_,
                                           Foam::IOobject::MUST_READ));
// pimple_ = new pimpleControl(*mesh_);//added here why don't add here?

 
  
  p_rgh_ = new volScalarField (IOobject ("p_rgh", runTime_->timeName(), *mesh_,    //
                                     IOobject::MUST_READ,                         //
                                     IOobject::AUTO_WRITE), *mesh_);              //
  U_ = new volVectorField (IOobject ("U", runTime_->timeName(), *mesh_,
                                     IOobject::MUST_READ,
                                     IOobject::AUTO_WRITE), *mesh_);
  phi_ = new surfaceScalarField (IOobject("phi",runTime_->timeName(), *mesh_,      //
                                    IOobject::READ_IF_PRESENT,                   //   
                                    IOobject::AUTO_WRITE),                        // 
                                    linearInterpolate(*U_) & mesh_->Sf());        //
  //---------------------------------------------------------------------------------------
  T_ = new volScalarField (IOobject ("T", runTime_->timeName(), *mesh_, //added here
                                     IOobject::MUST_READ,
                                     IOobject::AUTO_WRITE), *mesh_);     //end here
  //---------------------------------------------------------------------------------------
  immiscibleIncompressibleTwoPhaseMixture mixture(*U_,*phi_);//added here
  volScalarField& alpha1(mixture.alpha1());//
  volScalarField& alpha2(mixture.alpha2());//
  const dimensionedScalar& rho1 = mixture.rho1();//
  const dimensionedScalar& rho2 = mixture.rho2();//end here


  IOdictionary transportProperties(IOobject("transportProperties",
                                            runTime_->constant(),
                                            *mesh_,
                                            IOobject::MUST_READ_IF_MODIFIED,
                                            IOobject::NO_WRITE));
  //----------------------------------------------------------------------------------------
  k1_ = new dimensionedScalar("k1",
                              transportProperties.lookup("k1"));
  Cv1_ = new dimensionedScalar("Cv1",
                              transportProperties.lookup("Cv1"));
  k2_ = new dimensionedScalar("k2",
                              transportProperties.lookup("k2"));
  Cv2_ = new dimensionedScalar("Cv2",
                              transportProperties.lookup("Cv2"));

  //----------------------------------------------------------------------------------------
  nu1_ = new dimensionedScalar("nu1",                               //
                              transportProperties.lookup("nu1"));                 //
  nu2_ = new dimensionedScalar("nu2",                                //
                              transportProperties.lookup("nu2"));                 //
              
  nu_  = new volScalarField(IOobject("nu", runTime_->timeName(), *mesh_,          //
                                    IOobject::READ_IF_PRESENT),                   //
                                    alpha1*(*nu1_) + alpha2*(*nu2_),              //
                                    alpha1.boundaryField().types());              //
  //----------------------------------------------------------------------------------------
  Deff_ = new volScalarField(IOobject("Deff", runTime_->timeName(), *mesh_,          //
                                    IOobject::READ_IF_PRESENT),                   //
                                    alpha1*(*k1_)/(*Cv1_)+(1.0-alpha1)*(*k2_)/(*Cv2_),             //
                                    alpha1.boundaryField().types());
  //----------------------------------------------------------------------------------------
  rho_ = new volScalarField(IOobject("rho", runTime_->timeName(), *mesh_,         //
                                    IOobject::READ_IF_PRESENT),                   //
                                    alpha1*rho1 + alpha2*rho2,                    //*
                                    alpha1.boundaryField().types());rho_->oldTime();//*old time?
  rhoPhi_ = new surfaceScalarField(IOobject("rhoPhi", runTime_->timeName(), *mesh_,//
                                    IOobject::NO_READ,                            //
                                    IOobject::NO_WRITE),                          //
                                    fvc::interpolate(*rho_)*(*phi_));                   //
  
  n_ = new volScalarField(IOobject ("n", runTime_->timeName(), *mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE), *mesh_);
  f_ = new volVectorField (IOobject ("f", runTime_->timeName(), *mesh_,
                                     IOobject::MUST_READ,
                                     IOobject::AUTO_WRITE), *mesh_);

  autoPtr<incompressible::turbulenceModel> turbulence
  (
    incompressible::turbulenceModel::New(*U_, *phi_, mixture)
  );

  //#include "readGravitationalAcceleration.H"//
  //#include "readhRef.H"                     // 
  //#include "gh.H"                           //added here?
  
 //-------------------------------------------------------------------------include readGravitational

   g_ = new uniformDimensionedVectorField (IOobject ("g", runTime_->constant(),*mesh_,
                                         IOobject::MUST_READ,
                                         IOobject::AUTO_WRITE));
 //-------------------------------------------------------------------------include readRef.H

   hRef_ = new uniformDimensionedScalarField (IOobject("hRef",runTime_->constant(),*mesh_,
                                                   IOobject::READ_IF_PRESENT,
                                                   IOobject::NO_WRITE),dimensionedScalar("hRef", dimLength, 0));
                                                                  
//--------------------------------------------------------------------------include gh.H
   ghRef_ = new dimensionedScalar (mag(g_->value()) > SMALL 
                                  ? *g_ & (cmptMag(g_->value())/mag(g_->value()))*(*hRef_)
                                  : dimensionedScalar("ghRef", g_->dimensions()*dimLength, 0));
   gh_ = new volScalarField (("gh"), (*g_ & (mesh_->C())) - (*ghRef_));
   ghf_ = new surfaceScalarField (("ghf"), (*g_ & (mesh_->Cf())) - (*ghRef_));

//--------------------------------------------------------------------------


   p_ = new volScalarField(IOobject ("p", runTime_->timeName(), *mesh_,
                                   IOobject::NO_READ,
                                   IOobject::AUTO_WRITE),
                                   *p_rgh_ + (*rho_)*(*gh_));

  pRefCell_ = 0;                                                                             //
  pRefValue_ = 0.0;                                                                          //
  setRefCell(*p_, *p_rgh_ ,mesh_->solutionDict().subDict("pimple"), pRefCell_, pRefValue_);  //pimple dict
  
  if(p_rgh_->needReference()) { *p_ += dimensionedScalar("p",p_->dimensions(),                    //*****
                                                  pRefValue_ - getRefCellValue(*p_, pRefCell_));
                                                  *p_rgh_ = *p_ -(*rho_)*(*gh_);}
                                                        
  mesh_->setFluxRequired(p_rgh_->name());                                                    // 
  mesh_->setFluxRequired(alpha1.name());                                                   //
  alphaPhi_ = new surfaceScalarField(IOobject("alphaPhi", runTime_->timeName(), *mesh_,      //
                                    IOobject::READ_IF_PRESENT,                               //
                                    IOobject::AUTO_WRITE),                                   //
                                    (*phi_)*fvc::interpolate(alpha1));                           //
  
  gradp_ = new volVectorField(fvc::grad(*p_));                                                //
  
  cumulativeContErr_ = 0.0;                                                                   //
  tmp<surfaceScalarField> talphaPhiCorr0;                                                   //
}

demInterFoamBaseheat::~demInterFoamBaseheat() {
  Info << "cleaning up demInterFoamBaseheat" << endl;
  if (gradp_) delete gradp_;
  if (rhoPhi_) delete rhoPhi_;
  if (f_) delete f_;   //
  if (n_) delete n_;         
  if (U_) delete U_;
  if (p_) delete p_;
  if (rho_) delete rho_;
  if (phi_) delete phi_;          //
//-----------------------------------------
  if (k1_) delete k1_;
  if (k2_) delete k2_;
  if (Cv1_) delete Cv1_;
  if (Cv2_) delete Cv2_;
  if (Deff_) delete Deff_;
  if (T_) delete T_;
//-----------------------------------------
  if (nu1_) delete nu1_;          //
  if (nu2_) delete nu2_;          //
  if (nu_) delete nu_;          //modifying
  if (alphaPhi_) delete alphaPhi_;//
  if (p_rgh_) delete p_rgh_;           //
  if (mesh_) delete mesh_;
  if (runTime_) delete runTime_;
  if (args_) delete args_;
}

double demInterFoamBaseheat::flux_on_patch(char *patch_name)
{
  label inletPatchi = (*mesh_).boundaryMesh().findPatchID(patch_name);
  if (inletPatchi == -1)
    throw std::runtime_error("Cannot find boundary patch");
  scalar massFlux = sum((*rhoPhi_).boundaryField()[inletPatchi]);       //phi->rhophi
  return massFlux;
}

double demInterFoamBaseheat::cell_flux(int cell, int face) {
  label gface = mesh_->cells()[cell][face];
  if (mesh_->faceOwner()[gface]==cell) return (*rhoPhi_)[gface];         //phi->rhophi
  else return -(*rhoPhi_)[gface];
}

int demInterFoamBaseheat::cell_near(double x, double y, double z) {
  meshSearch search(*mesh_);
  return search.findNearestCell(point(x,y,z));
}
