/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interFoam

Group
    grpMultiphaseSolvers

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

\*---------------------------------------------------------------------------*/


//
//// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
//int main(int argc, char *argv[])
//{
//    #include "setRootCase.H"
//    #include "createTime.H"
//    #include "createMesh.H"
//
//    pimpleControl pimple(mesh);
//
//    #include "createTimeControls.H"
//    #include "createRDeltaT.H"
//    #include "initContinuityErrs.H"
//    #include "createFields.H"
//    #include "createMRF.H"
//    #include "createFvOptions.H"
//    #include "correctPhi.H"//???//
//
//    turbulence->validate();
//
//    if (!LTS)
//    {
//        #include "readTimeControls.H"
//        #include "CourantNo.H"
//        #include "setInitialDeltaT.H"
//    }
//
//    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//    Info<< "\nStarting time loop\n" << endl;
        //#include "pimpleControl.H"
        
        

//demInterFoamBase::demInterFoamBase(){

//	Info << "creating demInterFoam object" << nl <<endl;//
//	f_ = new volVectorField(IOobject ("f", runTime_->timeName(), *mesh_,
//	                                  IOobject::MUST_READ,
//	                                  IOobject::AUTO_WRITE), *mesh_);
//      n_ = new volScalarField(IOobject ("n", runTime_->timeName(), *mesh_,
//	                                  IOobject::MUST_READ,
//	                                  IOobject::AUTO_WRITE), *mesh_);
//	pimple_ = new pimpleControl(*mesh_);
//}

//demInterFoamBase::~demInterFoamBase(){
//	Info << "cleaning up demInterFoam" << nl << endl;
//	if (pimple_) delete pimple_;
//	if (f_) delete f_;
//       if (n_) delete n_;
//}
#include "demInterFoamheat.H"
demInterFoamheat::demInterFoamheat(){
 Info << "creating demInterFoamheat object" << nl << endl;
   pimple_ = new pimpleControl(*mesh_);
 }

demInterFoamheat::~demInterFoamheat(){
 Info << "cleaning up demInterFoamheat object" << nl << endl;
  if(pimple_) delete pimple_;
 }

void demInterFoamheat::run(double time_increment) {
//

//
        Foam::Time           &runTime = *runTime_;
	Foam::fvMesh         &mesh = *mesh_;
        pimpleControl        &pimple =*pimple_;
        volScalarField    &p_rgh = *p_rgh_;
        volVectorField    &U = *U_;
        surfaceScalarField &phi = *phi_;
         //--------------------------------
        immiscibleIncompressibleTwoPhaseMixture mixture(U,phi);//added her
        volScalarField& alpha1(mixture.alpha1());//
        volScalarField& alpha2(mixture.alpha2());//
        const dimensionedScalar& rho1 = mixture.rho1();//
        const dimensionedScalar& rho2 = mixture.rho2();//end here


        //---------------------------------
        dimensionedScalar &nu1 = *nu1_;
        dimensionedScalar &nu2 = *nu2_;
        volScalarField    &nu = *nu_;
        //---------------------------------
        dimensionedScalar &k1 = *k1_;
        dimensionedScalar &k2 = *k2_;
        dimensionedScalar &Cv1 = *Cv1_;
        dimensionedScalar &Cv2 = *Cv2_;
        volScalarField    &Deff = *Deff_;
        volScalarField    &T = *T_;
        //---------------------------------
	volScalarField    &rho = *rho_;
        surfaceScalarField &rhoPhi = *rhoPhi_;
        volVectorField    &f = *f_;
        volScalarField    &n = *n_;
        autoPtr<incompressible::turbulenceModel> turbulence
        (
         incompressible::turbulenceModel::New(U, phi, mixture)
         );
         //--------------------------------
        uniformDimensionedVectorField &g = *g_;
        uniformDimensionedScalarField &hRef = *hRef_;
        dimensionedScalar  &ghRef = *ghRef_;
        volScalarField     &gh = *gh_;
        surfaceScalarField &ghf = *ghf_; 
         //--------------------------------

        volScalarField    &p = *p_;
        surfaceScalarField   &alphaPhi = *alphaPhi_;   
        volVectorField    &gradp = *gradp_;
        label &pRefCell = pRefCell_;
        scalar &pRefValue = pRefValue_;
       
        tmp<surfaceScalarField> talphaPhiCorr0;
	//scalar &cumulativeContErr = cumulativeContEEr_;

      //add 
      // #include "createTimeControls.H"
       #include "createRDeltaT.H"
       #include "initContinuityErrs.H"
       #include "createMRF.H"
       #include "createFvOptions.H"
       #include "correctPhi.H"
       #include "turbulentTransportModel.H"
       turbulence->validate();
       

    runTime.setEndTime(runTime.value() + time_increment);

    while (runTime.loop())//runTime.loop())+runTime++= runTime.run
    {
        #include "CourantNo.H"//need remake
        #include "alphaCourantNo.H"
        //#include "setDeltaT.H"
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();
           //#include "TEqn"---------------------------------------------------------
           fvScalarMatrix TEqn
           (
             fvm::ddt(rho,T)
            +fvm::div(rhoPhi,T)
            -fvm::laplacian(Deff,T) 
           );
           TEqn.solve();
            //#include "TEqn"---------------------------------------------------------
           
           // #include "UEqn.H"------------------------------------------------------
            MRF.correctBoundaryVelocity(U);

           fvVectorMatrix UEqn
           (
              fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
            + MRF.DDt(rho, U)
            + turbulence->divDevRhoReff(rho, U)-rho*f/n//add body force into Ueqn
           ==
              fvOptions(rho, U)
            );

           UEqn.relax();

           fvOptions.constrain(UEqn);

           if (pimple.momentumPredictor())
           {
              solve
             (
                UEqn
              ==
                fvc::reconstruct
                (
                    (
                        mixture.surfaceTensionForce()
                      - ghf*fvc::snGrad(rho)
                      - fvc::snGrad(p_rgh)
                    ) * mesh.magSf()
                )
             );

             fvOptions.correct(U);
            }
            //#include "UEqn.H------------------------------------------------------
            // --- Pressure corrector loop
            while (pimple.correct())
            {
            //#include "pEqn.H"-----------------------------------------------------
              #include "pEqn.H"
            //#include "pEqn.H"-----------------------------------------------------
            }
            
           //#include "TEqn"---------------------------------------------------------
           //fvScalarMatrix TEqn
           //(
           //  fvm::ddt(rho,T)
           // +fvm::div(rhoPhi,T)
           // -fvm::laplacian(Deff,T) 
          // );
          // TEqn.solve();
            //#include "TEqn"---------------------------------------------------------

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    //return 0;
}

// ************************************************************************* //
