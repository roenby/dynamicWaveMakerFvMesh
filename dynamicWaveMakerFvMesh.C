/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "dynamicWaveMakerFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "interpolateSplineXY.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicWaveMakerFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicWaveMakerFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicWaveMakerFvMesh::dynamicWaveMakerFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    ),
    times_(dynamicMeshCoeffs_.lookup("times")),
    pistonPositions_(dynamicMeshCoeffs_.lookup("pistonPositions")),
    nPistons_(pistonPositions_.size()),
    xl_(readScalar(dynamicMeshCoeffs_.lookup("xLeft"))),
    xr_(readScalar(dynamicMeshCoeffs_.lookup("xRight"))),
    repetitions_(dynamicMeshCoeffs_.lookupOrDefault("repetitions",1.0)),
    timeInterpolation_(dynamicMeshCoeffs_.lookupOrDefault<word>("timeInterpolation","spline")),
    spaceInterpolation_(dynamicMeshCoeffs_.lookupOrDefault<word>("spaceInterpolation","linear")),
	yPistonCentres_(nPistons_),
	nMovingPoints_(sum(neg(points().component(vector::X)-xr_))),
	movingPoints_(nMovingPoints_),
	xOriginal_(nMovingPoints_),
	yOriginal_(nMovingPoints_),
	writePositionsToLogFile_(dynamicMeshCoeffs_.lookupOrDefault<bool>("positionsToLog",0)),
	amplificationFactor_(dynamicMeshCoeffs_.lookupOrDefault("amplification",1.0))
{
    Info<< "Creating dynamicWaveMakerFvMesh... " << endl;

	//Amplifying/attenuating signal
	pistonPositions_ = amplificationFactor_*pistonPositions_;
	
	//Finding moving mesh points
	const pointField& p = points();
	scalarField isMoving(neg(p.component(vector::X)-xr_));
	label iMov(-1);
	forAll(isMoving,ip)
	{
		if (isMoving[ip])
		{
			iMov++;
			movingPoints_[iMov] = ip;
			xOriginal_[iMov] = p[ip].component(vector::X);
			yOriginal_[iMov] = p[ip].component(vector::Y);
		}
	}
	
	//Reading piston y-positions
	if (dynamicMeshCoeffs_.found("yPiston")) 
	{
		yPistonCentres_ = dynamicMeshCoeffs_.lookup("yPiston");	
	} 
	else
	{
		Info << "yPiston not found. Distributing pistons uniformly between y = " 
			<< gMax(yOriginal_) << " and "
			<< gMin(yOriginal_) << endl;

		scalar pistonWidth = (gMax(yOriginal_)-gMin(yOriginal_))/nPistons_ ;
		forAll(yPistonCentres_,ip)
		{
			yPistonCentres_[ip] = gMin(yOriginal_) + (ip+0.5)*pistonWidth;
		}
	}
	if (writePositionsToLogFile_)
	{
		Info << "Piston y-positions: " << yPistonCentres_ << endl;
	}
	update();
	write();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicWaveMakerFvMesh::~dynamicWaveMakerFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicWaveMakerFvMesh::update()
{

	scalar timeBeforeMeshUpdate = time().elapsedCpuTime();
	
    //For times later than latest time in tabulated table the mesh is treated as static unless a number of repetitions is specified
	scalar t = time().value();
	scalar period = max(times_);
	if ( t <= repetitions_*period )
	{
			t = fmod(t,period);
	}

	//Copy of all mesh points to be modified and handed to fvMesh::movePoints below
	pointField p(points());
		
    if ( t <= max(times_) && t>=min(times_) ) 
    {
		//Finding piston positions xPist at current time by interpolation
		scalarField xPist(nPistons_);
		forAll(xPist,ip)
		{
			if ( timeInterpolation_ == "spline" )
			{
				xPist[ip] = Foam::interpolateSplineXY(t, times_, pistonPositions_[ip]);
			}
			else if (timeInterpolation_ == "linear" )
			{
				xPist[ip] = Foam::interpolateXY(t, times_, pistonPositions_[ip]);
			}
			else
			{
				Info << "Warning: unknown timeInterpolation method (options: linear and spline)" << endl;				
			}
		}
		
		if (writePositionsToLogFile_)
		{
			Info << "Piston x-positions: " << xPist << endl;
		}

		//Creating piston position x and y arrays for use in subsequent interpolation
		scalarField Xi(nPistons_+2), Yi(nPistons_+2);
		Xi[0] = xPist[0];
		Yi[0] = -1e10;
		forAll(xPist,ix) 
		{
			Xi[ix+1] = xPist[ix];
			Yi[ix+1] = yPistonCentres_[ix];
		}
		Xi[nPistons_+1] = xPist[nPistons_-1];
		Yi[nPistons_+1] = 1e10;

		scalarField XNEW(nMovingPoints_);
		if ( spaceInterpolation_ == "spline" )
		{
			XNEW = Foam::interpolateSplineXY(yOriginal_, Yi, Xi);
		}
		else if ( spaceInterpolation_ == "linear" )
		{
			XNEW = Foam::interpolateXY(yOriginal_, Yi, Xi);
		}
		else
		{
			Info << "Warning: unknown spaceInterpolation method (options: linear and spline)" << endl;
		}

		//Changing x-coordinate of all moving points
		scalar xmid(0.5*(xl_+xr_)), L(0.5*(xr_-xl_)), X(0.0);
		forAll(movingPoints_,ip)
		{	
			scalar x = xOriginal_[ip];
			if ( x < xl_ ) 
			{
				x = (x + XNEW[ip]);
			}
			else if ( x >= xl_ && x < xmid ) 
			{
				X = x - xl_;
				X = X*(1 - .5*(X/L)*(XNEW[ip]/L));
				x = X + xl_ + XNEW[ip];
			}
			else if ( x >= xmid && x < xr_ ) 
			{
				X = -(x - xr_);
				X = X*(1 - 0.5*(X/L)*(XNEW[ip]/L));
				x = xr_ - X;
			}
			p[movingPoints_[ip]].replace(vector::X, x);
		}
	}

	scalar newPositionCalculationTime = time().elapsedCpuTime() - timeBeforeMeshUpdate;

	fvMesh::movePoints(p);

	if (foundObject<volVectorField>("U"))
	{
		volVectorField& U =
			const_cast<volVectorField&>(lookupObject<volVectorField>("U"));
		U.correctBoundaryConditions();
	}
	if (foundObject<volVectorField>("alpha1"))
	{
		volScalarField& alpha1 =
			const_cast<volScalarField&>(lookupObject<volScalarField>("alpha1"));
		alpha1.correctBoundaryConditions();
	}

	Info<< "Time spent on mesh update: " <<  time().elapsedCpuTime() - timeBeforeMeshUpdate <<
			"s, hereof " << newPositionCalculationTime << "s on new point calculation." << endl;

    return true;
}


// ************************************************************************* //
