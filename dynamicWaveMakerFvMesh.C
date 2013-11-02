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
    yMin_(readScalar(dynamicMeshCoeffs_.lookup("yMin"))),
    yMax_(readScalar(dynamicMeshCoeffs_.lookup("yMax"))),
    repetitions_(dynamicMeshCoeffs_.lookupOrDefault("repetitions",1.0)),
    timeInterpolation_(dynamicMeshCoeffs_.lookupOrDefault<word>("timeInterpolation","spline")),
    spaceInterpolation_(dynamicMeshCoeffs_.lookupOrDefault<word>("spaceInterpolation","linear")),
	yPistonCentres_(nPistons_),
    pointsAssociatedWithPistons_(nPistons_),
    xOrg_(nPistons_)
{
    Info<< "Creating dynamicWaveMakerFvMesh... " << endl;
	
	//Calculating y coordinates of piston centre positions
	scalar pistonWidth = (yMax_-yMin_)/nPistons_;
	forAll(yPistonCentres_,iPist)
	{
		yPistonCentres_[iPist] = yMin_ + (.5+iPist)*pistonWidth;
	}

	//Identifying which of the moving points that belong to which piston
	const pointField& p = points();
	forAll(p,ip)
	{
		scalar xp = p[ip].component(vector::X);
		scalar yp = p[ip].component(vector::Y);
		if ( xp <= xr_ && (yp >= yMin_ && yp <= yMax_) ) //Points not satisfying this condition (typically most points) will be stationary.
		{
			scalar yp = p[ip].component(vector::Y) - 1e-6*(yMax_-yMin_);
			label pistonInd = max(0,floor((yp-yMin_)/pistonWidth)); //Together the small value subtracted from yp and the max(0,...) prevents pistonInd>nPistons_-1 when yp=yMax_. 
			pointsAssociatedWithPistons_[pistonInd].append(ip);
			xOrg_[pistonInd].append(xp);
		}
	}
	forAll(xOrg_,iPist) 
	{
		pointsAssociatedWithPistons_[iPist].shrink();
		xOrg_[iPist].shrink();
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
		
		//Creating piston position x and y arrays for use in subsequent interpolation
		scalarField Xi(nPistons_+2), Yi(nPistons_+2);
		Xi[0] = xPist[0];
		Yi[0] = yMin_;
		forAll(xPist,ix) 
		{
			Xi[ix+1] = xPist[ix];
			Yi[ix+1] = yPistonCentres_[ix];
		}
		Xi[nPistons_+1] = xPist[nPistons_-1];
		Yi[nPistons_+1] = yMax_;
		
		forAll(yPistonCentres_,iPist)
		{
			labelList iPoint = pointsAssociatedWithPistons_[iPist];
			scalarField Yp(iPoint.size());
			forAll(Yp,iy)
			{
				Yp[iy] = p[iPoint[iy]].component(vector::Y);
			}
			//The n'th element in Xp below will be the piston position in front of the moving mesh point with index iPoint[n].
			scalarField Xp(Yp.size());
			if ( spaceInterpolation_ == "spline" )
			{
				Xp = Foam::interpolateSplineXY(Yp, Yi, Xi);
			}
			else if ( spaceInterpolation_ == "linear" )
			{
				Xp = Foam::interpolateXY(Yp, Yi, Xi);
			}
			else
			{
				Info << "Warning: unknown spaceInterpolation method (options: linear and spline)" << endl;
			}
			
			//Changing x-coordinate of all moving points
			scalar xmid(0.5*(xl_+xr_)), L(0.5*(xr_-xl_)), X(0.0);
			forAll(iPoint,ip)
			{	
				scalar x = xOrg_[iPist][ip];
				if ( x < xl_ ) 
				{
					x = (x + Xp[ip]);
				}
				else if ( x >= xl_ && x < xmid ) 
				{
					X = x - xl_;
					X = X*(1 - .5*(X/L)*(Xp[ip]/L));
					x = X + xl_ + Xp[ip];
				}
				else if ( x >= xmid && x < xr_ ) 
				{
					X = -(x - xr_);
					X = X*(1 - 0.5*(X/L)*(Xp[ip]/L));
					x = xr_ - X;
				}
				p[iPoint[ip]].replace(vector::X, x);
			}
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
