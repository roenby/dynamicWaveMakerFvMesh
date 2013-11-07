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

Class
    Foam::dynamicWaveMakerFvMesh

Description
    Mesh motion for wave maker with prescribed positions of the left 
    wall (in the X-direction) to make the wall mimic 3D wave maker piston 
    motion.
     
    Based on dynamicInkJetFvMesh.

SourceFiles
    dynamicWaveMakerFvMesh.C

Author
    Johan Roenby, DHI.  All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef dynamicWaveMakerFvMesh_H
#define dynamicWaveMakerFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicWaveMakerFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicWaveMakerFvMesh
:
    public dynamicFvMesh
{
    // Private data
		
		//- Dictionary to hold contents from <case>/constant/dynamicMeshDict
        dictionary dynamicMeshCoeffs_;

		//- Array of times
		scalarField times_;

        //- Array of piston positions corresponding to times_
		List<scalarField> pistonPositions_;
		
		//- Number of pistons
		int nPistons_;

		//- Left bound of region with extracting/contracting cells
		scalar xl_;
		
		//- Right bound of region with extracting/contracting cells
		scalar xr_;

        //- Number of times to repeat tabulated piston motion (default = 1) 
		scalar repetitions_;

        //- Interpolation method for piston positions between listed times (default = spline) 
		word timeInterpolation_;

        //- Interpolation method for positions between listed piston positions (default = linear) 
		word spaceInterpolation_;

		//- Y coordinate of piston centres 
		scalarField yPistonCentres_;

		//- Number of moving mesh points 
		int nMovingPoints_;	
		
		//- Indices of moving points in points file
		labelList movingPoints_;
		
		//- Original x-coordinates of moving points
		scalarField xOriginal_;
		
		//- y-coordinate of moving points
		scalarField yOriginal_;

		//- Should piston positions be written to log (default: no) 
		bool writePositionsToLogFile_;

		//- Optional factor multiplied on piston positions to amplify aor attenuate them
		scalar amplificationFactor_;
		
		// Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicWaveMakerFvMesh(const dynamicWaveMakerFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicWaveMakerFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicWaveMakerFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicWaveMakerFvMesh(const IOobject& io);


    //- Destructor
    ~dynamicWaveMakerFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
