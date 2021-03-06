/***************************************************************************
 *   Copyright (C) 2006 by edouard severing   *
 *   edouard@localhost   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "defragM.h"


void traceM_clear( TtraceMatrix &M )
{
	if ( M.size() > 0 )
	{
		for ( int j = 0; j != M.size(); j++ )
			M[j].clear();
		M.clear();
	}
}
	
void traceM_build( TtraceMatrix &M, int Cols, int Rows )
{
	M.resize( Cols );
	for ( int i = 0; i != Cols; i++ )
		M[i].resize( Rows );
}

void traceM_init( TtraceMatrix &M)
{
	int xM = M.size();
	int yM = M[0].size();
	
	for ( int x = 0; x != xM; x++ )
	{
		for ( int y = 0; y != yM; y++ )
			M[x][y] = 0;
	}		
}

bool needsExtention( std::string &S )
{
	if ( S.length() > 0 )
		if ( S[S.length() - 1] == '0' )
			return false;
	return true;
}

TPos calcEndPos( Tpath& P )
{
	int i = P.Starti;
	int j = P.Startj;
	char w;
	
	for ( int k = 0; k != P.path.length(); k++ )
	{
		w = P.path[k] - 48;
			
		if ( ( w & 1 ) == 1 || ( w & 2 ) == 2 )
			i--;
		if ( ( w & 1 ) == 1 || ( w & 4 ) == 4 )
			j--;
	}
	
	return cPos( i, j );
}

Tpath traceM_buildSeed( generalmatrix &Mat )
{
	Tpath seed;
	
	if ( Mat.mode_localOnly() )
	{
		// were local hence seed = ""
		seed.Starti = Mat.locali();
		seed.Startj = Mat.localj(); 
		seed.path = "";
	}
	else
	{
		if ( Mat.mode_chargeOverHang() )
		{
			// easy start at the end of the matrix 
			seed.Starti = Mat.nCols() - 1;
			seed.Startj = Mat.nRows() - 1;
			seed.path = ""; // no start path 
		}
		else
		{
			if ( Mat.overhangi() < Mat.nCols() - 1 )
			{
				seed.Startj = Mat.nRows() - 1;
				
				if ( Mat.mode_showSkip() )
				{
					seed.path = stringOf( Mat.nCols() - Mat.overhangi() - 1, '2' );
					seed.Starti = Mat.nCols() - 1;
				}
				else
				{
					seed.Starti = Mat.overhangi();
					seed.path = "";
				}
			}
			else
				if ( Mat.overhangj() < Mat.nRows() - 1 )
				{
					seed.Starti = Mat.nCols() - 1;
						
					if ( Mat.mode_showSkip() )
					{
						seed.path = stringOf( Mat.nRows() - Mat.overhangj() - 1, '4' );
						seed.Startj = Mat.nRows() - 1; 
					}
					else
					{
						seed.Startj = Mat.overhangj();
						seed.path = "";
					}
				}
				else
				{
					seed.Starti = Mat.nCols() - 1;
					seed.Startj = Mat.nRows() - 1;
					seed.path = "";
				}
		}
	}
	
	return seed;
}
	
void traceM_construct( TmatrixPaths &paths, generalmatrix &Mat, int mx )
{
	paths.clear(); // make sure its empty 
	
	int i, j;
	std::string seed = "";
	
	if ( Mat.mode_localOnly() )
	{
		i = Mat.locali();
		j = Mat.localj();
	}
	else
		if ( !Mat.mode_chargeOverHang() )
		{
			
			if ( Mat.mode_showSkip() )
			{
				if ( Mat.overhangi() < Mat.nCols() - 1 )
				{
					seed = stringOf( Mat.nCols() - Mat.overhangi(), '2' );
					i = Mat.overhangi();
					j = Mat.nRows() - 1;
				}
				else
					if ( Mat.overhangj() < Mat.nRows() - 1 )
					{
						seed = stringOf( Mat.nRows() - Mat.overhangj(), '4' );
						i = Mat.nCols() - 1;
						j = Mat.overhangj(); 
					}
					else
					{
						i = Mat.nCols() - 1;
						j = Mat.nRows() - 1;
					}
			}
			else
			{
				i = Mat.nCols() - 1;
				j = Mat.nRows() - 1;
			} 				
		}
		else
		{
			i = Mat.nCols() - 1;
			j = Mat.nRows() - 1;
		}
		
	//seed+= 48 + Mat.getT( i , j ); 
	Tpath cPath;
	Tpath tmp;
	cPath.Starti = i;
	cPath.Startj = j;
	cPath.path = seed;
	
	
	// test
	cPath = traceM_buildSeed( Mat );
	paths.push_back( cPath );
	int nr = 0;
	int exs = 1;
	char T;
	TPos P;
	int mg;
	
	while ( exs != 0 )
	{
		exs = 0;
		nr = 0;
		while (nr != paths.size() )
		{
			cPath = paths[nr];
			if ( 	needsExtention( cPath.path ) )
			{
				exs++;
				P = calcEndPos( cPath );
				T = Mat.getT( P.x, P.y );
						
				mg = 0;
				
				if ( T == 0 )
					
				{
					paths[nr].path+= '0';
				}
				
				if ( ( T & 1 ) == 1 )
				{
					mg++;
					paths[nr].path+= '1';
				}
				if ( ( T & 2 ) == 2 )
				{
					if ( mg > 0 )
					{
						//cout << "doub" << endl;
						if ( paths.size() < mx )
						{
							tmp = cPath;
							tmp.path+=stringOf( Mat.getHgap( P.x, P.y ), '2' );
							paths.push_back( tmp );
							//cout << paths.size() << endl;
						}
					}
					else
						paths[nr].path+=stringOf( Mat.getHgap( P.x, P.y ), '2' );
					mg++;
				}
				if ( ( T & 4 ) == 4 )
				{
					if ( mg > 0 )
					{
						//cout << "doub" << endl;
						if ( paths.size() < mx )
						{
							
							tmp = cPath;
							tmp.path += stringOf( Mat.getVgap( P.x, P.y ), '4' );
							paths.push_back( tmp );
							//cout << paths.size() << endl;
						}
					}
					else
						paths[nr].path+= stringOf( Mat.getVgap( P.x, P.y), '4' );
						
				}
			}
			nr++;
		
		}
	}

}

int pointExtention( TtraceMatrix &M, int i, int j )
{
	// left extention
	int uC = 0;
	int lC = 0;
	int rC = 0;
	int dC = 0;
	int lU = 0;
	int rD = 0;
	
	
	int sI = i;
	int sJ = j;
	
	int mI = M[0].size();
	int mJ = M.size();
	
	if ( M[sI - 1][sJ - 1] > 0 )
		lU = 1;
	if ( M[sI + 1][sJ + 1] > 0 )
		rD = 1;
		 
	for ( int I = sI - 1; I != 0; I-- )
	{
		//if ( M[I][sJ] > 0 )
		//	lC = 1;
		//cout << "Z";
	}
	cout << "A";
	
	for ( int I = sI + 1; I != mI; I++ )
	{
		//if ( M[I][sJ] > 0 )
		//	rC = 1;
	}
	
	cout << "B";
	for ( int J = sJ - 1; J != 0; J-- )
	{ 
		if ( M[sI][J] > 0 )
			uC = 1;
	}
	cout << "C";
	for ( int J = sJ + 1; J != mJ; J++ )
	{
		if ( M[sI][J] > 0 )
			dC = 1;
	}
	cout << "D";
	
	int T = lC + rC + uC + dC;
	
	if ( T > 0 )
		T = 5;
	else
		T = 0;
		
	return T;

}
void traceM_plot( TtraceMatrix &M, Tpath &p )
{
	int i = p.Starti;
	int j = p.Startj;
	
	int z;
	
	for ( int k = 0; k != p.path.length(); k++ )
	{
		z = p.path[k] - 48;
		
		M[i][j] = M[i][j] | z; // exit code
		
		if ( z == 1	|| z == 2 )
			i--;
		if ( z == 1 || z == 4 )
			j--;
		
		M[i][j] = M[i][j] | ( z * 8 ); // incoming code
	}
	
}	
void traceM_fill( TtraceMatrix &M, TmatrixPaths &paths )
{
	for ( int k = 0; k != paths.size(); k++ )
	{
		traceM_plot( M, paths[k] );
	}
} 
	// from here strat single of trace
	
void traceM_Draw( TtraceMatrix &M, generalmatrix &Mat )
{
 //
 int AA = 0;
 
 vector< int > Iline;
 vector< int > Jline;
 
 Iline.resize( Mat.nCols() );
 Jline.resize( Mat.nRows() );
 
 for ( int i = 0; i != Mat.nCols(); i++ )
 {
 	Iline[i] = 0;
 }
  
 for ( int j = 0; j != Mat.nRows(); j++ )
 {
 	Jline[j] = 0;
 }
 
 for ( int i = 0; i != Mat.nCols(); i++ )
 {
 	for ( int j = 0; j != Mat.nRows(); j++ )
 	{
 		if ( M[i][j] > 0 )
 		{
 			Jline[j] = Jline[j] + 1;
 			Iline[i] = Iline[i] + 1;
 		}
 	}
 }
 
 int D = 0;
 
 for ( int j = 0; j != Mat.nRows(); j++)
 {
 	//cout << j << endl;
 	for ( int i = 0; i != Mat.nCols(); i++ )
 	{
 		if ( i > 1 && i < Mat.nCols()  -2 && j > 1 && j < Mat.nRows() -2)
 		{
 			//AA = pointExtention( M, i, j );
 			AA = 0;
 			if ( M[i-1][j-1] > 0 )
 				D = 1;
 			else
 				D = 0;
 				
 			if ( M[i+1][j+1] > 0 )
 				D = D + 1;
 				 
 			if ( Jline[j]  > 1 && Iline[i] > 1 )
 				AA = 5;
 			if ( D > 1 )
 			{
 				if ( Jline[j] > 1 || Iline[i] > 1 )
 					AA = 5;
 			}
 		}
 		else
 		{
 			AA = 0; 
 		}	
 		if ( M[i][j] == 0 )
 			cout << ' ';
 		else
 		{
 			if ( AA == 0 )
 			{
 				cout << Mat.sSequence1().iSequence[i - 1];
 			}
 			else
 				cout << "+";
 		}
 	}
 	cout << endl;
 }
}
