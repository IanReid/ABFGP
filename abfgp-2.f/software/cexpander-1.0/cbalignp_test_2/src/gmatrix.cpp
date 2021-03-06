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

#include "gmatrix.h"

void generalmatrix::clearmatrix()
{
	if ( Matrix.size() > 0 )
	{
		for ( int i=0; i != Matrix.size(); i++ )
			Matrix[i].clear();
		Matrix.clear();
	}
}

void generalmatrix::buildmatrix()
{
	iCols = Sequence1.iSequence.length() + 1;
	iRows = Sequence2.iSequence.length() + 1;
	
	Matrix.resize(iCols);
	for ( int i = 0; i != Matrix.size(); i++ )
		Matrix[i].resize( iRows ); //
}

int generalmatrix::getT( int i, int j )
{
	if ( i == 0 && j == 0 )
		return 0;
		
	if ( i == 0 )
	{
		if ( localOnly )
			return 0; // if added 22-10-2006
		else // else added 22-10-2006
			return 4; // -- original
	}
	if ( j == 0 )
	{
		if ( localOnly )
			return 0;
		else
			return 2; // same as for i 22-10-2006
	}
		
	int V = getV( i, j );
	int T = 0;
	
	if ( localOnly && V == 0 )
		return 0;
		
	if ( getM( i, j ) == V )
		T = T + 1; // diagonal
	if ( getIx( i, j ) == V )
		T = T + 2; // gap seq2
	if ( getIy( i, j ) == V )
		T = T + 4;  // gap Seq1
		
	return T;
}

void generalmatrix::initialfill()
{
	// nul point
	
	int w = -1000000;
	setM( 0, 0, 0 );
	setIx( 0, 0, 0 );
	setIy( 0, 0, 0 );
	setHgap( 0, 0, 0 );
	setVgap( 0, 0, 0 );
	
	// top  row
	
	for ( int i = 1; i != iCols; i++ )
	{
		if ( localOnly || !chargeOverHang )
		{
			setIx( i, 0, 0); //-(gapOpen + ( gapExtention * ( i - 1 ) ) ) );
			setM( i, 0, 0);//- (gapOpen + ( gapExtention * ( i - 1 ) ) ));
		}
		else
		{
			setM( i, 0, - (gapOpen + ( gapExtention * ( i - 1 ) ) ) );
			setIx( i, 0, - (gapOpen + ( gapExtention * ( i - 1 ) ) ) );
		}
		
		setIy( i, 0, w );
		setHgap(i,0, i ); // new
		setVgap(i, 0, 0 ); // new
		
	}
	
	for ( int j = 1; j != iRows; j++ )
	{
		if ( localOnly || !chargeOverHang )
		{
			setIy( 0, j, 0);//-(gapOpen + ( gapExtention * ( j - 1 ) ) ) );
			setM( 0, j, 0); //(gapOpen + ( gapExtention * ( j - 1 ) ) ) );
		}
		else
		{
			setM( 0, j, -(gapOpen + ( gapExtention * ( j - 1 ) ) ) );
			setIy( 0, j, -(gapOpen + ( gapExtention * ( j - 1 ) ) ) );
		}
		
		setIx( 0, j, w );
		setVgap( 0, j , j );
		setHgap( 0, j,  0 );
	}
}
		
void generalmatrix::matrixfill()
{
	int IxM,IxX;
	int IyM,IyY;
	
	int MM, MX, MY;
	
	int M, Ix, Iy;
	
	int besti = getV( iCols - 1, iRows - 1 );
	int bestj = besti;
	int pI = iCols - 1;
	int pJ = iRows - 1; 
	
	int localM = 0;
	int Q =0;
	
	for ( int j = 1; j != iRows; j++ )
		for ( int i = 1; i != iCols; i++ )
		{
			IxM = getM( i - 1, j ) - gapOpen;
			IxX = getIx( i - 1, j ) - gapExtention;
			Ix = max( IxM, IxX );
			 
			IyM = getM( i, j - 1 ) - gapOpen;
			IyY = getIy( i, j - 1 ) - gapExtention;
			Iy = max( IyM, IyY );
			
			MM = getM( i - 1, j - 1 );
			MX = getIx( i - 1, j - 1 );
			MY = getIy( i - 1, j - 1 );
			M = max( MM, max( MX, MY ) ) + score( i, j );
			
			if ( localOnly) // condition modified 22-10-2006
			{ 
				Q = max( M, (max (Iy, Ix )) );
				Q = max( Q, 0);
				if (Q == 0 )
				{
					M =0;
					Iy = 0;
					Ix = 0; 
				//M = max( M, 0 ); // original
				//Iy = max(Iy, 0); // tt
				//Ix = max(Ix, 0); //
						
				}
			} // end edit 22-10-2006 
			setIx( i, j, Ix );
			setIy( i, j, Iy );
			setM( i, j, M );
			
			if ( IxM >= IxX ) // >= instead of > 22-10-2006
				setHgap( i, j , 1 );
			else
				setHgap( i, j , getHgap( i - 1, j ) + 1 );
			
			if ( IyM >= IyY ) // >= instead of > 22-10-2006
				setVgap( i, j, 1 );
			else
				setVgap( i, j, getVgap( i, j - 1 ) + 1);
			
			// additional add 22-10-2006
			
			if ( localOnly && M == 0 && Iy == 0 && Ix ==0 )
			{
				setHgap( i, j, 0 );
				setVgap( i, j, 0 );
			}
			
			// end security addition 22-10-2006 
			
			if ( localOnly )
			{
				if ( getV( i, j ) > localM )
				{
					localM = getV( i, j );
					pI = i;
					pJ = j;
				}
			}
			else	
			if ( !chargeOverHang )
			{
				if ( i == iCols - 1 ) 
				{
					if ( getV( i, j ) > bestj )
					{
						bestj = getV( i, j );
						pJ = j;
					}
				}
				
				if ( j == iRows - 1 )
				{
					if ( getV( i, j ) > besti )
					{
						besti = getV( i, j );
						pI = i;
					}
				}
			}	
		
			if ( localOnly )
			{
				if ( getV( i,j ) > localM )
				{
					localM = getV( i, j );
					pI= i;
					pJ= j;
				}
			}
		}
		
		if ( localOnly )
		{
			local_i = pI;
			local_j = pJ;
		}
		else
		if ( !chargeOverHang )
		{
			if ( besti > bestj )
			{
				overhang_i = pI;
				overhang_j = iRows - 1;
			}
			else
			{
				overhang_i = iCols - 1;
				overhang_j = pJ;
			}
		}
}	

void generalmatrix::traceback()
{
	alignment.a1 = "";
	alignment.a2 = "";
	alignment.mL = "";
	
	int i, j;
	int T;
	int temp1, temp2;
	
	if ( localOnly )
	{
		alignment.score = getV( local_i, local_j );
		i = local_i;
		j = local_j;
		
		alignment.a1_end = i;
		alignment.a2_end = j;
	}
	else	
	{
		if ( showSkip )
		{
			i = iCols - 1;
			j = iRows - 1; 
		
			alignment.a1_end = i;
			alignment.a2_end = j;
		
			if ( !chargeOverHang )
			{
				if ( overhang_i < iCols - 1 )
				{
			
					while ( i != overhang_i )
					{
						alignment.a1+= Sequence1.iSequence[i - 1];
						alignment.a2+= '-';
						alignment.mL+= ' ';
						i--;
					}
				}
				else
				if ( overhang_j < iRows - 1 )
				{
							
				
					while ( j != overhang_j )
					{
						alignment.a1+='-';
						alignment.a2+= Sequence2.iSequence[ j - 1 ];
						alignment.mL+= ' ';
						j--;
					}
				}
			}
		}
		else
		{
			i = overhang_i;
			j = overhang_j;
			
			alignment.a1_end = i;
			alignment.a2_end = j;
		}
	}	
	T = getT( i, j );
	alignment.score = getV( i, j );
	while (T != 0 )
	{
			if ( ( T & 1 ) == 1 )
		{
			alignment.a1+= Sequence1.iSequence[i - 1];
			alignment.a2+= Sequence2.iSequence[j - 1];
			alignment.mL+= mType( i, j ); 
			i--;
			j--;	
		}
		else
		
		if ( ( T & 2 ) == 2 )
		{
			temp2 = 0;
			temp1 = getHgap( i, j );
			while ( temp2 != temp1 )
			{
				alignment.a1+= Sequence1.iSequence[ i - 1];
				alignment.a2+= '-';
				alignment.mL+= ' ';
				i--;
				temp2++;
			}
		}
		else
		if ( ( T & 4 ) == 4 )
		{
			temp2 = 0;
			temp1 = getVgap( i, j );
			while ( temp1 != temp2 )
			{
				alignment.a1+='-';
				alignment.a2+= Sequence2.iSequence[ j - 1 ];
				alignment.mL+=' ';
				j--;
				temp2++;
			}
		} 
		
		if ( !showSkip )
		{
			if ( i == 0 || j == 0 )
				T = 0;
			else T = getT( i, j );
		}
		else	T = getT( i, j );
	}

	alignment.a1 = revString( alignment.a1 );
	alignment.a2 = revString( alignment.a2 );
	alignment.mL = revString( alignment.mL );
	
	alignment.a1_start = 1 + alignment.a1_end - cCount_NNot( alignment.a1, '-' );
	alignment.a2_start = 1 + alignment.a2_end - cCount_NNot( alignment.a2, '-' ); 
	
	alignment.a1_alignedP = cCount_NNot( alignment.a1, '-' );
	alignment.a2_alignedP = cCount_NNot( alignment.a2, '-' );
}

void generalmatrix::setSequences( generalSequence S1, generalSequence S2 )
{
	//if ( S1.iSequence.length() > S2.iSequence.length() )
	//{
	//	Sequence1 = S1;
	//	Sequence2 = S2;
	//} 
	//else
	//{
	//	Sequence1 = S2;
	//	Sequence2 = S1;
	//}
	Sequence1 = S1; // new 6-10-2007 Just keep it the same
	Sequence2 = S2; 
}

void generalmatrix::doalignment()
{
	clearmatrix();
	buildmatrix();
	initialfill();
	matrixfill();
	traceback();
	calcScores(); // calculate identity positives gaps 
	calcselfsims();
}

void generalmatrix::setmodus( bool lonly, bool coverhang, bool showskip )
{
	localOnly = lonly;
	chargeOverHang = coverhang;
	showSkip = showskip;
}

void generalmatrix::dumpAlignment( int colW )
{
	cout << "Title:" << myTitle << endl;
	cout << endl;

	cout << "The Sequences:\n\n";
	cout << "Sequence1: ";
	cout << Sequence1.iHeader << endl;
	cout << "length1: ( " << Sequence1.iSequence.length() << " )" << endl;
	cout << endl;
	
	cout << "Sequence2: ";
	cout << Sequence2.iHeader << endl;
	cout << "length2: ( " << Sequence2.iSequence.length() << " )" << endl;
	cout << endl;
	
	cout << "The alignment:\n\n";
	
	int nc = alignment.a1.length() / colW;
	if ( alignment.a1.length() % colW != 0 )
		nc++; 
	
	int o = Idigs( max( alignment.a1_end, alignment.a2_end ) );
	
	string a1;
	string a2;
	
	int s1 = alignment.a1_start - 1;
	int s2 = alignment.a2_start - 1; 
	
	int l, p;
	int k = 0;
	
	while ( k < nc )
	{
		p = k * colW; 
		
		if ( p + colW <= alignment.a1.length() )
			l = colW;
		else
			l = alignment.a1.length() - p; 
		
		a1 = alignment.a1.substr( p, l );
		a2 = alignment.a2.substr( p, l );
		
		cout << "S1:" << iSpace(5);
		cout << iSpace( o - Idigs( s1 + 1 ) ) << s1 + 1 << iSpace(5);
		cout << a1 << iSpace(5) << s1 + cCount_NNot( a1, '-' ) << endl;
		
		cout << "Ml:" << iSpace(5); 
		cout << iSpace( o ) << iSpace(5);
		cout << alignment.mL.substr( p, l ) << endl;
		
		cout << "S2:" << iSpace(5);
		cout << iSpace( o - Idigs( s2 + 1 ) ) << s2 + 1 << iSpace(5);
		cout << a2 << iSpace(5) << s2 + cCount_NNot( a2, '-' ) << endl;
		cout << endl;
		
		s1+=cCount_NNot( a1, '-' );
		s2+=cCount_NNot( a2, '-' );
		k++; 
	}
	
	float Aid = 100 *  float(alignment.identities) / float(alignment.a1.length() ) ;
	float AiP = 100 *  float(alignment.positives) / float( alignment.a1.length() );
	float AiG = 100 *  float(alignment.gaps) / float( alignment.a1.length() ) ;
	cout << "Score: " << alignment.score << endl;
	cout << "Alignment length: " << alignment.a1.length() << endl;
	
	cout << "identities: " << alignment.identities << "/" << alignment.a1.length();
	cout << " (" << Aid << " ) " << endl;

	cout << "Positives: " << alignment.positives << "/" << alignment.a1.length();
	cout << " (" << AiP << " ) " << endl;

	cout << "Gaps: " << alignment.gaps << "/" << alignment.a1.length();
	cout << " (" << AiG << " ) " << endl;
	
	cout << endl;
	cout << "alignment type: ";
	
	if ( localOnly )
		cout << "local" << endl;
	else
	{
		cout << "global" << endl;
		cout << "overhang charge: ";
		if ( chargeOverHang )
			cout << "yes" << endl;
		else
			cout << "no" << endl;
	}
	
	cout << "Gap open: " << gapOpen << endl;
	cout << "Gap extention: " << gapExtention << endl;
	cout << endl;
	cout << "self_sim_1:" << selfSc_1 << endl;
	cout << "sRatio_1: " << sRatio_1 << endl;
	cout << endl;
	cout << "self_sim_2:" << selfSc_2 << endl;
	cout << "sRatio_2: " << sRatio_2 << endl;
	cout << "meanR: " << ( sRatio_1 + sRatio_2 ) / 2 << endl; 
	cout << endl;
	display_specific();

}

void generalmatrix::dumpAlignment_COMPUTERS()
{
	cout << "Title:" << myTitle << endl;
	cout << endl;
	cout << "The Sequences:\n\n";
	cout << "Sequence1: ";
	cout << Sequence1.iHeader << endl;
	cout << "length1:" << Sequence1.iSequence.length()  << endl;
	cout << endl;
	
	cout << "Sequence2: ";
	cout << Sequence2.iHeader << endl;
	cout << "length2:" << Sequence2.iSequence.length() << endl;
	cout << endl;
	cout << "alignment:" << endl << endl;
	cout << "a1:" << alignment.a1 << endl;
	cout << "mL:" << alignment.mL << endl;
	cout << "a2:" << alignment.a2 << endl;
	cout << endl;
	cout << "Specs:" << endl << endl;
	cout << "a1_start:" << alignment.a1_start << endl;
	cout << "a2_start:" << alignment.a2_start << endl;
	cout << "a1_end  :" << alignment.a1_end << endl;
	cout << "a2_end  :" << alignment.a2_end << endl;
	
	float wP1 = float(alignment.a1_alignedP) / float(Sequence1.iSequence.length() );
	float wP2 = float(alignment.a2_alignedP) / float(Sequence2.iSequence.length() );
	
	cout << "a1_P: " << alignment.a1_alignedP << "," << wP1 << endl;
	cout << "a2_P: " << alignment.a2_alignedP << "," << wP2 << endl; 
	
	cout << endl;
	
	
	float Aid = 100 *  float(alignment.identities) / float(alignment.a1.length() ) ;
	float AiP = 100 *  float(alignment.positives) / float( alignment.a1.length() );
	float AiG = 100 *  float(alignment.gaps) / float( alignment.a1.length() ) ;
	
	
	cout << "Scores: " << endl << endl;
	
	cout << "Score: " << alignment.score << endl;
	cout << "Alignment length: " << alignment.a1.length() << endl;
	
	cout << "identities: " << alignment.identities << "/" << alignment.a1.length();
	cout << "," << Aid << endl;

	cout << "Positives: " << alignment.positives << "/" << alignment.a1.length();
	cout << "," << AiP << endl;

	cout << "Gaps: " << alignment.gaps << "/" << alignment.a1.length();
	cout << "," << AiG << endl;
	cout << endl;
	cout << "self_sim_1:" << selfSc_1 << endl;
	cout << "sRatio_1: " << sRatio_1 << endl;
	cout << endl;
	cout << "self_sim_2:" << selfSc_2 << endl;
	cout << "sRatio_2: " << sRatio_2 << endl;
	cout << "meanR: " << ( sRatio_1 + sRatio_2 ) / 2 << endl; 
	cout << endl;
	cout << "alignment type: ";
	
	if ( localOnly )
		cout << "local" << endl;
	else
	{
		cout << "global" << endl;
		cout << "overhang charge: ";
		if ( chargeOverHang )
			cout << "yes" << endl;
		else
			cout << "no" << endl;
	}
	
	cout << endl;
	display_specific();

	
	
}

void generalmatrix::calcselfsims()
{
	// some adhoc measures are taken here 2- 10 - 2006
	
	selfSc_1 = 0;
	generalSequence tmp1 = Sequence1;
	generalSequence tmp2 = Sequence2;
	string tmp = Sequence2.iSequence;
	Sequence2.iSequence = tmp1.iSequence;
	for ( int i = alignment.a1_start ; i != alignment.a1_end + 1; i++ )
		selfSc_1+=score( i , i );
	selfSc_2 = 0;
	Sequence2.iSequence = tmp;
	tmp = Sequence1.iSequence;
	Sequence1.iSequence = tmp2.iSequence;
	Sequence2.iSequence = tmp2.iSequence;// Sequence2.iSequence;
	for ( int i = alignment.a2_start; i != alignment.a2_end + 1; i++ )
		selfSc_2+=score( i, i );
			
	sRatio_1 = alignment.score / selfSc_1;
	sRatio_2 = alignment.score / selfSc_2; 
	Sequence2.iSequence = tmp;
	Sequence1 = tmp1;
	Sequence2 = tmp2;
}
