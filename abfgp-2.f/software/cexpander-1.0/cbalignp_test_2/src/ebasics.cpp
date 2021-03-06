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

#include "ebasics.h"

vector< string > split( string& S, char c, bool keepEmpty )
{
	vector < string > theSubs;
	signed int lastPos = -1;
	signed int curPos = 0;
	signed int dst;
	signed int theLength = S.length();
	
	while ( curPos <= theLength )
	{
		if ( S[curPos] == c || curPos == theLength )
		{
			dst = ( curPos - lastPos ) - 1; 
			if ( dst < 1 )
			{
				if ( keepEmpty )
					theSubs.push_back("");
			}
			else
			{
				theSubs.push_back( S.substr( lastPos + 1, dst ) );
			}
			lastPos = curPos; 
		}
		curPos++;
	}
	return theSubs;
}


vector< string > split_low( string& S, bool keepEmpty )
{
	vector< string> theSubs;
	signed int lastPos = -1;
	signed int curPos = 0;
	signed int dst;
	signed int theLength = S.length();
	
	while ( curPos <= theLength )
	{
		if ( S[curPos] < 33 || curPos == theLength )
		{
			dst = ( curPos - lastPos ) - 1;
			if ( dst < 1 )
			{
				if ( keepEmpty )
					theSubs.push_back("");
			}
			else
			{
				theSubs.push_back( S.substr( lastPos + 1, dst ) );
			}
			lastPos = curPos; 
		}
		curPos++;
	}
	
	return theSubs;
} 

int scanLow( string& S, int pos, int step )
{
	while ( pos > -1 && pos < S.length() )
	{
		if ( S[pos] < 33 )
			pos = pos + step;
		else
			break;
	}
	return pos;
}

int scanFor( string& S, int pos, char c,int step )
{
	while ( pos > -1 && pos < S.length() )
	{
		if ( S[pos] != c ) 
			pos = pos + step;
		else
			break;
	} 
	return pos;
}

int scanNot( string& S, int pos, char c, int step ) 
{
	while ( pos > -1 && pos < S.length() )
	{
		if ( S[pos] == c )
			pos = pos + step;
		else
			break;
	}
	return pos;
}

string stripL_low( string S ) 
{
	int L = scanLow( S, 0, 1); // positive direction 
	
	if ( L != S.length() )
		S = S.substr(L, S.length() - L );
	return S;
}

string stripR_low( string S )
{
	int R = scanLow( S, S.length() - 1, - 1);
	
	if ( R != -1 )
		S = S.substr( 0, R + 1);
	return S; 
} 

string strip_low( string S)
{
	S = stripL_low( S );
	S = stripR_low( S );
	return S;
}

string stripL( string S, char c )
{
	int L = scanNot( S, 0, c, 1);
	
	if ( L != S.length() )
		S = S.substr( L, S.length() - L );
	return S;
}

string stripR( string S, char c ) 
{
	int R = scanNot(S, S.length() - 1, c, -1);
	
	if ( R != -1 )
		S = S.substr(0, R + 1);
	return S;
}

string strip( string S, char c )
{
	S = stripL( S, c );
	S = stripR( S, c );
	return S;
}

string revString( string S )
{
	int M = S.length() / 2;
	char tmp;
	int pos = 0;
	int rpos = S.length() - 1; 
	
	while ( pos != M )
	{
		tmp = S[pos];
		S[pos] = S[rpos];
		S[rpos] = tmp; // swapping 
		pos++;
		rpos--;
	}
	return S;
} 
	
int cCount(string& S, char ch)
{
	int c = 0;
	for ( int i = 0; i != S.length(); i++)
		if ( S[i] == ch )
			c++;
	return c;
}

int cCount_NNot(string& S, char ch)
{
	return S.length() - cCount(S,ch );
}

string stringOf(int N, char ch)
{
	string rt = "";
	
	if ( N > 0 )
	{
		char* u = new char[N+1];
		u[N] = 0;
		for (int k = 0; k != N;k++)
			u[k] = ch;
		rt=u;
		delete[] u;
	}
	return rt;
}
	
string iSpace( int N )
{
	return stringOf( N, ' ' );
}

int Idigs( int N )
{
	int w = 0;
	while ( N / 10 > 0 )
	{
		w++;
		N = N / 10;
	}
	w++;
	return w;
}

vector<generalSequence> readFastaFile(char* Fname)
{
	vector<generalSequence> Seqs;
	ifstream myFile; 
	
	myFile.open( Fname );
	
	if ( myFile.is_open() )
	{
		string Istring;
		generalSequence sequence;
		
		sequence.iHeader = "";
		sequence.iSequence = "";
		
		while ( !myFile.eof() )
		{
			getline( myFile, Istring );
			// check if it's a sequence-header 
			if ( Istring.size() > 0 )
			{ 
				if ( Istring[0] == '>' )
				{
					// do we need to store 
					if ( sequence.iHeader != "" && sequence.iSequence != "" )
						Seqs.push_back( sequence );
					
					sequence.iHeader = Istring; // reset ( new sequence ) 
					sequence.iSequence = ""; // reset 
				}
				else
					sequence.iSequence = sequence.iSequence + Istring; // elongate 
			}
		}
		
		// check ones more if we need to store 
		if ( sequence.iHeader != "" && sequence.iSequence != "" )
			Seqs.push_back( sequence ); 
		
		myFile.close();
	}
			
	return Seqs;
}

generalSequence nextFasta( ifstream& leStream )
{
  generalSequence current;
  
  current.iHeader = "";
  current.iSequence = "";
  
  while ( !leStream.eof() )
  {
      char ch = leStream.peek(); //test
      
      if ( ch == '>' )
      {
	if ( current.iHeader != "" )
	  break; // done
	else
	  getline( leStream, current.iHeader );
      }
      else
      {
	string Istring;
	
	getline( leStream, Istring );
	
	
	if ( Istring.size() > 0 && current.iHeader != "" )
	{
	  current.iSequence+= Istring;
	}
      }
  }
  
  return current;
}

	
clineItem clI( string com, string Idefault, bool switched ,int  IvDefault,short T )
{
	clineItem cc;
	cc.com = com;
	cc.Ivalue = IvDefault;
	cc.switched = switched;
	cc.Istring = Idefault;
	cc.T = T;
	return cc; 
} 

int gItemNr( clineItemList& L, string com )
{
	int fnd = -1;
	int p = 0;
	while ( p < L.size() && fnd == -1 )
	{
		if ( L[p].com == com )
			fnd = p;
		else
			p++;
	}
	return fnd; 
} 


string handleCommandLine( clineItemList& Ilist, int argc, char* argv[]  )
{
	cout << Ilist.size() << endl;
	int expect = 0;
	int p = 1;
	int F = 0;
	int T = 0;
	
	while ( p < argc )
	{
		if ( expect == 0 )
		{
			F = gItemNr( Ilist, argv[p] );
			if ( F < 0 )
			{
				string rt = argv[p];
				return "Unknown command: " + rt;
			}		 
			
			T = Ilist[F].T;
			Ilist[F].switched = true;
			
			if ( T > 0 )
				expect = 1;
		}
		else
		{
			 
			switch(T)
			{
				case 1:
					Ilist[F].Istring = argv[p]; 
					break;
				case 2:
					Ilist[F].Ivalue = atoi(argv[p] );
					break;
			} 
			expect = 0;
		}
		p++;
	}
	return "";
}
	

TPos cPos( int x, int y )
{
	TPos P;
	P.x = x;
	P.y = y;
	return P;
}
