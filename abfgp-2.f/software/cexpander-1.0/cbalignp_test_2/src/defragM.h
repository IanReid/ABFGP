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

#ifndef __defragMH
#define __defragMH

#include "gmatrix.h"
#include "ebasics.h"

typedef vector< vector < char > > TtraceMatrix;

struct Tpath{
	int Starti;
	int Startj;
	std::string path;
};

typedef vector < Tpath > TmatrixPaths;

void traceM_clear( TtraceMatrix &M );

void traceM_build( TtraceMatrix &M, int Cols, int Rows );

void traceM_init( TtraceMatrix &M);

bool needsExtention( std::string &S );

TPos calcEndPos( Tpath& P );

void traceM_construct( TmatrixPaths &paths, generalmatrix &Mat, int mx );

int pointExtention( TtraceMatrix &M, int i, int j );

void traceM_plot( TtraceMatrix &M, Tpath &p );

void traceM_fill( TtraceMatrix &M, TmatrixPaths &paths );
	
void traceM_Draw( TtraceMatrix &M, generalmatrix &Mat );



#endif 