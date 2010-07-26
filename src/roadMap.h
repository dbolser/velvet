/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#ifndef _ROADMAP_H_
#define _ROADMAP_H_

#include <stdio.h>

struct roadMapArray_st {
	RoadMap *array;
	Annotation *annotations;
	IDnum length;
	int WORDLENGTH;
	boolean double_strand;
	IDnum referenceCount;
};

////////////////////////////////////////////////////////////////////
//      Annotation stuff
////////////////////////////////////////////////////////////////////
IDnum getAnnotSequenceID(Annotation * annot);
Coordinate getFinish(Annotation * annot);
Coordinate getStart(Annotation * annot);
Coordinate getPosition(Annotation * annot);
Coordinate getAnnotationLength(Annotation * annot);
void incrementAnnotationCoordinates(Annotation * annot);

void setStartID(Annotation * annot, IDnum nodeID);
IDnum getStartID(Annotation * annot);
void setFinishID(Annotation * annot, IDnum nodeID);
IDnum getFinishID(Annotation * annot);

char *readAnnotation(Annotation * annot);

Annotation *getNextAnnotation(Annotation * annot);

////////////////////////////////////////////////////////////////////
//      RoadMap stuff
////////////////////////////////////////////////////////////////////
RoadMap *newRoadMap();

IDnum getAnnotationCount(RoadMap * rdmap);

RoadMap *getRoadMapInArray(RoadMapArray * array, IDnum index);

// Same thing but for the RoadMap file generated by the hash 
RoadMapArray *importRoadMapArray(char *filename);
void destroyRoadMapArray(RoadMapArray * rdmap); 
#endif
