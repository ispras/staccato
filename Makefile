################################################################################
#   Staccato: Disjoint Support Decompositions from BDDs                        #
#   Copyright (C) 2003-2010  University of Michigan                            #
#   http://www.eecs.umich.edu/staccato                                         #
#                                                                              #
#   Contributors include Stephen Plaza and Valeria Bertacco                    #
#                                                                              #
#   This library is free software; you can redistribute it and/or              #
#   modify it under the terms of the GNU Lesser General Public                 #
#   License as published by the Free Software Foundation; either               #
#   version 2.1 of the License, or (at your option) any later version.         #
#                                                                              #
#   This library is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          #
#   Lesser General Public License for more details.                            #
#                                                                              #
#   You should have received a copy of the GNU Lesser General Public           #
#   License along with this library; if not, write to the Free Software        #
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA              #
#   02110-1301  USA                                                            #
#                                                                              #
#   Stephen Plaza <splaza@umich.edu>                                           #
#   Valeria Bertacco <valeria@umich.edu>                                       #
#                                                                              #
#   University of Michigan                                                     #
#   Electrical Engineering and Computer Science Dept.                          #
#   2260 Hayward St.                                                           #
#   Ann Arbor, MI 48109-2121                                                   #
################################################################################

LIB      = lib/libdsd.a
SRC      = $(shell ls src/*.c)
OBJ      = $(SRC:%.c=%.o)

CUDD_INC = -I/path/to/cudd-2.4.2/include


CFLAGS   = -c -O3 -funroll-all-loops $(CUDD_INC)
CC       = g++ 

	
$(LIB): $(OBJ)
	ar rv $(LIB) $(OBJ) 		

clean: 
	rm -rf $(OBJ) $(LIB) sample

.PHONY: test clean

	
#To compile the program as mentioned in STACCATO Code Example	
CUDD_LIB = -L/path/to/cudd-2.4.2/cudd
EPD_LIB = -L/path/to/cudd-2.4.2/epd
ST_LIB = -L/path/to/cudd-2.4.2/st
UTIL_LIB = -L/path/to/cudd-2.4.2/util
MTR_LIB = -L/path/to/cudd-2.4.2/mtr
DSD_LIB = -L/path/to/STACCATO-1.2/lib
	
test:
	$(CC) src/sample.c $(CUDD_LIB) $(MTR_LIB) $(ST_LIB) $(EPD_LIB) $(UTIL_LIB) $(DSD_LIB) \
	-lcudd -lmtr -lst -lepd -lutil -ldsd $(CUDD_INC) -o sample