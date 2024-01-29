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

# Build static library by default
BUILD_TYPE ?= static

ifeq ($(BUILD_TYPE), static)
	LIB = lib/libdsd.a
else
	LIB = libSTACCATO.so
endif

SRC = $(shell ls src/*.c)
OBJ = $(SRC:%.c=%.o)

CFLAGS	= -c -O3 -funroll-all-loops -fPIC -I$(CUDD_INCLUDE)/util \
-I$(CUDD_INCLUDE)/cudd -I$(CUDD_INCLUDE) -I$(CUDD_INCLUDE)/st \
-I$(CUDD_INCLUDE)/mtr -I$(CUDD_INCLUDE)/epd

CFLAGS += $(SM)
CFLAGS += $(GC)

# Determining the path to the libcudd.so library
ifeq ($(CUDD_DIR),)
	CUDD_LIBRARY_PATH = /usr/local/lib
else
	CUDD_LIBRARY_PATH = $(CUDD_DIR)
endif

CC = gcc

ifeq ($(BUILD_TYPE), static)
$(LIB): $(OBJ)
	@mkdir -p lib
	ar rv $(LIB) $(OBJ)
else
$(LIB): $(OBJ) $(CUDD_LIBRARY_PATH)/libcudd.so
	$(CC) -shared -o $(LIB) $(OBJ) -L$(CUDD_LIBRARY_PATH) -lcudd
endif

clean:
	rm -rf $(OBJ) $(LIB) libSTACCATO.so lib sample

.PHONY: test clean

# To compile the program as mentioned in STACCATO Code Example
DSD_LIB = -L$(shell pwd)/lib
    
test:
	$(CC) src/sample.c $(DSD_LIB) -lcudd -ldsd \
	-I$(CUDD_INCLUDE)/util -I$(CUDD_INCLUDE)/cudd -I$(CUDD_INCLUDE) \
	-I$(CUDD_INCLUDE)/st -I$(CUDD_INCLUDE)/mtr -I$(CUDD_INCLUDE)/epd -o sample
