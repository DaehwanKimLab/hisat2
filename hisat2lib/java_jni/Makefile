#
# Copyright 2018, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#

JAVA_HOME := $(JAVA_HOME)

HISAT2_DIR = ../..

CC = gcc


MACOS = 0
ifneq (,$(findstring Darwin,$(shell uname)))
MACOS = 1
endif

INCLUDES = -I $(JAVA_HOME)/include
INCLUDES += -I $(JAVA_HOME)/include/linux

TARGET_LIB = libht2jni.so
SHARED_FLAG = -shared

ifeq (1,$(MACOS))
INCLUDES += -I $(JAVA_HOME)/include/darwin
TARGET_LIB = libht2jni.jnilib
SHARED_FLAG = -dynamiclib
endif

INCLUDES += -I $(HISAT2_DIR)/hisat2lib

CFLAGS = -fPIC $(INCLUDES)
#CFLAGS += -DDEBUG

SRCS = ht2module.c

OBJS = $(SRCS:.c=.o)

all: lib 

lib: $(TARGET_LIB)

$(TARGET_LIB): HT2Module.h $(OBJS)
	$(CC) $(SHARED_FLAG) -o $@ $(OBJS) $(HISAT2_DIR)/libhisat2lib.a -lstdc++

HT2Module.h HT2Module.class: HT2Module.java
	javac -h . HT2Module.java

example: HT2ModuleExample.class

HT2ModuleExample.class: HT2ModuleExample.java
	javac HT2ModuleExample.java

clean:
	rm -f *.class HT2Module.h HT2Module_HT2Position.h *.so *.o $(TARGET_LIB)

test: all example
	java -Djava.library.path=. HT2ModuleExample
