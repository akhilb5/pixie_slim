#!/bin/make
# GNUmakefile using implicit rules and standard definitions
SHELL=/bin/sh

# uncomment this line if processing Rev. D data
REVISIOND = 1
#------- instruct make to search through these
#------- directories to find files
vpath %.f scan/ 
vpath %.h include/
vpath %.cpp src/
vpath %.cxx src/

DIRA2=/opt/local/hhirf
#DIRB= /usr/acq2/lib
#DIRA2 = /usr/hhirf
DIRB  = /opt/local/hhirf
#DESTDIR=/home/pixie16/MTAS/RUN01
DESTDIR=/Users/akhilbhardwaj/I136_analysis/pixie_slim/bin

LIBS = $(DIRA2)/scanorlib.a $(DIRA2)/orphlib.a \
       $(DIRB)/acqlib.a  $(DIRB)/ipclib.a

OutPutOpt     = -o $(DESTDIR)/pixie_ldf_c_slim # keep whitespace after "-o"
ObjSuf        = o

#------- define file suffixes
fSrcSuf   = f
cSrcSuf   = c
c++SrcSuf = cpp
cxxSrcSuf = cxx

#------- define compilers
#FC        = g77 
#uncomment to compile with gfortran (>=4.2) if required for the hhirf libs
FC        = gfortran
GCC       = /opt/homebrew/bin/gcc-11
#GCC       = gcc
CXX       = /opt/homebrew/bin/g++-11
#CXX       = g++
LINK.o    = $(FC) $(LDFLAGS) $(TARGET_ARCH)

#------- define basic compiler flags, no warnings on code that is not my own
FFLAGS   += -g

GCCFLAGS += -g -fPIC $(CINCLUDEDIRS)
CXXFLAGS += -g -Wall -fPIC $(CINCLUDEDIRS) -Dnewreadout 
ifdef REVISIOND
CXXFLAGS += -DREVD
endif

#------- include directories for the pixie c files
CINCLUDEDIRS  = -Iinclude

#------- basic linking instructions
LDLIBS   += -lm -lstdc++ -lgcc 

ifdef BLINDED
CXXFLAGS += -DBLINDED
endif

ifeq ($(FC),gfortran)
FFLAGS	+= -fsecond-underscore
LDLIBS	+= -lgfortran
GCCFLAGS += -O3
CXXFLAGS += -O3 -std=c++11
else
LDFLAGS += -g77libs  
LDLIBS	+= -lg2c
endif

#-------- define file variables -----------------------

# objects from fortran
SET2CCO          = set2cc.$(ObjSuf)
MESSLOGO         = messlog.$(ObjSuf)
MILDATIMO        = mildatim.$(ObjSuf)
SCANORUXO        = scanorux.$(ObjSuf)

# objects from cpp
PIXIEO           = PixieStd.$(ObjSuf)
ACCUMULATORO     = StatsAccumulator.$(ObjSuf)
HISTOGRAMMERO    = DeclareHistogram.$(ObjSuf)
EVENTPROCESSORO  = EventProcessor.$(ObjSuf)
SCINTPROCESSORO  = ScintProcessor.$(ObjSuf)
GEPROCESSORO     = GeProcessor.$(ObjSuf)
MCPPROCESSORO    = McpProcessor.$(ObjSuf)
MTCPROCESSORO    = MtcProcessor.$(ObjSuf)
DSSDPROCESSORO   = DssdProcessor.$(ObjSuf)
SSDPROCESSORO    = SsdProcessor.$(ObjSuf)
MTASPROCESSORO    = MtasProcessor.$(ObjSuf)
TRACESUBO        = TraceAnalyzer.$(ObjSuf)
DETECTORDRIVERO  = DetectorDriver.$(ObjSuf)
CORRELATORO      = Correlator.$(ObjSuf)
RAWEVENTO        = RawEvent.$(ObjSuf)
ROOTPROCESSORO   = RootProcessor.$(ObjSuf)
RANDOMPOOLO      = RandomPool.$(ObjSuf)
STATSDATAO       = StatsData.$(ObjSuf)
WAVEFORMPROCESSORO = WaveformProcessor.$(ObjSuf)
PULSERPROCESSORO   = PulserProcessor.$(ObjSuf)
PSPMTPROCESSORO   = PspmtProcessor.$(ObjSuf)
MTASPSPMTPROCESSORO   = MtasPspmtProcessor.$(ObjSuf)
#VANDLEPROCESSORO   = VandleProcessor.$(ObjSuf)


ifdef USEROOT
PIXIE            = pixie_ldf_c_root$(ExeSuf)
else
PIXIE            = pixie_ldf_c$(ExeSuf)
endif

ifdef REVISIOND
READBUFFDATAO    = ReadBuffData.RevD.$(ObjSuf)
else
READBUFFDATAO    = ReadBuffData.$(ObjSuf)
endif

#----- list of objects
OBJS   = $(READBUFFDATAO) $(SET2CCO) $(DSSDSUBO) $(DETECTORDRIVERO) \
	$(MTCPROCESSORO) $(MCPPROCESSORO) $(CORRELATORO) \
	$(MESSLOGO) $(MILDATIMO) $(SCANORUXO) $(ACCUMULATORO) $(PIXIEO) \
	$(HISTOGRAMMERO) $(EVENTPROCESSORO) $(SCINTPROCESSORO) \
	$(GEPROCESSORO) $(SPLINEFITPROCESSORO) $(SPLINEPROCESSORO) \
	$(DSSDPROCESSORO) $(SSDPROCESSORO) $(RAWEVENTO) $(RANDOMPOOLO) \
	$(MTASPROCESSORO) $(STATSDATAO) \
	$(WAVEFORMPROCESSORO)  $(PULSERPROCESSORO) \
	$(TRACESUBO) $(PSPMTPROCESSORO) $(MTASPSPMTPROCESSORO)
#$(VANDLEPROCESSORO) $(PULSERPROCESSORO) \


ifdef USEROOT
OBJS  += $(ROOTPROCESSORO)
endif

PROGRAMS = $(PIXIE)

DISTTARGETS = src include scan manual Makefile Doxyfile map.txt cal.txt
DISTNAME = pixie_scan
DOCSTARGETS = html latex

#--------- Add to list of known file suffixes
.SUFFIXES: .$(cxxSrcSuf) .$(fSrcSuf) .$(c++SrcSuf) .$(cSrcSuf)

.phony: all clean dist distdocs
all:     $(PROGRAMS)

#----------- remove all objects, core and .so file
clean:
	@echo "Cleaning up..."
	@rm -f $(OBJS) $(PIXIE) core *~ src/*~ include/*~ scan/*~

dist:
	@mkdir $(DISTNAME)
	@cp -t $(DISTNAME) -r $(DISTTARGETS)
	@tar -czf $(DISTNAME)-`date +%d%m%y`.tgz --exclude=*~ $(DISTNAME)
	@$(RM) -r $(DISTNAME)

distdocs:
	@mkdir $(DISTNAME)
	@cp -t $(DISTNAME) -r $(DISTTARGETS) $(DOCSTARGETS)
	@tar -czf $(DISTNAME)-with-docs-`date +%d%m%y`.tgz --exclude=*~ $(DISTNAME)
	@$(RM) -r $(DISTNAME)


#----------- link all created objects together
#----------- to create pixie_ldf_c program

ifeq ($(FC),gfortran)
$(PIXIE): $(OBJS) $(LIBS)
	$(LINK.o) $^ -o $(DESTDIR)/$@_slim  $(LDLIBS)
else
$(PIXIE): $(OBJS) $(LIBS)
	$(LINK.o) $(LDLIBS) $^ -o $(DESTDIR)/$@_slim 
endif
