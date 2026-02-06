#ROOT
root_config = root-config
root_include = $(shell $(root_config) --cflags)
root_libs = $(shell $(root_config) --libs) -lMinuit -lEG -lSpectrum -lGeom -lEve -lGui
root_ldflag = $(shell $(root_config) --ldflags)

# Directories
bin_dir = bin
src_dir = src
dst_dir = dst
lib_dir = lib
link_dir = link
include_dir = include
src_build_dir = $(src_dir)/build
dst_build_dir = $(dst_dir)/build
dict_dir = $(dst_dir)/dict

CXX = g++
CFLAGS	= -g -O3 -Wall
CFLAGS	+= -Wno-sign-compare
CFLAGS	+= -std=c++17
DFLAGS += -DMemoryLeak
FLAGS = $(CFLAGS) $(DFLAGS) -I. $(root_include) $(addprefix -I,$(src_dir) $(include_dir))


lib_objs = $(lib_dir)/libAnalyzer.a
analib = -L$(lib_dir)/ -lAnalyzer
libs = $(root_libs)
#
core_srcs = $(wildcard $(src_dir)/*.cc)
headers = $(wildcard $(include_dir)/*.hh)
dst_srcs = $(wildcard $(dst_dir)/*.cc)

core_deps = $(core_srcs:$(src_dir)/%.cc=$(src_build_dir)/%.d)
dst_deps = $(dst_srcs:$(dst_dir)/%.cc=$(dst_build_dir)/%.d)

core_objs	= $(core_srcs:$(src_dir)/%.cc=$(src_build_dir)/%.o)
dst_objs	= $(dst_srcs:$(dst_dir)/%.cc=$(dst_build_dir)/%.o)

dst_target	= $(dst_srcs:$(dst_dir)/%.cc=$(bin_dir)/%)



rootcint = rootcint
rootdict_src = $(dict_dir)/DstEventDisplay_E42Dict.cpp 
rootdict_obj = $(dict_dir)/DstEventDisplay_E42Dict.o 

LinkDefFile = $(link_dir)/LinkDef.h  
src_files = $(dst_dir)/DstEventDisplay_E42.h 

ar	:= ar r
echo	:= /bin/echo -e
mkdir	:= mkdir -p
mv	:= mv -f
rm	:= rm -rfv
cp      := cp -fr


.PHONY: all dst lib clean distclean

all: lib dst

-include $(core_deps) $(dst_deps)
lib: $(lib_objs) $(core_objs)
dst: $(dst_target) $(dst_objs) $(rootdict_obj)


#Additional Link for EventDisplay
$(bin_dir)/DstEventDisplay_E42: $(dst_build_dir)/DstEventDisplay_E42.o $(rootdict_obj)
	@ $(mkdir) $(bin_dir)
	$(CXX) $(FLAGS) -o $@ $< $(rootdict_obj) $(analib) $(libs)

$(bin_dir)/%: $(dst_build_dir)/%.o $(lib_objs)
	@ $(mkdir) $(bin_dir)
	$(CXX) $(FLAGS) -o $@ $< $(analib) $(libs)

$(lib_objs): $(core_objs)
	@ $(mkdir) $(lib_dir)
	$(ar) $(lib_objs) $(core_objs)

$(src_build_dir)/%.o: $(src_dir)/%.cc
	@ $(mkdir) $(src_build_dir)
	$(CXX) $(FLAGS) -o $@ -MMD -c $<

$(dst_build_dir)/%.o: $(dst_dir)/%.cc
	@ $(mkdir) $(dst_build_dir)
	$(CXX) $(FLAGS) -o $@ -MMD -c $<



$(dict_dir)/DstEventDisplay_E42Dict.cpp: $(LinkDefFile) $(dst_dir)/DstEventDisplay_E42.h $(headers)
	@ $(mkdir) -p $(dict_dir)
	$(rootcint) -q -f $(dict_dir)/DstEventDisplay_E42Dict.cpp -c $(dst_dir)/DstEventDisplay_E42.h $(headers) $(LinkDefFile) -I$(include_dir)

$(rootdict_obj): $(dict_dir)/DstEventDisplay_E42Dict.cpp
	$(CXX) $(FLAGS) -o $(rootdict_obj) -c $(dict_dir)/DstEventDisplay_E42Dict.cpp

clean:
	@ $(rm) $(dst_objs) $(core_objs) $(lib_objs)
	@ $(rm) $(dst_deps) $(core_deps)
	@ $(rm) -rf $(dict_dir)/* 


distclean:
	@ $(rm) $(bin_dir)/* $(lib_dir)/* $(src_build_dir) $(dst_build_dir)
	@ $(rm) $(dict_dir)




