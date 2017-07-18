# the installation prefix where it will be installed
LIB_PREFIX	= ~/.local/lib/haskell
EXE_PREFIX	= ~/.local/bin

# the Haskell compiler
HC	= ghc

# the Haskell compiler flags
HFLAGS		= -O2 -Wall
HINCLUDE	= src/modules

# building routines
.PHONY: libs
all: libs

libs: data.chemistry.xyz

data.chemistry.xyz:
	cd src/modules/Data/Chemistry && $(HC) -i$(HINCLUDE) $(HFLAGS) --make XYZ.hs 

# cleaning
clean: clean_cabal clean_data.chemistry.xyz

clean_cabal:
	rm -rf dist

clean_data.chemistry.xyz:
	cd src/modules/Data/Chemistry && rm -f *.hi *.o

# installation
install: install_data.chemistry.xyz

install_data.chemistry.xyz:
	mkdir -p $(LIB_PREFIX)/Data/Chemistry
	cd src/modules/Data/Chemistry && cp XYZ.o $(LIB_PREFIX)/Data/Chemistry
