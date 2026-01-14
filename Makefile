-include make.def.*

COMPILER = mpicc
GENERAL_OPT_FLAGS = -std=c11 -O3
DEBUG_FLAGS = -Wall -Wextra
LIB_FLAGS = -lhdf5 -lz -lm

SRC = $(wildcard main.c src/*.c)

compile:
	$(COMPILER) $(SRC) -o lsiblbm $(GENERAL_OPT_FLAGS) $(ADDITIONAL_OPT_FLAGS) $(DEBUG_FLAGS) $(HDF5_FLAGS) $(LIB_FLAGS)

run_local: 
	mpirun -n $(n) lsiblbm

clean:
	rm -f lsiblbm

cleandata:
	rm -f *.h5

install_hdf5:
	mkdir -p hdf5
	tar -xvzf archives/hdf5* -C hdf5 --strip-components=2
	cd hdf5;\
	export ac_cv_lib_sz_SZ_BufftoBuffCompress=no;\
	export ac_cv_header_szlib_h=no;\
	CC=mpicc ./configure --enable-parallel --disable-shared --disable-szlib;\
	make;\
	make install;\
	cd ../;\
	mv hdf5/hdf5/include hdf5/;\
	mv hdf5/hdf5/lib hdf5/
