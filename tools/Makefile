all: h52ascii h5clustersize

h52ascii: h52ascii.c
	gcc -std=gnu99 -o h52ascii h52ascii.c -ldfh5 -lhdf5_hl 

h5clustersize: h5clustersize.c
	gcc -std=gnu99 -o h5clustersize h5clustersize.c -lhdf5 -lhdf5_hl -lm
