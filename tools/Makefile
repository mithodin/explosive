.PHONY: all clean

all: h52ascii h5clustersize

h52ascii: h52ascii.c
	gcc -std=gnu99 -o h52ascii h52ascii.c -lhdf5 -lhdf5_hl -O4

h5clustersize: h5clustersize.c
	gcc -std=gnu99 -o h5clustersize h5clustersize.c -lhdf5 -lhdf5_hl -lm -O4

clean:
	rm h52ascii h5clustersize
