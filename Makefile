default: all

getdist: ./source/*.*90
	cd ./source && make getdist

cosmomc: ./source/*.*90
	cd ./source && make

cosmomc_debug: ./source/*.*90
	cd ./source && make cosmomc_debug OUTPUT_DIR=Debug

camspec: ./source/*.*90
	cd ./source && make highL=../highL PLANCKLIKE=cliklike_CamSpec

clean:
	cd ./source && make clean

all: ./source/*.*90
	cd ./source && make all

