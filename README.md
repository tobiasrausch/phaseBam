phaseBam installation (using recursive clone)
---------------------------------------------

`git clone --recursive https://github.com/tobiasrausch/phaseBam.git`

`cd phaseBam/`

`make all`

phaseBam installation (using apt-get)
-------------------------------------

`apt-get update`

`apt-get install -y build-essential g++ git cmake zlib1g-dev ant libbz2-dev libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev`

`git clone https://github.com/samtools/htslib.git`

`cd htslib && make && make lib-static && cd ..`

`export BOOST_ROOT=/usr`

`export SEQTK_ROOT=<htslib_path>`

`git clone https://github.com/tobiasrausch/phaseBam.git`

`cd phaseBam/ && touch .htslib .boost && make all && cd ..`


Running phaseBam
----------------

`./src/phaseBam -s HG00731 -v HG00731.phased.bcf --hap1 HG00731.h1.bam --hap2 HG00731.h2.bam HG00731.bam`
