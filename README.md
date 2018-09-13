# LZKP

## Install

* Install [GSL](https://packages.ubuntu.com/xenial/libgsl-dev) library:
  * `sudo apt-get install libgsl-dev`
* Install [GMP](https://packages.ubuntu.com/xenial/libgmp-dev) library:
  * `sudo apt-get install libgmp-dev`
* Install [cryptoTools](https://github.com/ladnir/cryptoTools) library:
  * `cd ~`
  * `git clone --recursive https://github.com/ladnir/cryptoTools`
  * `cd cryptoTools/thirdparty/linux`
  * Enable program_options package for boost:
    * `vi boost.get`
    * Change line:
    
      `./b2 stage --with-system --with-thread link=static -mt`
      
      To:
      
      `./b2 stage --with-system --with-thread --with-program_options link=static -mt`
  * `bash all.get`
  * `cd ../..`
  * `cmake  .`
  * `make`
  * `sudo make install`
* Configure boost:
  * `echo export BOOST_ROOT=~/cryptoTools/thirdparty/linux/boost >> ~/.bashrc`
  * `source ~/.bashrc`
* Install LZKP:
  * `cd ~`
  * `git clone git@github.com:cryptobiu/LZKP.git`
  * `cd LZKP`
  * `cmake .`
  * `make`
  
  
