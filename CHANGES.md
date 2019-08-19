# NCBI External Developer Release:


## NCBI NGS Toolkit 2.10.0
**March 18, 2019**

  **ngs-tools**: updated for _sra-tools 2.10.0_


## NCBI NGS Toolkit 2.9.6
**March 18, 2019**

  **ngs-tools**: changed version to match that of _sra-tools_


## NCBI NGS Toolkit 2.9.4
**January 31, 2019**

  **ngs, ngs-tools**: dump-ref-fasta: added an option to skip non-local references  


## NCBI NGS Toolkit 2.9.3
**October 11, 2018**

  **kns**: added possibility to skip server's certificate validation  
  **kns**: expect to receive HTTP status 200 when sending range-request that includes the whole file  
  **vdb**: fixed a bug in accessing pagemap process request for cursors which do not have pagemap thread running  


## NCBI NGS Toolkit 2.9.2
**July 23, 2018**

  **kfg, vfs**: Introduced enhanced handling of download-only NGC files that lack read/decrypt permissions


## NCBI NGS Toolkit 2.9.1
**June 15, 2018**

  **fasterq-dump**: a tool to dump a whole run in fastq by using a simple query engine approach  
  **kfg, vdb-config**: name resolver service now makes use of fcgi  
  **kfg, vfs**: Fixed a bug that prevented decryption of objects encrypted with non-UTF8 text password keys  
  **kns**: Randomly select from multiple proxies in configuration  
  **ngs-tools**: all tools now report their version correctly  


## NGS Tools 2.9.0
**February 23, 2018**

  **build**: Fixed configure allowing to run it on Perl with version >= v5.26 that has "." removed from @INC  
  **build**: added "smoke tests"  
  **build**: now supports "./configure && make" build in ngs-tools/; "make" works in all subdirectories; gmake invokes CMake  
  **build, ngs-tools**: added make targets runtests and slowtests  
  **build, ngs-tools**: binaries are given versioned names, with corresponding symlinks  
  **build, ngs-tools**: library dependencies search is now based on configuration files  
  **build, sra-tools**: "make runtests" now invokes "make all"  
  **kfg**: added searching of configuration files in ../etc/ncbi/ relative to the binaries  
  **kfs**: fix to improve on windows  
  **klib**: Reverted KTimeMakeTime to use UTC  
  **kns**: Accept the same http_proxy specifications as wget  
  **kns**: Added possibility to report server's IP address after network error  
  **kns**: Ignore HTTP headers sent multiple times  
  **kns**: Improved reporting of network errors  
  **kns**: fixed generation of invalid error code in response to dropped connection  
  **ngs-engine**: improved performance when iterating through partially aligned and unaligned reads  
  **ngs-engine**: optimized filtered access to unaligned runs  
  **ngs-tools**: Added optional dependency on 'sra-tools' needed for some tests  
  **ngs-tools**: Created a tool to compute coverage for contigs  
  **ngs-tools**: added build instructions  
  **ngs-tools, sra-tools**: general-loader and pileup-stats have been moved from sra-tools to ngs-tools   
  **sra-search**: added option --fasta for output in FASTA format  
  **sra-search**: added option to display version number  
  **sra-search**: added option to search unaligned and partially aligned fragments only  
  **sra-search**: improved performance in reference-driven mode  
  **sra-search**: various efficiency/readability improvements in the code   
  **vdb**: An assert triggered by a rare condition was fixed  
  **vdb**: new api to estimate pileup-workload based on slice  
  **vdb**: new function to open HTTP-file with an arbitrary page-size  
  **vdb**: progressbar can now be created to output on stderr  
  **vfs**: Name resolving service was updated and switched to protocol version 3.0  


## NGS Tools 2.8.2
**March 6, 2017**

  **build**: Added ability to specify ncbi-vdb/configure --with-magic-prefix. Look for libraries in (lib lib64) when running "configure --with-...-prefix"  
  **build**: configure detects location of ngs libraries  
  **build**: configure was fixed to skip options unrecognized by gcc 4.4.7  
  **build**: created sra-toolkit Debian package  
  **build**: fixed a bug in 'configure' when in could not find source files in repository saved with non-standard name  
  **kns**: SRA tools respect standard set of environment variables for proxy specification  
  **kns**: updated mbedtls library to version 2.4.1  
  **ngs, sra-search**: now supports search on reference  
  **ngs-tools**: updated the NCBI download page to incorporate ngs versions into 3rd party package names  


## NGS Tools 2.8.0
**October 7, 2016**

### HTTPS-ENABLED RELEASE

  **version**: Moved version from 1.0.x to 2.8.0 to be in sync with sra-tools  
  **build, ngs-tools**: Now ngs-tools look for its dependencies using their normal build paths and does not reconfigure them  
  **build, ngs-tools**: Now ngs-tools use CMAKE_INSTALL_PREFIX for installation path  
  **kns**: All tools and libraries now support https  
  **ngs, ngs-tools, ref-variation**: added class ngs-vdb::VdbAlignment, featuring method IsFirst()  
  **ngs-tools**: Fixed Makefiles to keep supporting "./configure; make" build of sra-search, alongside CMake-based build.  
  **test-sra**: test-sra prints version of ncbi-vdb or ngs-sdk dynamic library  


## NGS Tools 1.0.2
**July 12, 2016**

  **ngs, search, sra-search**: sra-search was modified to support multiple threads.  
  **ngs-engine, ngs-tools, sra-tools, vfs**: The "auxiliary" nodes in configuration are now ignored  
  **ngs-engine**: Added support for blob-by-blob access to SEQUENCE table  
  **ngs-engine**: removed a potential memory leak in NGS_CursorMake()  
  **ngs-tools, sra-search**: now uses NGS-VDB extension of NGS  
  **ngs-tools**: NGS-VDB extension of NGS is now in place  
  **search**: now supports multi-threaded search  
  **sra-search**: now supports sorted output  

## NGS Tools 1.0.1
**May 25, 2016**

  **build, ngs-tools**: 'make test ' will work now when called from ngs-tools  
  **build, ngs-tools**: will create all required directories on the first build  
  **build**: MSVS 2013 toolset (12.0) is now supported across all repositories  
  **ngs-engine**: ncbi-ngs engine was updated: fixed a bug that made NGS read iterator return 0 reads on WGS accessions.  
  **ngs-tools**: now uses CMake as the primary build tool  
  **sra-search**: now supports near matches and gapped matches.  
  **vdb, ngs-engine**: Fixed a bound on memory cache that would never flush under certain access modes  

