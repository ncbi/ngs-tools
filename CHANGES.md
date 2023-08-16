# NCBI External Developer Release:


## NCBI NGS Toolkit 3.0.7
**August 29, 2023**

  **cloud, ngs-tools, sra-tools, vdb-config**: fixed use of AWS credentials  
  **kns, ngs, ngs-tools, sra-tools**: fixed a bug that caused failure during accession resolution while reading HTTP stream  


## NCBI NGS Toolkit 3.0.6
**July 10, 2023**

  **ngs-tools**: changed version to match that of _sra-tools_


## NCBI NGS Toolkit 3.0.5
**May 9, 2023**

  **kfg, ngs-tools, sra-tools, vfs**: stopped using old names resolver cgi  
  **ngs-tools**: updated aligns_to to use the latest libraries  


## NCBI NGS Toolkit 3.0.2
**December 12, 2022**

  **ngs-tools**: changed version to match that of _sra-tools_


## NCBI NGS Toolkit 3.0.1
**November 15, 2022**

  **build**: added support for overriding cmake and ctest commands  
  **build**: will use a system-provided libmbedtls, otherwise the copy included in the source code will be used  
  **ngs-tools**: fixed configure script  
  **ref-variation**: added libraries and tools: ngs-vdb, general-writer, ref-variation, sra-search, general-loader, pileup-stats  


## NCBI NGS Toolkit 3.0.0
**February 10, 2022**

  **ngs-tools**: updated for _sra-tools 3.0.0_


## NCBI NGS Toolkit 2.11.2
**October 7, 2021**

  **klib, ngs-tools, sra-tools**: status messages (-v) are printed to stderr rather than stdout  
  **kns, ngs-tools, sra-tools**: old verbose messages now happen at verbosity > 1  
  **ncbi-vdb, ngs-tools, sra-tools, vdb, vfs**: added  support of SRA Lite files with simplified base quality scores  
  **sra-search**: fixed a data race in reference-driven mode  


## NCBI NGS Toolkit 2.11.1
**August 17, 2021**

  **ncbi-vdb, ngs, ngs-tools, sra-tools**: configure prints the version of compiler  
  **sra-download**: sra-download-manager tool has been removed  


## NCBI NGS Toolkit 2.11.0
**March 15, 2021**

  **build, ncbi-vdb, ngs, ngs-tools**: introduced an additional external library, libncbi-ngs  
  **kfg, sra-tools, vfs, ngs-tools**: dropped support of protected repositories  
  **kns, sra-tools, ngs-tools**: fixed formatting of HTTP requests for proxy  
  **ncbi-vdb, ngs, ngs-tools, sra-tools, vdb**: added support for 64-bit ARM (AArch64, Apple Silicon)  
  **sra-search**: no longer issues unneeded warnings in reference mode  


## NCBI NGS Toolkit 2.10.9
**December 16, 2020**

  **build**: added configure option to produce build in output directory relative to sources  
  **build, ngs-tools**: fixed a c++11 build issue for older compilers  
  **kns, sra-tools, vdb**: added a loop to retry failed connections when fetching SRA files  
  **vfs**: allow to find local files when remote repository is disabled  


## NCBI NGS Toolkit 2.10.8
**June 29, 2020**

  **vfs, sra-tools, ngs-tools**: report an error when file was encrypted for a different ngc file  


## NCBI NGS Toolkit 2.10.7
**May 21, 2020**

  **kns, ngs-tools, sra-tools**: added new header to HTTP requests to communicate VDB version 


## NCBI NGS Toolkit 2.10.6
**MAY 18, 2020**

  **align, sra-tools, ngs-tools**: fixed fetching of reference sequences from cloud  
  **kns, sra-tools, ngs-tools**: added new header to HTTP requests to communicate SRA version  
  **kns, sra-tools, ngs-tools**: introduced a additional configurable network retry loop  


## NCBI NGS Toolkit 2.10.5
**April 1, 2020**

  **ncbi-vdb, ngs, ngs-tools, sra-tools**: all Linux builds now use g++ 7.3 (C++11 ABI) 


## NCBI NGS Toolkit 2.10.4
**February 26, 2020**

  **kns, ngs-tools:**: fixed errors when using ngc file


## NCBI NGS Toolkit 2.10.3
**February 18, 2020***

  **ngs-tools**: updated for _sra-tools 2.10.3_


## NCBI NGS Toolkit 2.10.2
**January 15, 2020**

  **ngs-tools**: updated tax tools to the latest version  


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

