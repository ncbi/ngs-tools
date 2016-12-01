# NCBI External Developer Release:

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

