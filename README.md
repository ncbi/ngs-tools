ngs-tools
===

## Quick instructions for building and installing **ngs-tools** from source.

1. Create an **ncbi** source directory if needed.
While this is not strictly required, our configuration scripts will benefit by
being able to locate related projects without asking for explicit paths:

`$ mkdir ncbi`

2. Check out the sources:

`$ cd ncbi`

`$ git clone https://github.com/ncbi/ngs-tools.git`

3. If you have not yet installed **ncbi-vdb**, please do so now:

`$ git clone https://github.com/ncbi/ncbi-vdb.git`

and follow directions at
[https://github.com/ncbi/ncbi-vdb/wiki/Building-and-Installing-from-Source]
(https://github.com/ncbi/ncbi-vdb/wiki/Building-and-Installing-from-Source)
.

4. Configure the build:

`$ cd ngs-tools`

`$ ./configure`

By default, configure will select a build-output directory under your `$HOME`
and will install under `/usr/local/ncbi/ngs-tools` on Linux.
The default settings can be changed, of course. For all options, you can run:

`$ ./configure --help`

5. Make the tools:

`$ make`

6. Install the tools as admin (you may be asked for a password):

`$ sudo make install`

At this point, the installation should be complete,
although you will probably have to login again before all changes take place.
If the installation is successful,
you should find executables installed and an update to shell variables
_(only AFTER logging in again)_. To verify update of your environment:

`$ echo $PATH`  # _should now have the path to your installed ngs-tools, and_

`$ which sra-search`  # _should return the location of this utility_.
