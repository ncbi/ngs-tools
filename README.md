ngs-tools
===

### IMPORTANT ANNOUNCEMENT
As was announced in [OMB memorandum M-15-13](https://www.whitehouse.gov/sites/default/files/omb/memoranda/2015/m-15-13.pdf) on June 8, 2015, NCBI and all Federal systems will be transitioning to using HTTPS-only protocols before the end of 2016. This change will affect any software that uses NCBI APIs such as the E-utilities or NCBI software toolkits such as `sra-tools`, `ncbi-vdb` or `ngs`.

The NLM and NCBI may implement the switch to HTTPS-only as early as September 30, 2016.

In particular, software products that depend on `sra-tools`, `ncbi-vdb` or `ngs` may not function as expected after September 30 unless they are properly updated from this site or by the software provider.

If you use software that accesses NCBI SRA data in any way, your software will likely be affected by this change. Please check with your software provider for recent udpates or patches, and be sure to acquire these before September 30.
 
If you develop software that relies on `sra-tools`, `ncbi-vdb` or `ngs` in any way, you will likely need to update your code so that it accesses NCBI using HTTPS.

We have released new tools with version 2.8.0 that are HTTPS compatible and `M-15-13` compliant as of October 7, 2016. Please be certain to [update all of your binaries](https://github.com/ncbi/sra-tools/wiki/Downloads) and configuration files.

#### Example installation of E-Utilities and NCBI software toolkits:

`
mkdir ~/SRC ; cd ~/SRC
git clone https://github.com/ncbi/ngs.git
git clone https://github.com/ncbi/ncbi-vdb.git
git clone https://github.com/ncbi/ngs-tools.git
cd ngs/ngs-sdk ; ./configure -p=~/NCBI && make
cd ../../ncbi-vdb ; ./configure -p=~/NCBI && make
cd ../ngs-tools ; ./configure -p=~/NCBI
`
