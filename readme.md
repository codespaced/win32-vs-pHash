win32-vs-pHash
==============

This is a (slightly) modified version of [pHash](http://www.phash.org/) and MVPTree compilable under Visual Studio (tested with Express 2012).  

All of the static libraries in lib were built against MSVCR110D.DLL which comes with VS.

After compilation
-----------------
###Imget

imget needs pHash.dll to be on your path (or in the same directory).

>
>imget command filename [directory] [radius] [k-nearest]
>
>  command    - command - e.g. 'add', 'query' or 'print'
>
>  filename   - file from which to read the tree
>
>  directory  - directory (for add and query)
>
>  radius     - radius for query operation (default = 21.0)
>
>  k-nearest  - max number of results (default = 5)
>

####Running imget.exe
    
    >imget add images.mvp ..\picturedir
 
    >imget print images.mvp
 
    >images query images.mvp ..\querydir 10.1

The default radius (21.0) is normally sufficient. Restricitng the radius further (closer to 0) will reduce the results you get, but expanding it (with small picture sets atleast) will not increase your results past a certain point.





##3rd Party Packages and licenses:

*[pHash 0.9.6 and MVPTree 1.0.0](http://www.phash.org/) - The open source perceptual hash library
 [GNU LGPL v3](http://www.gnu.org/licenses/gpl-3.0.html)


*[mman-win32](http://code.google.com/p/mman-win32/) - mman library for Windows
 [GNU GPL v2](http://www.gnu.org/licenses/gpl-2.0.html)


*[dirent-1.13](http://softagalleria.net/download/dirent/dirent-1.13.zip) - dirent API for Microsoft Visual Studio
 [custom, see dirent.h](shared/dirent.h)


*[CImg 1.5.7](http://cimg.sourceforge.net) - The C++ Template Image Processing Toolkit.
 [CeCILL-C](http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) or [CeCILL v2.0](http://www.cecill.info/licences/Licence_CeCILL_V2-en.html)


*[JPEG 9](http://www.ijg.org/) - The Independent JPEG Group (IJG)
 [Independent JPEG Group License](http://directory.fsf.org/wiki?title=License:JPEG)


*[libpng 1.6.3](http://www.libpng.org/pub/png/libpng.html) - The official PNG reference library
 [Open Source](http://www.libpng.org/pub/png/src/libpng-LICENSE.txt)


*[zlib 1.2.8](http://www.zlib.net/) - Massively Spiffy Yet Delicately Unobtrusive Compression Library
 [zlib license](http://www.zlib.net/zlib_license.html)


*[Pthreads-win32](http://www.sourceware.org/pthreads-win32/) - POSIX Threads Library for Win32
 [GNU LGPL v2.1](http://www.gnu.org/licenses/lgpl-2.1.html) [Details](http://www.sourceware.org/pthreads-win32/copying.html)

