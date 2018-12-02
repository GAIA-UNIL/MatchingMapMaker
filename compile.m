if ismac
    mex movsae2.cpp -I/opt/local/include/ -L/opt/local/lib -lfftw3f -lfftw3f_threads -lut
elseif isunix
    mex movsae2.cpp -I/usr/include/ -L/usr/lib -lfftw3f -lfftw3f_threads -lut
elseif ispc
    if(exist(fftw3')==0)
        mget(ftp('ftp.fftw.org'),'/pub/fftw/fftw-3.3.5-dll64.zip');
        unzip('pub/fftw/fftw-3.3.5-dll64.zip','fftw3');
    end 
    mex('movsae2.cpp','-Iffw3','-Lffw3','-lfftw3f-3','-lut');
else
    disp('Platform not supported')
end

% compilation with icpc

%icpc -std=c++11 -mkl -xhost -fopenmp -c  -I/opt/intel/compilers_and_libraries_2018.2.199/linux/mkl/include/fftw/ -I"/home/knl/nvme/MATLAB//R2018a/extern/include"  -DHBW_MALLOC  -fPIC -std=c++11 xcorr2_4Lags.cpp -o xcorr2_4Lags.o && icpc -std=c++11 -mkl -xhost -fopenmp  -shared xcorr2_4Lags.o  -lut -L/home/knl/nvme/MATLAB/R2018a/bin/glnxa64 -lmx -lmex -lmat -lm -Wl,-rpath-link,/usr/local/MATLAB/R2018a/bin/glnxa64 -lstdc++ -lirng -lirc -lintlc -lmemkind -o xcorr2_4Lags.mexa64
%export LD_PRELOAD=/opt/intel/compilers_and_libraries_2018.2.199/linux/compiler/lib/intel64/libirc.so:/usr/lib64/libmemkind.so:/opt/intel/compilers_and_libraries_2018.2.199/linux/mkl/lib/intel64_lin/libmkl_core.so:/opt/intel/compilers_and_libraries_2018.2.199/linux/mkl/lib/intel64_lin/libmkl_intel_thread.so
%matlab -nodesktop -nosplash -nodisplay
