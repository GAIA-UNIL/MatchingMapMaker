#define MATLAB
#define _USE_MATH_DEFINES
#include <cmath>
#include <mex.h>
#include <atomic>
#include <future>
#include <algorithm>
#include <cstring>
#if _OPENMP
	#include <omp.h>
#endif
#include <fftw3.h>

#include "complexMulti.hpp"

#ifdef __cplusplus 
	extern "C" bool utIsInterruptPending();
	extern "C" bool utSetInterruptPending(bool);
#else
	extern bool utIsInterruptPending();
	extern bool utSetInterruptPending(bool);
#endif


#ifdef HBW_MALLOC
	#include <hbwmalloc.h>
	#define mem_malloc hbw_malloc
	#define mem_free hbw_free
#else
	#define mem_malloc malloc
	#define mem_free free
#endif


#include "matrix.h"

void mexFunctionWork(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[],std::atomic<bool> &done){

	const mxArray* pattern=prhs[0];
	const mxArray* im1=prhs[1];
	const mxArray* im2=prhs[2];
	const mxArray* lags=prhs[3];

	float* kernel=(float*)mxGetPr(pattern);
	float* im1Ptr=(float*)mxGetPr(im1);
	float* im2Ptr=(float*)mxGetPr(im2);
	int* lagsPtr=(int*)mxGetPr(lags);

	mwSize patternDim = mxGetNumberOfDimensions(pattern);
	const mwSize *sizePattern = mxGetDimensions(pattern);
	int offsetPX=-sizePattern[0]/2;
	int offsetPY=-sizePattern[1]/2;

	if(patternDim!=2) mexErrMsgIdAndTxt("error", "pattern need to be 2D");

	mwSize imageDim1 = mxGetNumberOfDimensions(im1);
	mwSize imageDim2 = mxGetNumberOfDimensions(im2);
	const mwSize *sizeIm1 = mxGetDimensions(im1);
	const mwSize *sizeIm2 = mxGetDimensions(im2);

	if(imageDim1!=2) mexErrMsgIdAndTxt("error", "only 2D");
	if(imageDim2!=2) mexErrMsgIdAndTxt("error", "only 2D");

	for (int i = 0; i < 2; ++i)
	{
		//printf("%d, %d \n",sizeIm1[i],sizeIm2[i]);
		if(sizeIm1[i]!=sizeIm2[i]) mexErrMsgIdAndTxt("error", "both image need to have the same size");
	}

	mwSize imageDimLags = mxGetNumberOfDimensions(lags);
	const mwSize * sizeImLags = mxGetDimensions(lags);

	if(imageDimLags!=2) mexErrMsgIdAndTxt("error", "lags need to be two column X and Y");
	if(sizeImLags[1]!=2) mexErrMsgIdAndTxt("error", "lags need to be two column X and Y");

	//printf("Number of lags: %d\n ",sizeImLags[0]);
	mexEvalString("drawnow");

	mwSize m=sizeIm1[0];
	mwSize n=sizeIm1[1];

	//printf("alloc output memory\n ");
	mexEvalString("drawnow");

	mxArray *lagsIndex=mxCreateNumericMatrix(m, n, 
  		mxUINT32_CLASS, mxREAL);
	mxArray *quality=mxCreateNumericMatrix(m, n, 
  		mxSINGLE_CLASS, mxREAL);

	int* lagsIndexPtr=(int*)mxGetPr(lagsIndex);
	float* qualityPtr=(float*)mxGetPr(quality);

	//printf("end allocation output memory\n ");
	mexEvalString("drawnow");

	#pragma omp parallel for simd default(none) firstprivate(qualityPtr,m,n)
	for (int i = 0; i < m*n; ++i)
	{
		qualityPtr[i]=std::numeric_limits<float>::infinity();
	}

	#pragma omp parallel for simd default(none) firstprivate(im1Ptr,im2Ptr,m,n)
	for (int i = 0; i < m*n; ++i)
	{
		if(std::isnan(im1Ptr[i]))im1Ptr[i]=0;
		if(std::isnan(im2Ptr[i]))im2Ptr[i]=0;
		
	}

	//printf("end remove nan \n ");
	mexEvalString("drawnow");


	// do operation

	//printf("temp array memory allocation \n ");
	mexEvalString("drawnow");

	float* tmp1=(float*)mem_malloc(sizeof(fftwf_complex)*(m/2+1)*n);
	float* tmp2=(float*)mem_malloc(sizeof(fftwf_complex)*(m/2+1)*n);
	float* patternMem=(float*)malloc(sizeof(fftwf_complex)*(m/2+1)*n);

	//printf("init memory allocation \n ");
	mexEvalString("drawnow");


	#pragma omp parallel for simd default(none) firstprivate(tmp1,m,n)
	for (int i = 0; i < 2*(m/2+1)*n; ++i)
	{
		tmp1[i]=0;
	}

	for (int py = 0; py < sizePattern[1]; ++py)
	{
		for (int px = 0; px < sizePattern[0]; ++px)
		{
			tmp1[(px+offsetPX+m)%m+((py+offsetPY+n)%n)*m]=kernel[px+py*sizePattern[0]];
		}
	}
	//tmp1[m/2+(n/2)*m]=1;
	//printf("fft plan memory allocation \n ");
	mexEvalString("drawnow");

	//tmp1[0]=1;

	fftwf_init_threads();
	#if _OPENMP
	fftw_plan_with_nthreads(omp_get_max_threads());
	#endif
	fftwf_plan forwardPlan=fftwf_plan_dft_r2c_2d(n,m,tmp1,(fftwf_complex*)tmp2,FFTW_ESTIMATE);
	#if _OPENMP
	fftw_plan_with_nthreads(omp_get_max_threads());
	#endif
	fftwf_plan patternPlan=fftwf_plan_dft_r2c_2d(n,m,tmp1,(fftwf_complex*)patternMem,FFTW_ESTIMATE);
	#if _OPENMP
	fftw_plan_with_nthreads(omp_get_max_threads());
	#endif
	fftwf_plan backwardPlan=fftwf_plan_dft_c2r_2d(n,m,(fftwf_complex*)tmp1,tmp2,FFTW_ESTIMATE);
	fftwf_execute(patternPlan);
	fftwf_destroy_plan(patternPlan);

	//printf("start resolution \n ");
	mexEvalString("drawnow");

	for (int localLagId = 0; localLagId < sizeImLags[0]; ++localLagId)
	{
		if(done) break;
		//memset(tmp1,0,sizeof(float)*m*n);
		auto begin = std::chrono::high_resolution_clock::now();
		#pragma omp parallel for simd default(none) firstprivate(tmp1,m,n)
		for (int i = 0; i < 2*(m/2+1)*n; ++i)
		{
			tmp1[i]=0;
		}
		int deltaX=lagsPtr[0*sizeImLags[0]+localLagId];
		int deltaY=lagsPtr[1*sizeImLags[0]+localLagId];
		printf("%d,%d. -- %d /%d \n", deltaX, deltaY, localLagId, sizeImLags[0]);
		mexEvalString("drawnow");

		int absolutLag=deltaY*m+deltaX;

		#pragma omp parallel for simd firstprivate(im1Ptr,im2Ptr,tmp1,absolutLag,sizePattern,offsetPX,offsetPY,kernel,m)
		for (int index = 0; index < m*n; ++index)
		{
			tmp1[index]=std::abs(im1Ptr[index]-im2Ptr[(index+absolutLag)%(m*n)]);
		}

		fftwf_execute(forwardPlan);
		
		#pragma omp parallel for simd firstprivate(tmp1, tmp2)
		for (int i = 0; i < 2*(m/2+1)*n; ++i)
		{
			tmp1[i]=tmp2[i];
		}
		
		//memset(tmp1,0,sizeof(fftwf_complex)*(m/2+1)*n);
		#pragma omp parallel for simd default(none) firstprivate(tmp1,m,n)
		for (int i = 0; i < 2*(m/2+1)*n; ++i)
		{
			tmp1[i]=0;
		}

		g2s::complexAddAlphaxCxD(tmp1, tmp2, patternMem, 1.0f, (m/2+1)*n);

		fftwf_execute(backwardPlan);

		#pragma omp parallel for simd default(none) firstprivate(tmp2, qualityPtr, lagsIndexPtr, localLagId, m, n)
		for (int i = 0; i < m*n; ++i)
		{
			if(tmp2[i]<qualityPtr[i]){
				qualityPtr[i]=tmp2[i];
				lagsIndexPtr[i]=localLagId+1;
			}
		}

		auto end = std::chrono::high_resolution_clock::now();
		double time = 1.0e-6 * std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
		printf("compuattion time: %7.2f s\n", time/1000);
		mexEvalString("drawnow");
	}

	fftwf_destroy_plan(forwardPlan);
	fftwf_destroy_plan(backwardPlan);
	fftwf_cleanup_threads();
	fftwf_cleanup();

	mem_free(tmp1);
	mem_free(tmp2);
	free(patternMem);

	#pragma omp parallel for simd default(none) firstprivate(qualityPtr, m, n)
	for (int i = 0; i < m*n; ++i)
	{
		qualityPtr[i]/=(m*n);
	}

	//finish
	if(nlhs>0) plhs[0]=lagsIndex;
	if(nlhs>1) plhs[1]=quality;

	done=true;
}


void testIfInterupted(std::atomic<bool> &done){
	while (!done){
		std::this_thread::sleep_for(std::chrono::milliseconds(300));
		if(utIsInterruptPending()){
			done=true;
		}
	}
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	std::atomic<bool> done(false);
	auto myFuture = std::async(std::launch::async, testIfInterupted,std::ref(done));
	mexFunctionWork(nlhs, plhs,  nrhs, prhs,std::ref(done));
	myFuture.wait();
	mexEvalString("drawnow");
}