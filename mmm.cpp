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
#include "mexInterrupt.hpp"


#ifdef HBW_MALLOC
	#include <hbwmalloc.h>
	#define mem_malloc hbw_malloc
	#define mem_free hbw_free
#else
	#define mem_malloc malloc
	#define mem_free free
#endif


#include "matrix.h"

enum matchingMode {
	SAE,
	MAE,
	SSE,
	MSE,
	CC,
	NMAE,
	NMSE,
	NCC
};

inline float sqr(float x){return x*x;}

void mexFunctionWork(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[],std::atomic<bool> &done){

	matchingMode mode=SAE;
	bool needNormalized=false;
	bool skipPatternNomailize=false;
	bool addAbsolut=false;
	bool lookForMax=false;
	bool needDeviation=false; 
	if(nrhs<4){
		done=true;
		mexErrMsgTxt("not enough argument\n(pattern, image1, image2, lagsVector)\n(pattern, image1, image2, lagsVector, algo)");
	}

	const mxArray* pattern=prhs[0];
	const mxArray* im1=prhs[1];
	const mxArray* im2=prhs[2];
	const mxArray* lags=prhs[3];
	if((nrhs>4) && mxIsChar(prhs[4])){
		char algo[4];
		if(!mxGetString(prhs[4], algo, 4))
		{
			if(!strcmp("SAE",algo))mode=SAE;
			if(!strcmp("MAE",algo))mode=MAE;
			if(!strcmp("SSE",algo))mode=SSE;
			if(!strcmp("MSE",algo))mode=MSE;
			if(!strcmp("CC",algo))mode=CC;
			if(!strcmp("NMAE",algo))mode=NMAE;
			if(!strcmp("NMSE",algo))mode=NMSE;
			if(!strcmp("NCC",algo))mode=NCC;
		}
	}

	if(!mxIsSingle(pattern) || !mxIsSingle(im1) || !mxIsSingle(im2) || !mxIsInt32(lags))
	{
		done=true;
		mexErrMsgTxt("pattren, image1 and image2 need to be single and lagVector int32");
	}

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

	switch(mode){
		case SAE:
		case SSE:
			skipPatternNomailize=true;
			break;
		case CC:
			//addAbsolut=true;
			lookForMax=true;
			break;
		case NMAE:
		case NMSE:
			needNormalized=true;
			break;
		case NCC:
			needNormalized=true;
			//addAbsolut=true;
			lookForMax=true;
			needDeviation=true;
			break;
		default:
			break;
	}

	float initValue=0;

	if(lookForMax)
		initValue=-std::numeric_limits<float>::infinity();
	else
		initValue=std::numeric_limits<float>::infinity();

	#pragma omp parallel for simd default(none) firstprivate(initValue,qualityPtr,m,n)
	for (int i = 0; i < m*n; ++i)
	{
		qualityPtr[i]=initValue;
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
	float* im1Mean=nullptr;
	float* im2Mean=nullptr;
	float *stdDev1=nullptr;
	float *stdDev2=nullptr;

	if(needNormalized || needDeviation){
		im1Mean=(float*)malloc(sizeof(fftwf_complex)*(m/2+1)*n);
		im2Mean=(float*)malloc(sizeof(fftwf_complex)*(m/2+1)*n);
	}
	if(needDeviation){
		stdDev1=(float*)malloc(sizeof(fftwf_complex)*(m/2+1)*n);
		stdDev2=(float*)malloc(sizeof(fftwf_complex)*(m/2+1)*n);
	}

	//printf("init memory allocation \n ");
	mexEvalString("drawnow");

	#pragma omp parallel for simd default(none) firstprivate(tmp1,m,n)
	for (int i = 0; i < 2*(m/2+1)*n; ++i)
	{
		tmp1[i]=0;
	}

	float sumKernel=1.0f;
	if(!skipPatternNomailize)
	{
		sumKernel=0.0f;
		for (int py = 0; py < sizePattern[1]; ++py)
		{
			for (int px = 0; px < sizePattern[0]; ++px)
			{
				sumKernel+=kernel[px+py*sizePattern[0]];
			}
		}
	}

	for (int py = 0; py < sizePattern[1]; ++py)
	{
		for (int px = 0; px < sizePattern[0]; ++px)
		{
			tmp1[(px+offsetPX+m)%m+((py+offsetPY+n)%n)*m]=kernel[px+py*sizePattern[0]]/sumKernel;
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
	fftwf_plan backwardPlan=fftwf_plan_dft_c2r_2d(n,m,(fftwf_complex*)tmp1,tmp2,FFTW_ESTIMATE);
	
	fftwf_execute_dft_r2c(forwardPlan,tmp1,(fftwf_complex*)patternMem);

	//compute mean
	if(needNormalized || needDeviation){
		memset(tmp2,0,2*(m/2+1)*n*4);
		fftwf_execute_dft_r2c(forwardPlan,im1Ptr,(fftwf_complex*)tmp1);
		g2s::complexAddAlphaxCxD(tmp2, tmp1, patternMem, 1.0f, (m/2+1)*n);
		fftwf_execute_dft_c2r(backwardPlan,(fftwf_complex*)tmp2,im1Mean);
		memset(tmp2,0,2*(m/2+1)*n*4);
		fftwf_execute_dft_r2c(forwardPlan,im1Ptr,(fftwf_complex*)tmp1);
		g2s::complexAddAlphaxCxD(tmp2, tmp1, patternMem, 1.0f, (m/2+1)*n);
		fftwf_execute_dft_c2r(backwardPlan,(fftwf_complex*)tmp2,im2Mean);
		#pragma omp parallel for simd default(none) firstprivate(im1Mean,im2Mean,m,n)
		for (int i = 0; i < m*n; ++i)
		{
			im1Mean[i]/=m*n;
			im2Mean[i]/=m*n;
		}
	}

	if(needDeviation){
		#pragma omp parallel for simd default(none) firstprivate(tmp2,im1Ptr,m,n)
		for (int i = 0; i < m*n; ++i)
		{
			tmp2[i]=sqr(im1Ptr[i]);
		}
		fftwf_execute_dft_r2c(forwardPlan,tmp2,(fftwf_complex*)tmp1);
		memset(tmp2,0,2*(m/2+1)*n*4);
		g2s::complexAddAlphaxCxD(tmp2, tmp1, patternMem, 1.0f, (m/2+1)*n);
		fftwf_execute_dft_c2r(backwardPlan,(fftwf_complex*)tmp2,stdDev1);

		#pragma omp parallel for simd default(none) firstprivate(tmp2,im2Ptr,m,n)
		for (int i = 0; i < m*n; ++i)
		{
			tmp2[i]=sqr(im2Ptr[i]);
		}
		fftwf_execute_dft_r2c(forwardPlan,tmp2,(fftwf_complex*)tmp1);
		memset(tmp2,0,2*(m/2+1)*n*4);
		g2s::complexAddAlphaxCxD(tmp2, tmp1, patternMem, 1.0f, (m/2+1)*n);
		fftwf_execute_dft_c2r(backwardPlan,(fftwf_complex*)tmp2,stdDev2);

		#pragma omp parallel for simd default(none) firstprivate(stdDev1,stdDev2,im1Mean,im2Mean,m,n)
		for (int i = 0; i < m*n; ++i)
		{	
			stdDev1[i]=sqrt(stdDev1[i]/(m*n)-sqr(im1Mean[i]));
			stdDev2[i]=sqrt(stdDev2[i]/(m*n)-sqr(im2Mean[i]));	
		}
	}

	//printf("start resolution \n ");
	mexEvalString("drawnow");

	for (int localLagId = 0; localLagId < sizeImLags[0]; ++localLagId)
	{
		if(done){
			break;
		}
		//memset(tmp1,0,sizeof(float)*m*n);
		auto begin = std::chrono::high_resolution_clock::now();
		#pragma omp parallel for simd default(none) firstprivate(tmp1,m,n)
		for (int i = 0; i < 2*(m/2+1)*n; ++i)
		{
			tmp1[i]=0;
		}
		int deltaX=lagsPtr[0*sizeImLags[0]+localLagId];
		int deltaY=lagsPtr[1*sizeImLags[0]+localLagId];
		printf("%d,%d -- %d /%d \n", deltaX, deltaY, localLagId, sizeImLags[0]);
		mexEvalString("drawnow");

		int absolutLag=deltaY*m+deltaX;

		switch(mode){
			case SAE:
			case MAE:
			#pragma omp parallel for simd firstprivate(im1Ptr,im2Ptr,tmp1,absolutLag,sizePattern,offsetPX,offsetPY,kernel,m)
			for (int index = 0; index < m*n; ++index)
			{
				tmp1[index]=std::abs(im1Ptr[index]-im2Ptr[(index+absolutLag)%(m*n)]);
			}
			break;
			case SSE:
			case MSE:
			#pragma omp parallel for simd firstprivate(im1Ptr,im2Ptr,tmp1,absolutLag,sizePattern,offsetPX,offsetPY,kernel,m)
			for (int index = 0; index < m*n; ++index)
			{
				tmp1[index]=sqr(im1Ptr[index]-im2Ptr[(index+absolutLag)%(m*n)]);
			}
			break;
			case CC:
			#pragma omp parallel for simd firstprivate(im1Ptr,im2Ptr,tmp1,absolutLag,sizePattern,offsetPX,offsetPY,kernel,m)
			for (int index = 0; index < m*n; ++index)
			{
				tmp1[index]=(im1Ptr[index])*(im2Ptr[(index+absolutLag)%(m*n)]);
			}
			break;
			case NMAE:
			#pragma omp parallel for simd firstprivate(im1Ptr,im2Ptr,tmp1,absolutLag,sizePattern,offsetPX,offsetPY,kernel,m)
			for (int index = 0; index < m*n; ++index)
			{
				tmp1[index]=std::abs((im1Ptr[index]-im1Mean[index])-(im2Ptr[(index+absolutLag)%(m*n)]-im2Mean[(index+absolutLag)%(m*n)]));
			}
			break;
			case NMSE:
			#pragma omp parallel for simd firstprivate(im1Ptr,im2Ptr,tmp1,absolutLag,sizePattern,offsetPX,offsetPY,kernel,m)
			for (int index = 0; index < m*n; ++index)
			{
				tmp1[index]=sqr((im1Ptr[index]-im1Mean[index])-(im2Ptr[(index+absolutLag)%(m*n)]-im2Mean[(index+absolutLag)%(m*n)]));
			}
			break;
			case NCC:
			#pragma omp parallel for simd firstprivate(im1Ptr,im2Ptr,tmp1,absolutLag,sizePattern,offsetPX,offsetPY,kernel,m)
			for (int index = 0; index < m*n; ++index)
			{
				tmp1[index]=(im1Ptr[index])*(im2Ptr[(index+absolutLag)%(m*n)]);
			}
			break;
		}

		fftwf_execute(forwardPlan);
		
		
		#pragma omp parallel for simd default(none) firstprivate(tmp1,m,n)
		for (int i = 0; i < 2*(m/2+1)*n; ++i)
		{
			tmp1[i]=0;
		}

		g2s::complexAddAlphaxCxD(tmp1, tmp2, patternMem, 1.0f, (m/2+1)*n);
		fftwf_execute(backwardPlan);


		if(mode==NCC){
			for (int i = 0; i < m*n; ++i)
			{
				tmp2[i]=(tmp2[i]/(n*m)-im1Mean[i]*im2Mean[i])/(stdDev1[i]*stdDev2[i]);
			}
		}

		if(addAbsolut)
		{
			#pragma omp parallel for simd default(none) firstprivate(tmp2,m,n)
			for (int i = 0; i < m*n; ++i)
			{
				tmp2[i]=std::abs(tmp2[i]);
			}
		}

		if(lookForMax){
			#pragma omp parallel for simd default(none) firstprivate(tmp2, qualityPtr, lagsIndexPtr, localLagId, m, n)
			for (int i = 0; i < m*n; ++i)
			{
				if(tmp2[i]>qualityPtr[i]){
					qualityPtr[i]=tmp2[i];
					lagsIndexPtr[i]=localLagId+1;
				}
			}
		}else{
			#pragma omp parallel for simd default(none) firstprivate(tmp2, qualityPtr, lagsIndexPtr, localLagId, m, n)
			for (int i = 0; i < m*n; ++i)
			{
				if(tmp2[i]<qualityPtr[i]){
					qualityPtr[i]=tmp2[i];
					lagsIndexPtr[i]=localLagId+1;
				}
			}
		}

		auto end = std::chrono::high_resolution_clock::now();
		double time = 1.0e-6 * std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
		printf("compuattion time: %7.2f s\n", time/1000);
		mexEvalString("drawnow");
	}

	if(mode!=NCC){
		#pragma omp parallel for simd default(none) firstprivate(qualityPtr, m, n)
		for (int i = 0; i < m*n; ++i)
		{
			qualityPtr[i]/=(m*n);
		}
	}

	fftwf_destroy_plan(forwardPlan);
	fftwf_destroy_plan(backwardPlan);
	fftwf_cleanup_threads();
	fftwf_cleanup();

	mem_free(tmp1);
	mem_free(tmp2);
	free(patternMem);
	if(needNormalized || needDeviation){
		free(im1Mean);
		free(im2Mean);
	}
	if(needDeviation){
		free(stdDev1);
		free(stdDev2);
	}

	//finish
	if(nlhs>0) plhs[0]=lagsIndex;
	if(nlhs>1) plhs[1]=quality;
	if(done){
		done=true;
		mexErrMsgIdAndTxt("User:Interruption", "Ctrl + C interruption");
	}
	done=true;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	std::atomic<bool> done(false);
	auto myFuture = mexInterrupt::startInterruptCheck(done);
	mexFunctionWork(nlhs, plhs,  nrhs, prhs,done);
	myFuture.wait();
	mexEvalString("drawnow");
}