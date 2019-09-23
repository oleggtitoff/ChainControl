#define _CRT_SECURE_NO_WARNINGS

#include <math.h>

#include "cheetah.h"

#include "FilesOperations.h"


#define CONFIG_FILE_NAME "Config.txt"

#define SAMPLE_RATE	 48000
#define EQ_BANDS_NUM 10

#define PI	  3.14159265358979323846
#define SQRT2 1.41421356237309504880


typedef enum {
	isBypassID = 100,
	fadeTimeID = 101
} CrossFadeParamID;

typedef enum {
	noiseGateIsActiveID		= 201,
	expanderIsActiveID		= 202,
	compressorIsActiveID	= 203,
	limiterIsActiveID		= 204,

	noiseThrID				= 205,
	expanderHighThrID		= 206,
	compressorLowThrID		= 207,
	limiterThrID			= 208,

	expanderRatioID			= 209,
	compressorRatioID		= 210,

	envelopeAttackTimeID	= 211,
	envelopeReleaseTimeID	= 212,
	expanderAttackTimeID	= 213,
	expanderReleaseTimeID	= 214,
	compressorAttackTimeID	= 215,
	compressorReleaseTimeID = 216
} AmplitudeProcParamID;

typedef enum {
	isActiveID	 = 1,
	filterTypeID = 2,
	FcID		 = 3,
	Qid			 = 4,
	peakGain	 = 5
} BiquadParamID;

typedef enum {
	lowpass   = 1,
	highpass  = 2,
	bandpass  = 3,
	notch	  = 4,
	peak	  = 5,
	lowShelf  = 6,
	highShelf = 7
} BiquadType;


typedef struct {
	int sampleRate;
	int8_t isBypass;
	double fadeTime;				// in ms
} CrossFadeParams;

typedef struct {
	int sampleRate;

	int8_t noiseGateIsActive;
	int8_t expanderIsActive;
	int8_t compressorIsActive;
	int8_t limiterIsActive;

	double noiseThr;				// in dB
	double expanderHighThr;			// in dB
	double compressorLowThr;		// in dB
	double limiterThr;				// in dB

	double expanderRatio;
	double compressorRatio;

	double envelopeAttackTime;		// in ms
	double envelopeReleaseTime;		// in ms
	double expanderAttackTime;		// in ms
	double expanderReleaseTime;		// in ms
	double compressorAttackTime;	// in ms
	double compressorReleaseTime;	// in ms
} AmplitudeProcParams;

typedef struct {
	int sampleRate;
	int8_t isActive;
	int8_t type;
	double Fc;
	double Q;
	double peakGain;
} BiquadParams;

typedef struct {
	int8_t isActive;
	BiquadParams biquadParams[EQ_BANDS_NUM];
} EQParams;

typedef struct {
	double inputGain;					// in dB
	double outputGain;					// in dB
	double balance;
	CrossFadeParams crossFadeParams;
	AmplitudeProcParams amplitudeProcParams;
	EQParams eqParams;
} Params;


typedef struct {
	int32_t targetGain;				// Q31
	int32_t fadeAlpha;				// Q31
} CrossFadeCoeffs;

typedef struct {
	int32_t alphaAttack;			// Q31
	int32_t alphaRelease;			// Q31
} EnvelopeCoeffs;

typedef struct {
	int32_t isActive;
	int32_t threshold;				// Q27
} LimiterCoeffs;
typedef LimiterCoeffs NoiseGateCoeffs;

typedef struct {
	int32_t isActive;
	int32_t threshold;				// Q27
	int32_t C1;						// Q27
	int32_t C2;						// Q27
	int32_t alphaAttack;			// Q31
	int32_t alphaRelease;			// Q31
} CompressorCoeffs;
typedef CompressorCoeffs ExpanderCoeffs;

typedef struct {
	EnvelopeCoeffs envelope;
	NoiseGateCoeffs noiseGate;
	LimiterCoeffs limiter;
	CompressorCoeffs compressor;
	ExpanderCoeffs expander;
} AmplitudeProcCoeffs;

typedef struct {
	int32_t isActive;
	int32_t a[3];					// Q27
	int32_t b[2];					// Q27
} BiquadCoeffs;

typedef struct {
	int32_t isActive;
	BiquadCoeffs biquadCoeffs[EQ_BANDS_NUM];
} EQCoeffs;

typedef struct {
	int32_t inputGain;				// Q31
	int32_t outputGain;				// Q31
	int32_t balanceL;				// Q31
	int32_t balanceR;				// Q31
	CrossFadeCoeffs crossFadeCoeffs;
	AmplitudeProcCoeffs amplitudeProcCoeffs;
	EQCoeffs eqCoeffs;
} Coeffs;


void init(Params *params, Coeffs *coeffs);
void readConfig(FILE *configFilePtr, Params *params);
void calcCoeffs(Params *params, Coeffs *coeffs);

s32 SPIinitDevice(Cheetah *device);
void writeCoeffs(Coeffs *coeffs, Cheetah *device);


int main()
{
	FILE *configFilePtr = openFile(CONFIG_FILE_NAME, read);

	Params params;
	Coeffs coeffs;

	init(&params, &coeffs);
	readConfig(configFilePtr, &params);
	calcCoeffs(&params, &coeffs);

	Cheetah device;
	SPIinitDevice(&device);
	writeCoeffs(&coeffs, &device);

	ch_close(device);
	return 0;
}


//
// *** INITIALIZATION ***
//
// PARAMS INITIALIZATION
//
void initCrossFadeParams(CrossFadeParams *params)
{
	params->sampleRate = SAMPLE_RATE;
	params->isBypass   = 1;
	params->fadeTime   = 7;
}

void initAmplitudeProcParams(AmplitudeProcParams *params)
{
	params->sampleRate			  = SAMPLE_RATE;

	params->noiseGateIsActive	  = 0;
	params->expanderIsActive	  = 0;
	params->compressorIsActive	  = 0;
	params->limiterIsActive		  = 0;

	params->noiseThr			  = 0;
	params->expanderHighThr		  = 0;
	params->compressorLowThr	  = 0;
	params->limiterThr			  = 0;

	params->expanderRatio		  = 1;
	params->compressorRatio		  = 1;

	params->envelopeAttackTime    = 0;
	params->envelopeReleaseTime   = 0;
	params->expanderAttackTime    = 0;
	params->expanderReleaseTime   = 0;
	params->compressorAttackTime  = 0;
	params->compressorReleaseTime = 0;
}

void initBiquadParams(BiquadParams *params)
{
	params->sampleRate = SAMPLE_RATE;
	params->isActive   = 0;
	params->type	   = 1;
	params->Fc		   = 0;
	params->Q		   = 0;
	params->peakGain   = 0;
}

void initEQParams(EQParams *params)
{
	params->isActive = 0;

	int8_t i;

	for (i = 0; i < EQ_BANDS_NUM; i++)
	{
		initBiquadParams(&params->biquadParams[i]);
	}
}

void initParams(Params *params)
{
	params->inputGain  = 0;
	params->outputGain = 0;
	params->balance	   = 0;

	initCrossFadeParams(&params->crossFadeParams);
	initAmplitudeProcParams(&params->amplitudeProcParams);
	initEQParams(&params->eqParams);
}

//
// COEFFS INIITIALIZATION
//
void initCrossFadeCoeffs(CrossFadeCoeffs *coeffs)
{
	coeffs->targetGain = 0x7fffffff;
	coeffs->fadeAlpha  = 0;
}

void initEnvelopeCoeffs(EnvelopeCoeffs *coeffs)
{
	coeffs->alphaAttack	 = 0;
	coeffs->alphaRelease = 0;
}

void initNoiseGateCoeffs(NoiseGateCoeffs *coeffs)
{
	coeffs->isActive  = 0;
	coeffs->threshold = 0;
}

void initLimiterCoeffs(LimiterCoeffs *coeffs)
{
	coeffs->isActive  = 0;
	coeffs->threshold = 0;
}

void initCompressorCoeffs(CompressorCoeffs *coeffs)
{
	coeffs->isActive	 = 0;
	coeffs->threshold	 = 0;

	coeffs->C1			 = 0;
	coeffs->C2			 = 0;
	
	coeffs->alphaAttack	 = 0;
	coeffs->alphaRelease = 0;
}

void initExpanderCoeffs(ExpanderCoeffs *coeffs)
{
	coeffs->isActive	 = 0;
	coeffs->threshold	 = 0;

	coeffs->C1			 = 0;
	coeffs->C2			 = 0;

	coeffs->alphaAttack  = 0;
	coeffs->alphaRelease = 0;
}

void initAmplitudeProcCoeffs(AmplitudeProcCoeffs *coeffs)
{
	initEnvelopeCoeffs(&coeffs->envelope);
	initNoiseGateCoeffs(&coeffs->noiseGate);
	initLimiterCoeffs(&coeffs->limiter);
	initCompressorCoeffs(&coeffs->compressor);
	initExpanderCoeffs(&coeffs->expander);
}

void initBiquadCoeffs(BiquadCoeffs *coeffs)
{
	coeffs->isActive = 0;

	coeffs->a[0]	 = 0;
	coeffs->a[1]	 = 0;
	coeffs->a[2]	 = 0;

	coeffs->b[0]	 = 0;
	coeffs->b[1]	 = 0;
}

void initEQCoeffs(EQCoeffs *coeffs)
{
	coeffs->isActive = 0;

	int8_t i;

	for (i = 0; i < EQ_BANDS_NUM; i++)
	{
		initBiquadCoeffs(&coeffs->biquadCoeffs[i]);
	}
}

void initCoeffs(Coeffs *coeffs)
{
	coeffs->inputGain  = 0x7fffffff;
	coeffs->outputGain = 0x7fffffff;
	coeffs->balanceL   = 0x7fffffff;
	coeffs->balanceR   = 0x7fffffff;

	initCrossFadeCoeffs(&coeffs->crossFadeCoeffs);
	initAmplitudeProcCoeffs(&coeffs->amplitudeProcCoeffs);
	initEQCoeffs(&coeffs->eqCoeffs);
}

//
// GENERAL INITIALIZATION
//
void init(Params *params, Coeffs *coeffs)
{
	initParams(params);
	initCoeffs(coeffs);
}
//
// *** INITIALIZATION END ***
//


double dBtoGain(const double dB)
{
	return pow(10, dB / 20.0);
}

int32_t doubleToFixed31(double x)
{
	if (x >= 1)
	{
		return INT32_MAX;
	}
	else if (x < -1)
	{
		return INT32_MIN;
	}

	return (int32_t)(x * (double)(1LL << 31));
}

double alphaCalc(const int sampleRate, const double time)
{
	// calculates alpha coefficient from sample rate and time in ms

	return (double)1 - exp((double)-1 / (sampleRate * time / 1000));
}


//
// *** COEFFS CALCULATION ***
//
void calcCrossFadeCoeffs(const CrossFadeParams *params, CrossFadeCoeffs *coeffs)
{
	coeffs->targetGain = doubleToFixed31(params->isBypass);
	coeffs->fadeAlpha  = doubleToFixed31(alphaCalc(params->sampleRate, params->fadeTime));
}

double calcC1(const double ratio)
{
	return (double)1 - ((double)1 / ratio);
}

double calcC2(const double threshold, const double C1)
{
	return pow(dBtoGain(threshold), C1);
}

void calcAmplitudeProcCoeffs(AmplitudeProcParams *params, AmplitudeProcCoeffs *coeffs)
{
	// envelope coeffs
	coeffs->envelope.alphaAttack  = doubleToFixed31(alphaCalc(params->sampleRate, 
															  params->envelopeAttackTime));
	coeffs->envelope.alphaRelease = doubleToFixed31(alphaCalc(params->sampleRate,
															  params->envelopeReleaseTime));

	// noise gate coeffs
	coeffs->noiseGate.isActive  = params->noiseGateIsActive;
	coeffs->noiseGate.threshold = doubleToFixed31(dBtoGain(params->noiseThr) / 16);

	// limiter coeffs
	coeffs->limiter.isActive  = params->limiterIsActive;
	coeffs->limiter.threshold = doubleToFixed31(dBtoGain(params->limiterThr) / 16);

	// compressor coeffs
	coeffs->compressor.isActive		= params->compressorIsActive;
	coeffs->compressor.threshold	= doubleToFixed31(dBtoGain(params->compressorLowThr) / 16);

	double C1						= calcC1(params->compressorRatio);
	coeffs->compressor.C1			= doubleToFixed31(C1 / 16);
	coeffs->compressor.C2			= doubleToFixed31(calcC2(params->compressorLowThr, C1) / 16);

	coeffs->compressor.alphaAttack	= doubleToFixed31(alphaCalc(params->sampleRate, 
													  params->compressorAttackTime));
	coeffs->compressor.alphaRelease = doubleToFixed31(alphaCalc(params->sampleRate, 
													  params->compressorReleaseTime));

	// expander coeffs
	coeffs->expander.isActive	  = params->expanderIsActive;
	coeffs->expander.threshold	  = doubleToFixed31(dBtoGain(params->expanderHighThr) / 16);

	C1							  = calcC1(params->expanderRatio);
	coeffs->expander.C1			  = doubleToFixed31(C1 / 16);
	coeffs->expander.C2			  = doubleToFixed31(calcC2(params->expanderHighThr, C1) / 16);

	coeffs->expander.alphaAttack  = doubleToFixed31(alphaCalc(params->sampleRate, 
															  params->expanderAttackTime));
	coeffs->expander.alphaRelease = doubleToFixed31(alphaCalc(params->sampleRate, 
															  params->expanderReleaseTime));
}

void calcBiquadCoeffs(BiquadParams *params, BiquadCoeffs *coeffs)
{
	coeffs->isActive = params->isActive;
	
	double norm;
	double Q = params->Q;
	double V = pow(10, fabs(params->peakGain) / 20);
	double K = tan(PI * params->Fc / params->sampleRate);
	double Kpow2 = K * K;
	double a[3];
	double b[2];

	switch (params->type)
	{
	case lowpass:
		norm = 1 / (1 + K / Q + Kpow2);
		a[0] = Kpow2 * norm;
		a[1] = 2 * a[0];
		a[2] = a[0];
		b[0] = 2 * (Kpow2 - 1) * norm;
		b[1] = (1 - K / Q + Kpow2) * norm;
		break;

	case highpass:
		norm = 1 / (1 + K / Q + Kpow2);
		a[0] = norm;
		a[1] = -2 * a[0];
		a[2] = a[0];
		b[0] = 2 * (Kpow2 - 1) * norm;
		b[1] = (1 - K / Q + Kpow2) * norm;
		break;

	case bandpass:
		norm = 1 / (1 + K / Q + Kpow2);
		a[0] = K / Q * norm;
		a[1] = 0;
		a[2] = -a[0];
		b[0] = 2 * (Kpow2 - 1) * norm;
		b[1] = (1 - K / Q + Kpow2) * norm;
		break;

	case notch:
		norm = 1 / (1 + K / Q + Kpow2);
		a[0] = (1 + Kpow2) * norm;
		a[1] = 2 * (Kpow2 - 1) * norm;
		a[2] = a[0];
		b[0] = a[1];
		b[1] = (1 - K / Q + Kpow2) * norm;
		break;

	case peak:
		if (params->peakGain >= 0)
		{
			norm = 1 / (1 + 1 / Q * K + Kpow2);
			a[0] = (1 + V / Q * K + Kpow2) * norm;
			a[1] = 2 * (Kpow2 - 1) * norm;
			a[2] = (1 - V / Q * K + Kpow2) * norm;
			b[0] = a[1];
			b[1] = (1 - 1 / Q * K + Kpow2) * norm;
		}
		else
		{
			norm = 1 / (1 + V / Q * K + Kpow2);
			a[0] = (1 + 1 / Q * K + Kpow2) * norm;
			a[1] = 2 * (Kpow2 - 1) * norm;
			a[2] = (1 - 1 / Q * K + Kpow2) * norm;
			b[0] = a[1];
			b[1] = (1 - V / Q * K + Kpow2) * norm;
		}
		break;

	case lowShelf:
		if (params->peakGain >= 0)
		{
			norm = 1 / (1 + SQRT2 * K + Kpow2);
			a[0] = (1 + sqrt(2 * V) * K + V * Kpow2) * norm;
			a[1] = 2 * (V * Kpow2 - 1) * norm;
			a[2] = (1 - sqrt(2 * V) * K + V * Kpow2) * norm;
			b[0] = 2 * (Kpow2 - 1) * norm;
			b[1] = (1 - SQRT2 * K + Kpow2) * norm;
		}
		else
		{
			norm = 1 / (1 + sqrt(2 * V) * K + V * Kpow2);
			a[0] = (1 + SQRT2 * K + Kpow2) * norm;
			a[1] = 2 * (Kpow2 - 1) * norm;
			a[2] = (1 - SQRT2 * K + Kpow2) * norm;
			b[0] = 2 * (V * Kpow2 - 1) * norm;
			b[1] = (1 - sqrt(2 * V) * K + V * Kpow2) * norm;
		}
		break;

	case highShelf:
		if (params->peakGain >= 0)
		{
			norm = 1 / (1 + SQRT2 * K + Kpow2);
			a[0] = (V + sqrt(2 * V) * K + Kpow2) * norm;
			a[1] = 2 * (Kpow2 - V) * norm;
			a[2] = (V - sqrt(2 * V) * K + Kpow2) * norm;
			b[0] = 2 * (Kpow2 - 1) * norm;
			b[1] = (1 - SQRT2 * K + Kpow2) * norm;
		}
		else
		{
			norm = 1 / (V + sqrt(2 * V) * K + Kpow2);
			a[0] = (1 + SQRT2 * K + Kpow2) * norm;
			a[1] = 2 * (Kpow2 - 1) * norm;
			a[2] = (1 - SQRT2 * K + Kpow2) * norm;
			b[0] = 2 * (Kpow2 - V) * norm;
			b[1] = (V - sqrt(2 * V) * K + Kpow2) * norm;
		}
		break;
	}

	coeffs->a[0] = doubleToFixed31(a[0] / 16);
	coeffs->a[1] = doubleToFixed31(a[1] / 16);
	coeffs->a[2] = doubleToFixed31(a[2] / 16);
	coeffs->b[0] = doubleToFixed31(b[0] / 16);
	coeffs->b[0] = doubleToFixed31(b[1] / 16);
}

void calcEQCoeffs(EQParams *params, EQCoeffs *coeffs)
{
	coeffs->isActive = params->isActive;

	int8_t i;

	for (i = 0; i < EQ_BANDS_NUM; i++)
	{
		calcBiquadCoeffs(&params->biquadParams[i], &coeffs->biquadCoeffs[i]);
	}
}

void calcCoeffs(Params *params, Coeffs *coeffs)
{
	coeffs->inputGain  = doubleToFixed31(dBtoGain(params->inputGain));
	coeffs->outputGain = doubleToFixed31(dBtoGain(params->outputGain));
	
	double balance = params->balance / 10;

	if (balance < 0)
	{
		coeffs->balanceL = 0x7fffffff;
		coeffs->balanceR = doubleToFixed31(1 + balance);
	}
	else if (balance > 0)
	{
		coeffs->balanceL = doubleToFixed31(1 - balance);
		coeffs->balanceR = 0x7fffffff;
	}

	calcCrossFadeCoeffs(&params->crossFadeParams, &coeffs->crossFadeCoeffs);
	calcAmplitudeProcCoeffs(&params->amplitudeProcParams, &coeffs->amplitudeProcCoeffs);
	calcEQCoeffs(&params->eqParams, &coeffs->eqCoeffs);
}
//
// *** COEFFS CALCULATION END ***
//


void parseConfigString(char *str, uint16_t *paramId, double *paramValue)
{
	uint8_t inputIndex = 0;
	uint8_t idIndex = 0;
	uint8_t paramIndex = 0;
	uint8_t isIndex = 1;
	char idStr[6] = { 0 };
	char paramStr[13] = { 0 };

	while (str[inputIndex])
	{
		if (str[inputIndex] == '/' ||
			((str[inputIndex] == ' ' || str[inputIndex] == '	') && paramIndex > 0))
		{
			break;
		}
		else if (str[inputIndex] == ':')
		{
			isIndex = 0;
		}
		else if (str[inputIndex] != ' ' && str[inputIndex] != '	' && str[inputIndex] != '\n')
		{
			if (isIndex)
			{
				idStr[idIndex] = str[inputIndex];
				idIndex++;
			}
			else
			{
				paramStr[paramIndex] = str[inputIndex];
				paramIndex++;
			}
		}

		inputIndex++;
	}

	*paramId = atoi(idStr);
	*paramValue = atof(paramStr);
}


//
// *** SETTING PARAMS ***
//
void crossFadeSetParam(CrossFadeParams *params, const uint16_t id, const double value)
{
	switch (id)
	{
	case isBypassID:
		params->isBypass = (int8_t)value;
		break;

	case fadeTimeID:
		params->fadeTime = value;
		break;
	}
}

void amplitudeProcSetParam(AmplitudeProcParams *params, const uint16_t id, double value)
{
	switch (id)
	{
	case noiseGateIsActiveID:
		params->noiseGateIsActive = (int8_t)value;
		break;

	case expanderIsActiveID:
		params->expanderIsActive = (int8_t)value;
		break;

	case compressorIsActiveID:
		params->compressorIsActive = (int8_t)value;
		break;

	case limiterIsActiveID:
		params->limiterIsActive = (int8_t)value;
		break;

	case noiseThrID:
		if (value > -70)
		{
			value = -70;
		}

		params->noiseThr = value;
		break;

	case expanderHighThrID:
		params->expanderHighThr = value;
		break;

	case compressorLowThrID:
		params->compressorLowThr = value;
		break;

	case limiterThrID:
		params->limiterThr = value;
		break;

	case expanderRatioID:
		params->expanderRatio = value;
		break;

	case compressorRatioID:
		params->compressorRatio = value;
		break;

	case envelopeAttackTimeID:
		params->envelopeAttackTime = value;
		break;

	case envelopeReleaseTimeID:
		params->envelopeReleaseTime = value;
		break;

	case expanderAttackTimeID:
		params->expanderAttackTime = value;
		break;

	case expanderReleaseTimeID:
		params->expanderReleaseTime = value;
		break;

	case compressorAttackTimeID:
		params->compressorAttackTime = value;
		break;

	case compressorReleaseTimeID:
		params->compressorReleaseTime = value;
		break;
	}
}

void biquadSetParam(BiquadParams *params, const uint16_t id, const double value)
{
	switch (id)
	{
	case isActiveID:
		params->isActive = (int8_t)value;
		break;

	case filterTypeID:
		params->type = (int8_t)value;
		break;

	case FcID:
		params->Fc = value;
		break;

	case Qid:
		params->Q = value;
		break;

	case peakGain:
		params->peakGain = value;
		break;
	}
}

void eqSetParam(EQParams *params, const int8_t band, const uint16_t id, const double value)
{
	if (band < 0)
	{
		params->isActive = (int)value;
	}
	else
	{
		biquadSetParam(&params->biquadParams[band], id, value);
	}
}
//
// *** SETTING PARAMS END ***
//


void readConfig(FILE *configFilePtr, Params *params)
{
	char str[100] = { 0 };
	uint16_t paramId;
	double paramValue;

	while (fgets(str, 100, configFilePtr))
	{
		parseConfigString(str, &paramId, &paramValue);

		if (paramId < 150)
		{
			crossFadeSetParam(&params->crossFadeParams, paramId, paramValue);
		}
		else if (paramId == 151)
		{
			params->inputGain = paramValue;
		}
		else if (paramId == 152)
		{
			params->outputGain = paramValue;
		}
		else if (paramId == 153)
		{
			params->balance = paramValue;
		}
		else if (paramId < 500)
		{
			amplitudeProcSetParam(&params->amplitudeProcParams, paramId, paramValue);
		}
		else
		{
			eqSetParam(&params->eqParams, paramId / 1000 - 1, paramId % 1000, paramValue);
		}
	}
}


//
// *** SPI CONTROL ***
//
s32 SPIinitDevice(Cheetah *device)
{
	s32 status;
	s32 bitrate;
	u16 devsPortNums[1] = { 0 };
	int devsNum = ch_find_devices(1, devsPortNums);
	*device = ch_open(devsPortNums[0]);

	status = ch_spi_configure(*device,
					 CH_SPI_POL_RISING_FALLING,
					 CH_SPI_PHASE_SAMPLE_SETUP,
					 CH_SPI_BITORDER_MSB,
					 0x0);
	bitrate = ch_spi_bitrate(*device, 20000);

	if (bitrate < 0)
		return bitrate;

	return status;
}

s32 SPIshift(Cheetah *device, u08 *dataOut, u32 nBytesOut,
	u08 *dataIn, u32 nBytesIn)
{
	s32 status;

	status = ch_spi_queue_clear(*device);
	if (status < 0)
		return status;

	status = ch_spi_queue_oe(*device, 1);
	if (status < 0)
		return status;

	status = ch_spi_queue_ss(*device, 0);
	if (status < 0)
		return status;

	status = ch_spi_queue_ss(*device, 1);
	if (status < 0)
		return status;

	status = ch_spi_queue_array(*device, nBytesOut, dataOut);
	if (status < 0)
		return status;

	status = ch_spi_queue_ss(*device, 0);
	if (status < 0)
		return status;

	status = ch_spi_queue_oe(*device, 0);
	if (status < 0)
		return status;

	status = ch_spi_batch_shift(*device, nBytesIn, dataIn);

	return status;
}

u32 reverse32(u32 x)
{
	return (x >> 24) |
		   (x << 24) |
		   ((x & 0x00ff0000) >> 8) |
		   ((x & 0x0000ff00) << 8);
}

s32 SPIwrite32(Cheetah *device, u32 address, u32 *data)
{
	u08 *dataArr;
	u08 CMD = 0x03;
	u08 dummy = 0;
	u32 dataSize = sizeof(CMD) + sizeof(address) + sizeof(dummy) + sizeof(u32);
	u32 offset = 0;
	s32 status;
	u32 ackBytes = 0;

	if (!(dataArr = (u08*)malloc(dataSize)))
	{
		return -1;
	}

	dataArr[offset++] = CMD;

	*(u32*)&dataArr[offset] = reverse32(address);
	offset += sizeof(address);

	*(u32*)&dataArr[offset] = reverse32(*data);
	offset += sizeof(*data);

	dataArr[offset++] = dummy;

	status = SPIshift(device, dataArr, dataSize, dataArr, dataSize);
	if (status < 0)
		return status;

	if (0x80 == dataArr[dataSize - 1])
		ackBytes = 4;

	free(dataArr);
	return ackBytes;
}

s32 SPIread32(Cheetah *device, u32 address, u32 *data)
{
	u08 *dataArr;
	u08 CMD = 0x02;
	u32 dummy = 0;
	u32 dataSize = sizeof(CMD) + sizeof(address) + sizeof(dummy) + sizeof(u32);
	u32 offset = 0;
	s32 status;
	u32 *pData;

	if (!(dataArr = (u08*)malloc(dataSize)))
	{
		return -1;
	}

	dataArr[offset++] = CMD;
	
	*(u32*)&dataArr[offset] = reverse32(address);
	offset += sizeof(address);

	*(u32*)&dataArr[offset] = dummy;
	offset += sizeof(dummy);

	status = SPIshift(device, dataArr, dataSize, dataArr, dataSize);
	if (status < 0)
		return status;

	pData = (u32*)&(dataArr[9]);
	*data = reverse32(pData[0]);
	
	free(dataArr);
	return 4;
}

void writeCoeffs(Coeffs *coeffs, Cheetah *device)
{
	u08 i;
	u32 address = 0x60003904;

	SPIwrite32(device, address, &coeffs->inputGain);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->outputGain);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->balanceL);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->balanceR);
	address += 0x4;
	
	// CrossFade
	SPIwrite32(device, address, &coeffs->crossFadeCoeffs.targetGain);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->crossFadeCoeffs.fadeAlpha);
	address += 0x4;
	
	// Envelope
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.envelope.alphaAttack);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.envelope.alphaRelease);
	address += 0x4;
	
	// NoiseGate	
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.noiseGate.isActive);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.noiseGate.threshold);
	address += 0x4;
	
	// Limiter	
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.limiter.isActive);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.limiter.threshold);
	address += 0x4;

	// Compressor
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.compressor.isActive);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.compressor.threshold);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.compressor.C1);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.compressor.C2);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.compressor.alphaAttack);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.compressor.alphaRelease);
	address += 0x4;

	// Expander
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.expander.isActive);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.expander.threshold);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.expander.C1);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.expander.C2);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.expander.alphaAttack);
	address += 0x4;
	SPIwrite32(device, address, &coeffs->amplitudeProcCoeffs.expander.alphaRelease);
	address += 0x4;

	// EQ
	SPIwrite32(device, address, &coeffs->eqCoeffs.isActive);
	address += 0x4;

	for (i = 0; i < EQ_BANDS_NUM; i++)
	{
		SPIwrite32(device, address, &coeffs->eqCoeffs.biquadCoeffs[i].isActive);
		address += 0x4;

		SPIwrite32(device, address, &coeffs->eqCoeffs.biquadCoeffs[i].a[0]);
		address += 0x4;

		SPIwrite32(device, address, &coeffs->eqCoeffs.biquadCoeffs[i].a[1]);
		address += 0x4;

		SPIwrite32(device, address, &coeffs->eqCoeffs.biquadCoeffs[i].a[2]);
		address += 0x4;

		SPIwrite32(device, address, &coeffs->eqCoeffs.biquadCoeffs[i].b[0]);
		address += 0x4;

		SPIwrite32(device, address, &coeffs->eqCoeffs.biquadCoeffs[i].b[1]);
		address += 0x4;
	}

	// isReadyToUpdate
	u32 isReadyToUpdate = 0x1;
	SPIwrite32(device, 0x60003900, &isReadyToUpdate);
}
//
// *** SPI CONTROL END ***
//