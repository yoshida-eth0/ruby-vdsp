#include <ruby.h>
#include <Accelerate/Accelerate.h>


// VdspArray

typedef union {
  void *ptr;
  double *d;
} VdspArrayValue;

typedef struct {
  VdspArrayValue v;
  unsigned long length;
} VdspArrayNativeResource;

typedef struct {
  VdspArrayNativeResource *res0;
  VdspArrayNativeResource *res1;
  long offset;
  vDSP_Stride stride;
} VdspArrayParam;


// VdspBiquad

typedef union {
  void *ptr;
  struct vDSP_biquad_SetupStructD *d;
} VdspBiquadSetup;

typedef struct {
  char type;
  union {
    void *ptr;
    double *d;
  } coefs;
  union {
    void *ptr;
    double *d;
  } delay;
  VdspBiquadSetup setup;
  unsigned long sections;
  unsigned long alloc_sections;
} VdspBiquadNativeResource;


// VdspFFT

typedef union {
  void *value;
  FFTSetupD d;
} VdspFFTSetup;

typedef struct {
  char type;
  VdspFFTSetup setup;
  unsigned long length;
  unsigned long halflength;
  vDSP_Length log2n;
  FFTRadix radix;
} VdspFFTNativeResource;


// VdspSplitComplex

typedef struct {
  VdspArrayNativeResource *real;
  VdspArrayNativeResource *imag;
} VdspSplitComplexNativeResource;


extern VALUE rb_double_array_plus(VALUE self, VALUE other);
extern VALUE rb_double_array_mul(VALUE self, VALUE other);
extern void double_array_resize(VdspArrayNativeResource *_a, unsigned long len);

extern void Init_vdsp();
