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


extern VALUE rb_double_array_plus(VALUE self, VALUE other);
extern VALUE rb_double_array_mul(VALUE self, VALUE other);

extern void Init_vdsp();
