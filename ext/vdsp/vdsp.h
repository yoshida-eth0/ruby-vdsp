#include <ruby.h>
#include <Accelerate/Accelerate.h>

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

extern VALUE rb_double_array_plus(VALUE self, VALUE other);
extern VALUE rb_double_array_mul(VALUE self, VALUE other);

extern void Init_vdsp();
