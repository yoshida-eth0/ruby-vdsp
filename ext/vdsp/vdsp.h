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

extern VALUE rb_double_array_plus(VALUE self, VALUE other);
extern VALUE rb_double_array_mul(VALUE self, VALUE other);

extern void Init_vdsp();
