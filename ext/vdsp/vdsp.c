#include "vdsp.h"

VALUE rb_mVdsp;

void
Init_vdsp(void)
{
  rb_mVdsp = rb_define_module("Vdsp");
}
