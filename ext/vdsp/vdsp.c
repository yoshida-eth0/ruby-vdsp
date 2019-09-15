#include "vdsp.h"


VALUE rb_mVdsp;

VALUE rb_mVdspScalar;
VALUE rb_mVdspArray;

VALUE rb_cDoubleScalar;
VALUE rb_cDoubleArray;
VALUE rb_mDouble;


// Native resource

void vdsp_array_native_resource_delete(VdspArrayNativeResource * p)
{
  if (p->v.ptr) {
    free(p->v.ptr);
  }
  free(p);
}

VdspArrayNativeResource* get_vdsp_array_native_resource(VALUE va)
{
  if (!rb_obj_is_kind_of(va, rb_mVdspArray)) {
    rb_raise(rb_eArgError, "Vdsp::Array required");
  }

  VALUE resource = rb_iv_get(va, "native_resource");
  if (resource==Qnil) {
    return NULL;
  }

  VdspArrayNativeResource *p;
  Data_Get_Struct(resource, VdspArrayNativeResource, p);

  return p;
}

VdspArrayValue get_vdsp_array_value(VALUE va)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(va);
  return p->v;
}

double* get_double_array_value(VALUE va)
{
  return get_vdsp_array_value(va).d;
}

long get_vdsp_array_length(VALUE va)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(va);
  return _a->length;
}

void array_param(VdspArrayParam *param, VALUE arr, VALUE offset, VALUE stride)
{
  param->res = get_vdsp_array_native_resource(arr);
  param->offset = NUM2LONG(offset);
  param->stride = NUM2LONG(stride);
}


// Vdsp::DoubleScalar

VALUE rb_double_scalar_initialize(VALUE self, VALUE val)
{
  val = rb_funcall(val, rb_intern("to_f"), 0);
  rb_iv_set(self, "val", val);
  return self;
}

VALUE rb_double_scalar_plus(VALUE self, VALUE other)
{
  assert(rb_obj_is_kind_of(other, rb_mVdspArray));

  VALUE val = rb_iv_get(self, "val");
  return rb_double_array_plus(other, val);
}

VALUE rb_double_scalar_minus(VALUE self, VALUE other)
{
  assert(rb_obj_is_kind_of(other, rb_mVdspArray));

  other = rb_double_array_mul(other, DBL2NUM(-1.0));
  VALUE val = rb_iv_get(self, "val");
  return rb_double_array_plus(other, val);
}

VALUE rb_double_scalar_mul(VALUE self, VALUE other)
{
  assert(rb_obj_is_kind_of(other, rb_mVdspArray));

  VALUE val = rb_iv_get(self, "val");
  return rb_double_array_mul(other, val);
}

VALUE rb_double_scalar_div(VALUE self, VALUE other)
{
  assert(rb_obj_is_kind_of(other, rb_mVdspArray));

  VALUE val = rb_iv_get(self, "val");
  double _a = NUM2DBL(val);

  VdspArrayNativeResource *_b = get_vdsp_array_native_resource(other);
  VALUE lenv = LONG2NUM(_b->length);

  VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

  vDSP_svdivD(&_a, _b->v.d, 1, _c->v.d, 1, _b->length);

  return c;
}


// Vdsp::Array

VALUE rb_vdsp_array_length(VALUE self)
{
  return LONG2NUM(get_vdsp_array_length(self));
}


// Vdsp::DoubleArray

VALUE rb_double_array_initialize(VALUE self, VALUE length)
{
  if (!FIXNUM_P(length)) {
    rb_raise(rb_eArgError, "Integer required: length");
  }
  long _length = FIX2LONG(length);

  VdspArrayNativeResource *p = ALLOC(VdspArrayNativeResource);
  p->v.ptr = NULL;
  p->length = 0;

  VALUE resource = Data_Wrap_Struct(CLASS_OF(self), 0, vdsp_array_native_resource_delete, p);
  rb_iv_set(self, "native_resource", resource);

  p->v.ptr = calloc(_length, sizeof(double));
  p->length = _length;

  return self;
}

VALUE rb_double_array_set_values(VALUE self, VALUE ary)
{
  if (!RB_TYPE_P(ary, T_ARRAY)) {
    rb_raise(rb_eArgError, "Array required");
  }

  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);

  const VALUE *ary_p = RARRAY_CONST_PTR(ary);
  unsigned long ary_len = RARRAY_LEN(ary);

  if (p->length!=ary_len) {
    rb_raise(rb_eArgError, "Array length error: self.length=%ld src.length=%ld", p->length, ary_len);
  }

  double *d = p->v.d;

  for (unsigned long i=0; i<p->length; i++) {
    d[i] = NUM2DBL(ary_p[i]);
  }

  return self;
}

VALUE rb_double_array_get_values(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double *d = p->v.d;

  VALUE ret = rb_ary_new2(p->length);
  for (unsigned long i=0; i<p->length; i++) {
    rb_ary_push(ret, DBL2NUM(d[i]));
  }

  return ret;
}

VALUE rb_double_array_create(VALUE cls, VALUE ary)
{
  if (rb_obj_is_kind_of(ary, rb_mVdspArray)) {
    return rb_funcall(ary, rb_intern("to_da"), 0);
  }
  if (!RB_TYPE_P(ary, T_ARRAY)) {
    rb_raise(rb_eArgError, "Array required");
  }

  VALUE len = LONG2NUM(RARRAY_LEN(ary));
  VALUE obj = rb_class_new_instance(1, &len, rb_cDoubleArray);
  rb_double_array_set_values(obj, ary);

  return obj;
}

VALUE rb_double_array_to_da(VALUE self)
{
  return self;
}

VALUE rb_double_array_plus(VALUE self, VALUE other)
{
  if (rb_obj_is_kind_of(other, rb_mVdspArray)) {
    other = rb_funcall(other, rb_intern("to_da"), 0);

    VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
    VdspArrayNativeResource *_b = get_vdsp_array_native_resource(other);

    vDSP_Length len = MIN(_a->length, _b->length);
    VALUE lenv = LONG2NUM(len);
    VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
    VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

    vDSP_vaddD(_a->v.d, 1, _b->v.d, 1, _c->v.d, 1, len);
    return c;

  } else if (rb_obj_is_kind_of(other, rb_cNumeric)) {
    other = rb_funcall(other, rb_intern("to_f"), 0);
    double _b = NUM2DBL(other);

    VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
    VALUE lenv = LONG2NUM(_a->length);
    VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
    VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

    vDSP_vsaddD(_a->v.d, 1, &_b, _c->v.d, 1, _a->length);
    return c;

  } else {
    return rb_num_coerce_bin(self, other, '+');
  }
}

VALUE rb_double_array_minus(VALUE self, VALUE other)
{
  if (rb_obj_is_kind_of(other, rb_mVdspArray)) {
    other = rb_funcall(other, rb_intern("to_da"), 0);

    VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
    VdspArrayNativeResource *_b = get_vdsp_array_native_resource(other);

    vDSP_Length len = MIN(_a->length, _b->length);
    VALUE lenv = LONG2NUM(len);
    VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
    VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

    vDSP_vsubD(_b->v.d, 1, _a->v.d, 1, _c->v.d, 1, len);
    return c;

  } else if (rb_obj_is_kind_of(other, rb_cNumeric)) {
    other = rb_funcall(other, rb_intern("to_f"), 0);
    double _b = -NUM2DBL(other);

    VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
    VALUE lenv = LONG2NUM(_a->length);
    VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
    VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

    vDSP_vsaddD(_a->v.d, 1, &_b, _c->v.d, 1, _a->length);
    return c;

  } else {
    return rb_num_coerce_bin(self, other, '-');
  }
}

VALUE rb_double_array_mul(VALUE self, VALUE other)
{
  if (rb_obj_is_kind_of(other, rb_mVdspArray)) {
    other = rb_funcall(other, rb_intern("to_da"), 0);

    VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
    VdspArrayNativeResource *_b = get_vdsp_array_native_resource(other);

    vDSP_Length len = MIN(_a->length, _b->length);
    VALUE lenv = LONG2NUM(len);
    VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
    VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

    vDSP_vmulD(_a->v.d, 1, _b->v.d, 1, _c->v.d, 1, len);
    return c;

  } else if (rb_obj_is_kind_of(other, rb_cNumeric)) {
    other = rb_funcall(other, rb_intern("to_f"), 0);
    double _b = NUM2DBL(other);

    VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
    VALUE lenv = LONG2NUM(_a->length);
    VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
    VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

    vDSP_vsmulD(_a->v.d, 1, &_b, _c->v.d, 1, _a->length);
    return c;

  } else {
    return rb_num_coerce_bin(self, other, '*');
  }
}

VALUE rb_double_array_div(VALUE self, VALUE other)
{
  if (rb_obj_is_kind_of(other, rb_mVdspArray)) {
    other = rb_funcall(other, rb_intern("to_da"), 0);

    VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
    VdspArrayNativeResource *_b = get_vdsp_array_native_resource(other);

    vDSP_Length len = MIN(_a->length, _b->length);
    VALUE lenv = LONG2NUM(len);
    VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
    VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

    vDSP_vdivD(_a->v.d, 1, _b->v.d, 1, _c->v.d, 1, len);
    return c;

  } else if (rb_obj_is_kind_of(other, rb_cNumeric)) {
    other = rb_funcall(other, rb_intern("to_f"), 0);
    double _b = NUM2DBL(other);

    VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
    VALUE lenv = LONG2NUM(_a->length);
    VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
    VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

    vDSP_vsdivD(_a->v.d, 1, &_b, _c->v.d, 1, _a->length);
    return c;

  } else {
    return rb_num_coerce_bin(self, other, '/');
  }
}

VALUE rb_double_array_aref(VALUE self, VALUE i)
{
  long _i = NUM2LONG(i);

  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  if (0<=_i && (unsigned long)_i<p->length) {
    return DBL2NUM(p->v.d[_i]);
  } else {
    rb_raise(rb_eIndexError, "Index out of range: %ld", _i);
  }
}

VALUE rb_double_array_aset(VALUE self, VALUE i, VALUE val)
{
  long _i = NUM2LONG(i);

  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  if (0<=_i && (unsigned long)_i<p->length) {
    val = rb_funcall(val, rb_intern("to_f"), 0);
    p->v.d[_i] = NUM2DBL(val);
    return val;
  } else {
    rb_raise(rb_eIndexError, "Index out of range: %ld", _i);
  }
}

VALUE rb_double_array_coerce(VALUE self, VALUE other)
{
  other = rb_class_new_instance(1, &other, rb_cDoubleScalar);
  return rb_assoc_new(other, self);
}


// Vdsp static method

// c[i] = a[i] + b
VALUE rb_double_vsadd(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  double _b = NUM2DBL(b);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vsaddD(
    _a.res->v.d+_a.offset, _a.stride,
    &_b,
    _c.res->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

// c[i] = a[i] + b[i]
VALUE rb_double_vadd(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vaddD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _b.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

// c[i] = a[i] - b[i]
VALUE rb_double_vsub(
  VALUE cls,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _b;
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_b, b, b_offset, b_stride);
  array_param(&_a, a, a_offset, a_stride);
  array_param(&_c, c, c_offset, c_stride);
  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vsubD(
    _b.res->v.d+_b.offset, _b.stride,
    _a.res->v.d+_a.offset, _a.stride,
    _c.res->v.d+_c.offset, _c.stride,
     _n);

  return c;
}

// c[i] = a[i] * b
VALUE rb_double_vsmul(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  double _b = NUM2DBL(b);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vsmulD(
    _a.res->v.d+_a.offset, _a.stride,
    &_b,
    _c.res->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

// c[i] = a[i] * b[i]
VALUE rb_double_vmul(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vmulD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _b.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

// c[i] = a[i] / b
VALUE rb_double_vsdiv(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  double _b = NUM2DBL(b);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vsdivD(
    _a.res->v.d+_a.offset, _a.stride,
    &_b,
    _c.res->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

// c[i] = a / b[i]
VALUE rb_double_svdiv(
  VALUE cls,
  VALUE a,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _b;
  VdspArrayParam _c;

  double _a = NUM2DBL(a);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_svdivD(
    &_a,
    _b.res->v.d+_b.offset, _b.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

// o0[i] = i1[i] + i0[i]
// o1[i] = i1[i] - i0[i]
VALUE rb_double_vaddsub(
  VALUE cls,
  VALUE i0, VALUE i0_offset, VALUE i0_stride,
  VALUE i1, VALUE i1_offset, VALUE i1_stride,
  VALUE o0, VALUE o0_offset, VALUE o0_stride,
  VALUE o1, VALUE o1_offset, VALUE o1_stride,
  VALUE n)
{
  VdspArrayParam _i0;
  VdspArrayParam _i1;
  VdspArrayParam _o0;
  VdspArrayParam _o1;

  array_param(&_i0, i0, i0_offset, i0_stride);
  array_param(&_i1, i1, i1_offset, i1_stride);
  array_param(&_o0, o0, o0_offset, o0_stride);
  array_param(&_o1, o1, o1_offset, o1_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vaddsubD(
    _i0.res->v.d+_i0.offset, _i0.stride,
    _i1.res->v.d+_i1.offset, _i1.stride,
    _o0.res->v.d+_o0.offset, _o0.stride,
    _o1.res->v.d+_o1.offset, _o1.stride,
    _n);

  return rb_assoc_new(o0, o1);
}

// c[i] = a[i] / b[i]
VALUE rb_double_vdiv(
  VALUE cls,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _b;
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_b, b, b_offset, b_stride);
  array_param(&_a, a, a_offset, a_stride);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vdivD(
    _b.res->v.d+_b.offset, _c.stride,
    _a.res->v.d+_a.offset, _a.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

// d[i] = (a[i] + b[i]) * c
VALUE rb_double_vasm(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  double _c = NUM2DBL(c); 
  array_param(&_d, d, d_offset, d_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vasmD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _b.stride,
    &_c,
    _d.res->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

// d[i] = (a[i] + b[i]) * c[i]
VALUE rb_double_vam(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);
  array_param(&_d, d, d_offset, d_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vamD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _c.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _d.res->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

// d[i] = (a[i] - b[i]) * c
VALUE rb_double_vsbsm(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  double _c = NUM2DBL(c); 
  array_param(&_d, d, d_offset, d_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vsbsmD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _b.stride,
    &_c,
    _d.res->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

// d[i] = (a[i] - b[i]) * c[i]
VALUE rb_double_vsbm(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);
  array_param(&_d, d, d_offset, d_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vsbmD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _c.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _d.res->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

// d[i] = (a[i] * b[i]) + c
VALUE rb_double_vmsa(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  double _c = NUM2DBL(c); 
  array_param(&_d, d, d_offset, d_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vmsaD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _b.stride,
    &_c,
    _d.res->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

// d[i] = (a[i] * b) + c[i]
VALUE rb_double_vsma(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  double _b = NUM2DBL(b); 
  array_param(&_c, c, c_offset, c_stride);
  array_param(&_d, d, d_offset, d_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vsmaD(
    _a.res->v.d+_a.offset, _a.stride,
    &_b,
    _c.res->v.d+_c.offset, _c.stride,
    _d.res->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

// d[i] = (a[i] * b[i]) + c[i]
VALUE rb_double_vma(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);
  array_param(&_d, d, d_offset, d_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vmaD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _c.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _d.res->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

// d[i] = (a[i] * b[i]) - c[i]
VALUE rb_double_vmsb(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);
  array_param(&_d, d, d_offset, d_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vmsbD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _c.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _d.res->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

// e[i] = (a[i] * b) + (c[i] * d)
VALUE rb_double_vsmsma(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE d,
  VALUE e, VALUE e_offset, VALUE e_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;
  VdspArrayParam _e;

  array_param(&_a, a, a_offset, a_stride);
  double _b = NUM2DBL(b);
  array_param(&_c, c, c_offset, c_stride);
  double _d = NUM2DBL(d);
  array_param(&_e, e, e_offset, e_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vsmsmaD(
    _a.res->v.d+_a.offset, _a.stride,
    &_b,
    _c.res->v.d+_c.offset, _c.stride,
    &_d,
    _e.res->v.d+_e.offset, _e.stride,
    _n);

  return e;
}

// e[i] = (a[i] + b[i]) * (c[i] + d[i])
VALUE rb_double_vaam(int argc, VALUE *argv, VALUE cls)
{
  if (argc!=16) {
    rb_raise(rb_eArgError, "wrong number of arguments (given %d, expected 16)", argc);
  }
  VALUE a = argv[0];
  VALUE a_offset = argv[1];
  VALUE a_stride = argv[2];
  VALUE b = argv[3];
  VALUE b_offset = argv[4];
  VALUE b_stride = argv[5];
  VALUE c = argv[6];
  VALUE c_offset = argv[7];
  VALUE c_stride = argv[8];
  VALUE d = argv[9];
  VALUE d_offset = argv[10];
  VALUE d_stride = argv[11];
  VALUE e = argv[12];
  VALUE e_offset = argv[13];
  VALUE e_stride = argv[14];
  VALUE n = argv[15];

  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;
  VdspArrayParam _d;
  VdspArrayParam _e;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);
  array_param(&_d, d, d_offset, d_stride);
  array_param(&_e, e, e_offset, e_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vaamD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _b.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _d.res->v.d+_d.offset, _d.stride,
    _e.res->v.d+_e.offset, _e.stride,
    _n);

  return e;
}

// e[i] = (a[i] * b[i]) - (c[i] * d[i])
VALUE rb_double_vmmsb(int argc, VALUE *argv, VALUE cls)
{
  if (argc!=16) {
    rb_raise(rb_eArgError, "wrong number of arguments (given %d, expected 16)", argc);
  }
  VALUE a = argv[0];
  VALUE a_offset = argv[1];
  VALUE a_stride = argv[2];
  VALUE b = argv[3];
  VALUE b_offset = argv[4];
  VALUE b_stride = argv[5];
  VALUE c = argv[6];
  VALUE c_offset = argv[7];
  VALUE c_stride = argv[8];
  VALUE d = argv[9];
  VALUE d_offset = argv[10];
  VALUE d_stride = argv[11];
  VALUE e = argv[12];
  VALUE e_offset = argv[13];
  VALUE e_stride = argv[14];
  VALUE n = argv[15];

  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;
  VdspArrayParam _d;
  VdspArrayParam _e;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);
  array_param(&_d, d, d_offset, d_stride);
  array_param(&_e, e, e_offset, e_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vmmsbD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _b.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _d.res->v.d+_d.offset, _d.stride,
    _e.res->v.d+_e.offset, _e.stride,
    _n);

  return e;
}

// e[i] = (a[i] - b[i]) * (c[i] - d[i])
VALUE rb_double_vsbsbm(int argc, VALUE *argv, VALUE cls)
{
  if (argc!=16) {
    rb_raise(rb_eArgError, "wrong number of arguments (given %d, expected 16)", argc);
  }
  VALUE a = argv[0];
  VALUE a_offset = argv[1];
  VALUE a_stride = argv[2];
  VALUE b = argv[3];
  VALUE b_offset = argv[4];
  VALUE b_stride = argv[5];
  VALUE c = argv[6];
  VALUE c_offset = argv[7];
  VALUE c_stride = argv[8];
  VALUE d = argv[9];
  VALUE d_offset = argv[10];
  VALUE d_stride = argv[11];
  VALUE e = argv[12];
  VALUE e_offset = argv[13];
  VALUE e_stride = argv[14];
  VALUE n = argv[15];

  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;
  VdspArrayParam _d;
  VdspArrayParam _e;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);
  array_param(&_d, d, d_offset, d_stride);
  array_param(&_e, e, e_offset, e_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vsbsbmD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _b.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _d.res->v.d+_d.offset, _d.stride,
    _e.res->v.d+_e.offset, _e.stride,
    _n);

  return e;
}

// e[i] = (a[i] + b[i]) * (c[i] + d[i])
VALUE rb_double_vasbm(int argc, VALUE *argv, VALUE cls)
{
  if (argc!=16) {
    rb_raise(rb_eArgError, "wrong number of arguments (given %d, expected 16)", argc);
  }
  VALUE a = argv[0];
  VALUE a_offset = argv[1];
  VALUE a_stride = argv[2];
  VALUE b = argv[3];
  VALUE b_offset = argv[4];
  VALUE b_stride = argv[5];
  VALUE c = argv[6];
  VALUE c_offset = argv[7];
  VALUE c_stride = argv[8];
  VALUE d = argv[9];
  VALUE d_offset = argv[10];
  VALUE d_stride = argv[11];
  VALUE e = argv[12];
  VALUE e_offset = argv[13];
  VALUE e_stride = argv[14];
  VALUE n = argv[15];

  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;
  VdspArrayParam _d;
  VdspArrayParam _e;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);
  array_param(&_d, d, d_offset, d_stride);
  array_param(&_e, e, e_offset, e_stride);

  vDSP_Stride _n = NUM2LONG(n);

  vDSP_vasbmD(
    _a.res->v.d+_a.offset, _a.stride,
    _b.res->v.d+_b.offset, _b.stride,
    _c.res->v.d+_c.offset, _c.stride,
    _d.res->v.d+_d.offset, _d.stride,
    _e.res->v.d+_e.offset, _e.stride,
    _n);

  return e;
}


// Init

void Init_vdsp()
{
  // Vdsp
  rb_mVdsp = rb_define_module("Vdsp");

  // Vdsp::Scalar
  rb_mVdspScalar = rb_define_module_under(rb_mVdsp, "Scalar");

  //Vdsp::DoubleScalar
  rb_cDoubleScalar = rb_define_class_under(rb_mVdsp, "DoubleScalar", rb_cObject);
  rb_include_module(rb_cDoubleScalar, rb_mVdspScalar);
  rb_define_private_method(rb_cDoubleScalar, "initialize", rb_double_scalar_initialize, 1);
  rb_define_method(rb_cDoubleScalar, "+", rb_double_scalar_plus, 1);
  rb_define_method(rb_cDoubleScalar, "-", rb_double_scalar_minus, 1);
  rb_define_method(rb_cDoubleScalar, "*", rb_double_scalar_mul, 1);
  rb_define_method(rb_cDoubleScalar, "/", rb_double_scalar_div, 1);

  // Vdsp::Array
  rb_mVdspArray = rb_define_module_under(rb_mVdsp, "Array");
  rb_define_method(rb_mVdspArray, "length", rb_vdsp_array_length, 0);
  rb_define_alias(rb_mVdspArray,  "size", "length");

  // Vdsp::DoubleArray
  rb_cDoubleArray = rb_define_class_under(rb_mVdsp, "DoubleArray", rb_cObject);
  rb_include_module(rb_cDoubleArray, rb_mVdspArray);
  rb_define_method(rb_cDoubleArray, "initialize", rb_double_array_initialize, 1);
  rb_define_singleton_method(rb_cDoubleArray, "create", rb_double_array_create, 1);
  rb_define_method(rb_cDoubleArray, "to_da", rb_double_array_to_da, 0);
  rb_define_method(rb_cDoubleArray, "+", rb_double_array_plus, 1);
  rb_define_method(rb_cDoubleArray, "-", rb_double_array_minus, 1);
  rb_define_method(rb_cDoubleArray, "*", rb_double_array_mul, 1);
  rb_define_method(rb_cDoubleArray, "/", rb_double_array_div, 1);
  rb_define_method(rb_cDoubleArray, "[]", rb_double_array_aref, 1);
  rb_define_method(rb_cDoubleArray, "[]=", rb_double_array_aset, 2);
  rb_define_method(rb_cDoubleArray, "to_a", rb_double_array_get_values, 0);
  rb_define_method(rb_cDoubleArray, "coerce", rb_double_array_coerce, 1);

  // Vdsp::Double
  rb_mDouble = rb_define_module_under(rb_mVdsp, "Double");

  // Vector-based Arithmetic
  rb_define_singleton_method(rb_mDouble, "vsadd", rb_double_vsadd, 8);
  rb_define_singleton_method(rb_mDouble, "vadd", rb_double_vadd, 10);
  rb_define_singleton_method(rb_mDouble, "vsub", rb_double_vsub, 10);
  rb_define_singleton_method(rb_mDouble, "vsmul", rb_double_vsmul, 8);
  rb_define_singleton_method(rb_mDouble, "vmul", rb_double_vmul, 10);
  rb_define_singleton_method(rb_mDouble, "vsdiv", rb_double_vsdiv, 8);
  rb_define_singleton_method(rb_mDouble, "svdiv", rb_double_svdiv, 8);
  rb_define_singleton_method(rb_mDouble, "vaddsub", rb_double_vaddsub, 13);
  rb_define_singleton_method(rb_mDouble, "vdiv", rb_double_vdiv, 10);
  rb_define_singleton_method(rb_mDouble, "vasm", rb_double_vasm, 11);
  rb_define_singleton_method(rb_mDouble, "vam", rb_double_vam, 13);
  rb_define_singleton_method(rb_mDouble, "vsbsm", rb_double_vsbsm, 11);
  rb_define_singleton_method(rb_mDouble, "vsbm", rb_double_vsbm, 13);
  rb_define_singleton_method(rb_mDouble, "vmsa", rb_double_vmsa, 11);
  rb_define_singleton_method(rb_mDouble, "vsma", rb_double_vsma, 11);
  rb_define_singleton_method(rb_mDouble, "vma", rb_double_vma, 13);
  rb_define_singleton_method(rb_mDouble, "vmsb", rb_double_vmsb, 13);
  rb_define_singleton_method(rb_mDouble, "vsmsma", rb_double_vsmsma, 12);
  rb_define_singleton_method(rb_mDouble, "vaam", rb_double_vaam, -1);
  rb_define_singleton_method(rb_mDouble, "vmmsb", rb_double_vmmsb, -1);
  rb_define_singleton_method(rb_mDouble, "vsbsbm", rb_double_vsbsbm, -1);
  rb_define_singleton_method(rb_mDouble, "vasbm", rb_double_vasbm, -1);
}
