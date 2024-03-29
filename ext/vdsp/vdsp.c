#include "vdsp.h"


VALUE rb_mVdsp;

VALUE rb_mVdspScalar;
VALUE rb_mVdspArray;
VALUE rb_mVdspBiquad;
VALUE rb_mVdspFFT;

VALUE rb_cDoubleScalar;
VALUE rb_cDoubleArray;
VALUE rb_mUnsafeDouble;
VALUE rb_cDoubleBiquad;
VALUE rb_cDoubleFFT;


// Array native resource

void vdsp_array_native_resource_delete(VdspArrayNativeResource *p)
{
  if (p->v.ptr) {
    free(p->v.ptr);
    p->v.ptr = NULL;
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
    rb_raise(rb_eArgError, "Vdsp::Array has no native resource");
  }

  VdspArrayNativeResource *p;
  Data_Get_Struct(resource, VdspArrayNativeResource, p);

  return p;
}

void array_param(VdspArrayParam *param, VALUE arr0, VALUE offset, VALUE stride)
{
  param->res0 = get_vdsp_array_native_resource(arr0);
  param->offset = NUM2LONG(offset);
  param->stride = NUM2LONG(stride);
}

void array_param2(VdspArrayParam *param, VALUE arr0, VALUE arr1, VALUE offset, VALUE stride)
{
  param->res0 = get_vdsp_array_native_resource(arr0);
  param->res1 = get_vdsp_array_native_resource(arr1);
  param->offset = NUM2LONG(offset);
  param->stride = NUM2LONG(stride);
}


// Biquad native resource

void vdsp_biquad_native_resource_delete(VdspBiquadNativeResource *p)
{
  if (p->setup.ptr) {
    if (p->type=='d') {
      vDSP_biquad_DestroySetupD(p->setup.d);
      p->setup.ptr = NULL;
    }
  }
  if (p->delay.ptr) {
    free(p->delay.ptr);
    p->delay.ptr = NULL;
  }
  free(p);
}

VdspBiquadNativeResource* get_vdsp_biquad_native_resource(VALUE vb)
{
  if (!rb_obj_is_kind_of(vb, rb_mVdspBiquad)) {
    rb_raise(rb_eArgError, "Vdsp::Biquad required");
  }

  VALUE resource = rb_iv_get(vb, "native_resource");
  if (resource==Qnil) {
    rb_raise(rb_eArgError, "Vdsp::Biquad has no native resource");
  }

  VdspBiquadNativeResource *p;
  Data_Get_Struct(resource, VdspBiquadNativeResource, p);

  return p;
}


// FFT native resource

void vdsp_fft_native_resource_delete(VdspFFTNativeResource *p)
{
  if (p->setup.value) {
    if (p->type=='d') {
      vDSP_destroy_fftsetupD(p->setup.d);
    }
  }
  free(p);
}

VdspFFTNativeResource* get_vdsp_fft_native_resource(VALUE vf)
{
  if (!rb_obj_is_kind_of(vf, rb_mVdspFFT)) {
    rb_raise(rb_eArgError, "Vdsp::FFT required");
  }

  VALUE resource = rb_iv_get(vf, "native_resource");
  if (resource==Qnil) {
    rb_raise(rb_eArgError, "Vdsp::FFT has no native resource");
  }

  VdspFFTNativeResource *p;
  Data_Get_Struct(resource, VdspFFTNativeResource, p);

  return p;
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
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  return LONG2NUM(p->length);
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

VALUE rb_double_array_initialize_copy(VALUE self, VALUE orig)
{
  VALUE length = rb_vdsp_array_length(orig);
  rb_double_array_initialize(self, length);

  VdspArrayNativeResource *dst = get_vdsp_array_native_resource(self);
  VdspArrayNativeResource *src = get_vdsp_array_native_resource(orig);
  memcpy(dst->v.ptr, src->v.ptr, sizeof(double) * src->length);

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

VALUE rb_double_array_each(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double *d = p->v.d;

  RETURN_ENUMERATOR(self, 0, 0);

  for (unsigned long i=0; i<p->length; i++) {
      rb_yield(DBL2NUM(d[i]));
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

VALUE rb_double_array_entry(VALUE self, long offset)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
  long len = _a->length;

  if (offset < 0) {
    offset += len;
    if (offset < 0) return Qnil;
  }
  else if (len <= offset) {
    return Qnil;
  }
  return DBL2NUM(_a->v.d[offset]);
}

VALUE rb_double_array_subseq(VALUE self, long beg, long len)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
  long alen = _a->length;

  if (beg > alen) return Qnil;
  if (beg < 0 || len < 0) return Qnil;

  if (alen < len || alen < beg + len) {
    len = alen - beg;
  }

  VALUE lenv = LONG2NUM(len);
  VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);
  memcpy(_c->v.d, _a->v.d+beg, sizeof(double) * len);

  return c;
}

VALUE rb_double_array_aref2(VALUE self, VALUE b, VALUE e)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);

  long beg = NUM2LONG(b);
  long len = NUM2LONG(e);

  if (beg < 0) {
    beg += _a->length;
  }

  return rb_double_array_subseq(self, beg, len);
}

VALUE rb_double_array_aref1(VALUE self, VALUE arg)
{
  long beg, len;

  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
  long alen = _a->length;

  if (FIXNUM_P(arg)) {
    return rb_double_array_entry(self, FIX2LONG(arg));
  }
  switch (rb_range_beg_len(arg, &beg, &len, alen, 0)) {
    case Qfalse:
      break;
    case Qnil:
      return Qnil;
    default:
      return rb_double_array_subseq(self, beg, len);
  }
  return rb_double_array_entry(self, FIX2LONG(arg));
}

VALUE rb_double_array_aref(int argc, const VALUE *argv, VALUE self)
{
  rb_check_arity(argc, 1, 2);

  if (argc == 2) {
    return rb_double_array_aref2(self, argv[0], argv[1]);
  }
  return rb_double_array_aref1(self, argv[0]);
}

VALUE double_array_aset1(VALUE self, long idx, VALUE val)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
  double _val = NUM2DBL(val);

  long len = _a->length;

  if (idx < 0) {
    idx += len;
    if (idx < 0) {
      rb_raise(rb_eIndexError, "index %ld too small for array; minimum: %ld", idx - len, -len);
    }
  }

  if (idx >= len) {
    //rb_raise(rb_eIndexError, "Index out of range: %ld", idx);
    double_array_resize(_a, idx+1);
    _a->length = idx+1;
  }

  _a->v.d[idx] = _val;

  return self;
}

VALUE double_array_aset2(VALUE self, long beg, long len, VALUE val)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
  long alen = _a->length;

  if (len<0) {
    len += alen;
  }
  if (len<beg) {
    len = beg;
  }
  if (alen < len || alen < beg + len) {
    len = alen - beg;
  }

  if (RB_TYPE_P(val, T_ARRAY)) {
    val = rb_funcall(rb_cDoubleArray, rb_intern("create"), 1, val);
  }
  if (rb_obj_is_kind_of(val, rb_mVdspArray)) {
    val = rb_funcall(val, rb_intern("to_da"), 0);
    if (self==val) {
      val = rb_funcall(val, rb_intern("clone"), 0);
    }

    VdspArrayNativeResource *_b = get_vdsp_array_native_resource(val);
    long new_len = alen - len + _b->length;

    bool after_resize = false;
    bool rear_move = true;
    if (alen==new_len) {
      rear_move = false;
    } else if (alen<new_len) {
      double_array_resize(_a, new_len);
      _a->length = new_len;
    } else {
      after_resize = true;
    }

    if (rear_move) {
      memmove(_a->v.d+beg+_b->length, _a->v.d+beg+len, sizeof(double) * (alen-beg-len));
    }
    memcpy(_a->v.d+beg, _b->v.d, sizeof(double) * _b->length);

    if (after_resize) {
      double_array_resize(_a, new_len);
      _a->length = new_len;
    }

    return self;
  } else {
    long new_len = alen - len + 1;

    bool after_resize = false;
    if (alen==new_len) {
    } else if (alen<new_len) {
      double_array_resize(_a, new_len);
      _a->length = new_len;
    } else {
      after_resize = true;
    }
    
    memmove(_a->v.d+beg+1, _a->v.d+len, sizeof(double) * (alen-beg-len));
    _a->v.d[beg] = NUM2DBL(val);

    if (after_resize) {
      double_array_resize(_a, new_len);
      _a->length = new_len;
    }

    return self;
  }
}

VALUE rb_double_array_aset(int argc, const VALUE *argv, VALUE self)
{
  long offset, beg, len;
  VALUE val;

  if (argc == 3) {
    beg = NUM2LONG(argv[0]);
    len = NUM2LONG(argv[1]);
    val = argv[2];
    return double_array_aset2(self, beg, len, val);
  }

  rb_check_arity(argc, 2, 2);

  if (FIXNUM_P(argv[0])) {
    offset = FIX2LONG(argv[0]);
    val = argv[1];
    return double_array_aset1(self, offset, val);
  }

  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
  if (rb_range_beg_len(argv[0], &beg, &len, _a->length, 1)) {
    val = argv[1];
    return double_array_aset2(self, beg, len, val);
  }

  offset = NUM2LONG(argv[0]);
  val = argv[1];
  return double_array_aset1(self, offset, val);
}

VALUE rb_double_array_coerce(VALUE self, VALUE other)
{
  other = rb_class_new_instance(1, &other, rb_cDoubleScalar);
  return rb_assoc_new(other, self);
}

void double_array_resize(VdspArrayNativeResource *_a, unsigned long len)
{
  if (_a->length!=len) {
    void *new_v = realloc(_a->v.ptr, (_a->length + len) * sizeof(double));
    if (!new_v) {
      rb_raise(rb_eRuntimeError, "memory allocation error");
    }
    _a->v.ptr = new_v;

    if (_a->length<len) {
      memset(_a->v.d+_a->length, 0, (len-_a->length) * sizeof(double));
    }
  }
}

VALUE rb_double_array_resize(VALUE self, VALUE size)
{
  unsigned long len = NUM2LONG(size);

  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
  double_array_resize(_a, len);
  _a->length = len;

  return self;
}

VALUE rb_double_array_concat(int argc, const VALUE *argv, VALUE self)
{
  rb_check_frozen(self);
  rb_check_arity(argc, 1, UNLIMITED_ARGUMENTS);

  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);
  VdspArrayNativeResource *_argv[argc];

  unsigned long len = _a->length;
  unsigned long new_len = len;
  for (long i=0; i<argc; i++) {
    VALUE arg = rb_funcall(argv[i], rb_intern("to_da"), 0);
    _argv[i] = get_vdsp_array_native_resource(arg);
    new_len += _argv[i]->length;
  }

  double_array_resize(_a, new_len);

  new_len = len;
  for (long i=0; i<argc; i++) {
    VdspArrayNativeResource *_b = _argv[i];
    memcpy(_a->v.d+new_len, _b->v.d, _b->length * sizeof(double));
    new_len += _b->length;
  }
  _a->length = new_len;

  return self;
}

VALUE rb_double_array_abs(VALUE self)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);

  VALUE lenv = LONG2NUM(_a->length);
  VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

  vDSP_vabsD(_a->v.d, 1, _c->v.d, 1, _a->length);

  return c;
}

VALUE rb_double_array_abs_bang(VALUE self)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);

  vDSP_vabsD(_a->v.d, 1, _a->v.d, 1, _a->length);

  return self;
}

VALUE rb_double_array_nabs(VALUE self)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);

  VALUE lenv = LONG2NUM(_a->length);
  VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

  vDSP_vnabsD(_a->v.d, 1, _c->v.d, 1, _a->length);

  return c;
}

VALUE rb_double_array_nabs_bang(VALUE self)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);

  vDSP_vnabsD(_a->v.d, 1, _a->v.d, 1, _a->length);

  return self;
}

VALUE rb_double_array_negative(VALUE self)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);

  VALUE lenv = LONG2NUM(_a->length);
  VALUE c = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

  vDSP_vnegD(_a->v.d, 1, _c->v.d, 1, _a->length);

  return c;
}

VALUE rb_double_array_negative_bang(VALUE self)
{
  VdspArrayNativeResource *_a = get_vdsp_array_native_resource(self);

  vDSP_vnegD(_a->v.d, 1, _a->v.d, 1, _a->length);

  return self;
}

VALUE rb_double_array_vramp(VALUE cls, VALUE a, VALUE b, VALUE n)
{
  double _a = NUM2DBL(a);
  double _b = NUM2DBL(b);
  vDSP_Length _n = NUM2LONG(n);

  VALUE c = rb_class_new_instance(1, &n, rb_cDoubleArray);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

  vDSP_vrampD(&_a, &_b, _c->v.d, 1, _n);

  return c;
}

VALUE rb_double_array_vgen(VALUE cls, VALUE a, VALUE b, VALUE n)
{
  double _a = NUM2DBL(a);
  double _b = NUM2DBL(b);
  vDSP_Length _n = NUM2LONG(n);

  VALUE c = rb_class_new_instance(1, &n, rb_cDoubleArray);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

  vDSP_vgenD(&_a, &_b, _c->v.d, 1, _n);

  return c;
}

VALUE rb_double_array_blkman_window(int argc, const VALUE *argv, VALUE cls)
{
  rb_check_arity(argc, 1, 2);
  VALUE n = argv[0];

  vDSP_Length _n = NUM2LONG(n);
  int flag = 0;
  if (argc==2) {
    flag = (int)NUM2LONG(argv[1]);
  }

  VALUE c = rb_class_new_instance(1, &n, rb_cDoubleArray);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

  vDSP_blkman_windowD(_c->v.d, _n, flag);

  return c;
}

VALUE rb_double_array_hamm_window(int argc, const VALUE *argv, VALUE cls)
{
  rb_check_arity(argc, 1, 2);
  VALUE n = argv[0];

  vDSP_Length _n = NUM2LONG(n);
  int flag = 0;
  if (argc==2) {
    flag = (int)NUM2LONG(argv[1]);
  }

  VALUE c = rb_class_new_instance(1, &n, rb_cDoubleArray);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

  vDSP_hamm_windowD(_c->v.d, _n, flag);

  return c;
}

VALUE rb_double_array_hann_window(int argc, const VALUE *argv, VALUE cls)
{
  rb_check_arity(argc, 1, 2);
  VALUE n = argv[0];

  vDSP_Length _n = NUM2LONG(n);
  int flag = 0;
  if (argc==2) {
    flag = (int)NUM2LONG(argv[1]);
  }

  VALUE c = rb_class_new_instance(1, &n, rb_cDoubleArray);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

  vDSP_hann_windowD(_c->v.d, _n, flag);

  return c;
}

VALUE rb_double_array_clear_bang(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  vDSP_vclrD(p->v.d, 1, p->length);
  return self;
}

VALUE rb_double_array_fill_bang(VALUE self, VALUE a)
{
  double _a = NUM2DBL(a);
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  vDSP_vfillD(&_a, p->v.d, 1, p->length);
  return self;
}

VALUE rb_double_array_maxv(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_maxvD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_maxmgv(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_maxmgvD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_minv(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_minvD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_minmgv(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_minmgvD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_meanv(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_meanvD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_meamgv(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_meamgvD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_measqv(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_measqvD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_mvessq(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_mvessqD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_rmsqv(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_rmsqvD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_sve(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_sveD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_svemg(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_svemgD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_svesq(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_svesqD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}

VALUE rb_double_array_sve_svesq(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _sum;
  double _sum_of_squares;
  vDSP_sve_svesqD(p->v.d, 1, &_sum, &_sum_of_squares, p->length);
  return rb_assoc_new(DBL2NUM(_sum), DBL2NUM(_sum_of_squares));
}

VALUE rb_double_array_svs(VALUE self)
{
  VdspArrayNativeResource *p = get_vdsp_array_native_resource(self);
  double _c;
  vDSP_svsD(p->v.d, 1, &_c, p->length);
  return DBL2NUM(_c);
}


// Vdsp::Biquad

VALUE rb_vdsp_biquad_sections(VALUE self)
{
  VdspBiquadNativeResource *p = get_vdsp_biquad_native_resource(self);
  return LONG2NUM(p->sections);
}

VALUE rb_vdsp_biquad_alloc_sections(VALUE self)
{
  VdspBiquadNativeResource *p = get_vdsp_biquad_native_resource(self);
  return LONG2NUM(p->alloc_sections);
}


// Vdsp::DoubleBiquad

VALUE rb_double_biquad_initialize(VALUE self, VALUE alloc_sections)
{
  VdspBiquadNativeResource *p = ALLOC(VdspBiquadNativeResource);
  p->type = 'd';
  p->coefs.ptr = NULL;
  p->delay.ptr = NULL;
  p->setup.ptr = NULL;
  p->sections = 0;
  p->alloc_sections = 0;

  VALUE resource = Data_Wrap_Struct(CLASS_OF(self), 0, vdsp_biquad_native_resource_delete, p);
  rb_iv_set(self, "native_resource", resource);

  long _alloc_sections = NUM2LONG(alloc_sections);
  if (_alloc_sections<1) {
    _alloc_sections = 1;
  }
  p->coefs.ptr = calloc(_alloc_sections*5, sizeof(double));
  p->delay.ptr = calloc(_alloc_sections*2+2, sizeof(double));
  p->alloc_sections = _alloc_sections;

  return self;
}

VALUE rb_double_biquad_set_coefficients(VALUE self, VALUE coefficients)
{
  VdspBiquadNativeResource *p = get_vdsp_biquad_native_resource(self);

  coefficients = rb_ary_new3(1, coefficients);
  coefficients = rb_funcall(coefficients, rb_intern("flatten"), 0);

  VALUE rb_cBiquadCoefficient = rb_const_get(rb_mVdspBiquad, rb_intern("Coefficient"));

  unsigned long sections = RARRAY_LEN(coefficients);
  if (p->alloc_sections<sections) {
    rb_raise(rb_eArgError, "over sections: sections=%ld alloc_sections=%ld", sections, p->alloc_sections);
  }

  for (unsigned long i=0; i<sections; i++) {
    VALUE coefficient = RARRAY_AREF(coefficients, i);
    if (!rb_obj_is_kind_of(coefficient, rb_cBiquadCoefficient)) {
      rb_raise(rb_eArgError, "Vdsp::Biquad::Coefficient required");
    }
  }

  for (unsigned long i=0; i<sections; i++) {
    VALUE coefficient = RARRAY_AREF(coefficients, i);
    p->coefs.d[i*5+0] = NUM2DBL(rb_funcall(coefficient, rb_intern("b0"), 0));
    p->coefs.d[i*5+1] = NUM2DBL(rb_funcall(coefficient, rb_intern("b1"), 0));
    p->coefs.d[i*5+2] = NUM2DBL(rb_funcall(coefficient, rb_intern("b2"), 0));
    p->coefs.d[i*5+3] = NUM2DBL(rb_funcall(coefficient, rb_intern("a1"), 0));
    p->coefs.d[i*5+4] = NUM2DBL(rb_funcall(coefficient, rb_intern("a2"), 0));
  }

  if (p->setup.d) {
    vDSP_biquad_DestroySetupD(p->setup.d);
  }
  p->setup.d = vDSP_biquad_CreateSetupD(p->coefs.d, sections);
  p->sections = sections;

  return self;
}

VALUE rb_double_biquad_get_coefficients(VALUE self)
{
  VdspBiquadNativeResource *p = get_vdsp_biquad_native_resource(self);
  VALUE ret = rb_ary_new2(p->sections);

  VALUE rb_cBiquadCoefficient = rb_const_get(rb_mVdspBiquad, rb_intern("Coefficient"));

  for (unsigned long i=0; i<p->sections; i++) {
    VALUE argv[5] = {
      DBL2NUM(p->coefs.d[i*5+0]),
      DBL2NUM(p->coefs.d[i*5+1]),
      DBL2NUM(p->coefs.d[i*5+2]),
      DBL2NUM(p->coefs.d[i*5+3]),
      DBL2NUM(p->coefs.d[i*5+4])
    };
    VALUE coef = rb_class_new_instance(5, argv, rb_cBiquadCoefficient);
    rb_ary_push(ret, coef);
  }

  return ret;
}

VALUE rb_double_biquad_create(VALUE cls, VALUE coefficients)
{
  coefficients = rb_ary_new3(1, coefficients);
  coefficients = rb_funcall(coefficients, rb_intern("flatten"), 0);

  VALUE len = LONG2NUM(RARRAY_LEN(coefficients));
  VALUE obj = rb_class_new_instance(1, &len, rb_cDoubleBiquad);
  rb_double_biquad_set_coefficients(obj, coefficients);

  return obj;
}

VALUE rb_double_biquad_apply(VALUE self, VALUE x)
{
  //x = rb_funcall(x, rb_intern("to_da"), 0);
  x = rb_funcall(rb_cDoubleArray, rb_intern("create"), 1, x);
  VdspArrayNativeResource *_x = get_vdsp_array_native_resource(x);

  VALUE lenv = LONG2NUM(_x->length);
  VALUE y = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VdspArrayNativeResource *_y = get_vdsp_array_native_resource(y);

  VdspBiquadNativeResource *_b = get_vdsp_biquad_native_resource(self);

  vDSP_biquadD(_b->setup.d, _b->delay.d, _x->v.d, 1, _y->v.d, 1, _x->length);

  return y;
}


// Vdsp::DoubleFFT

VALUE rb_double_fft_initialize(int argc, const VALUE *argv, VALUE self)
{
  rb_check_arity(argc, 1, 2);
  long length = NUM2LONG(argv[0]);
  FFTRadix radix = kFFTRadix2;
  if (argc==2) {
    radix = (FFTRadix)NUM2LONG(argv[1]);
  }

  VdspFFTNativeResource *p = ALLOC(VdspFFTNativeResource);
  p->type = 'd';
  p->length = 0;
  p->log2n = 0;
  p->radix = 0;
  p->setup.value = NULL;

  VALUE resource = Data_Wrap_Struct(CLASS_OF(self), 0, vdsp_fft_native_resource_delete, p);
  rb_iv_set(self, "native_resource", resource);

  p->length = length;
  p->halflength = length / 2;
  p->log2n = log2(length);
  p->radix = radix;
  p->setup.d = vDSP_create_fftsetupD(p->log2n, radix);

  return self;
}

VALUE double_fft_forward_r(VALUE fft, VALUE va)
{
  VdspFFTNativeResource *_vf = get_vdsp_fft_native_resource(fft);
  VdspArrayNativeResource *_va = get_vdsp_array_native_resource(va);

  if (_vf->length!=_va->length) {
    rb_raise(rb_eArgError, "wrong length: Vdsp::FFT=%ld Vdsp::Array(real)=%ld", _vf->length, _va->length);
  }

  VALUE lenv = LONG2NUM(_vf->halflength);
  VALUE real = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VALUE imag = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VdspArrayNativeResource *_real = get_vdsp_array_native_resource(real);
  VdspArrayNativeResource *_imag = get_vdsp_array_native_resource(imag);

  DSPDoubleSplitComplex z;
  z.realp = _real->v.d;
  z.imagp = _imag->v.d;

  vDSP_ctozD((DSPDoubleComplex *)_va->v.d, 2, &z, 1, _vf->halflength);
  vDSP_fft_zripD(_vf->setup.d, &z, 1, _vf->log2n, FFT_FORWARD);

  double scale = 1.0 / (_vf->length * 2);
  vDSP_vsmulD(z.realp, 1, &scale, z.realp, 1, _vf->halflength);
  vDSP_vsmulD(z.imagp, 1, &scale, z.imagp, 1, _vf->halflength);

  return rb_assoc_new(real, imag);
}

VALUE double_fft_forward_z(VALUE fft, VALUE real, VALUE imag)
{
  VdspFFTNativeResource *_vf = get_vdsp_fft_native_resource(fft);
  VdspArrayNativeResource *_real = get_vdsp_array_native_resource(real);
  VdspArrayNativeResource *_imag = get_vdsp_array_native_resource(imag);

  if (_vf->length!=_real->length || _real->length!=_imag->length) {
    rb_raise(rb_eArgError, "wrong length: Vdsp::FFT=%ld Vdsp::Array(real)=%ld Vdsp::Array(imag)=%ld", _vf->length, _real->length, _imag->length);
  }

  real = rb_funcall(real, rb_intern("clone"), 0);
  imag = rb_funcall(imag, rb_intern("clone"), 0);
  _real = get_vdsp_array_native_resource(real);
  _imag = get_vdsp_array_native_resource(imag);

  DSPDoubleSplitComplex z;
  z.realp = _real->v.d;
  z.imagp = _imag->v.d;

  vDSP_fft_zipD(_vf->setup.d, &z, 1, _vf->log2n, FFT_FORWARD);

  double scale = 1.0 / _vf->length;
  vDSP_vsmulD(z.realp, 1, &scale, z.realp, 1, _vf->length);
  vDSP_vsmulD(z.imagp, 1, &scale, z.imagp, 1, _vf->length);

  return rb_assoc_new(real, imag);
}

VALUE rb_double_fft_forward(int argc, const VALUE *argv, VALUE self)
{
  rb_check_arity(argc, 1, 2);
  if (argc==2 && argv[1]!=Qnil) {
    return double_fft_forward_z(self, argv[0], argv[1]);
  }
  return double_fft_forward_r(self, argv[0]);
}

VALUE double_fft_inverse_r(VALUE fft, VALUE real, VALUE imag)
{
  real = rb_funcall(real, rb_intern("clone"), 0);
  imag = rb_funcall(imag, rb_intern("clone"), 0);

  VdspFFTNativeResource *_vf = get_vdsp_fft_native_resource(fft);
  VdspArrayNativeResource *_real = get_vdsp_array_native_resource(real);
  VdspArrayNativeResource *_imag = get_vdsp_array_native_resource(imag);

  DSPDoubleSplitComplex z;
  z.realp = _real->v.d;
  z.imagp = _imag->v.d;

  vDSP_fft_zripD(_vf->setup.d, &z, 1, _vf->log2n, FFT_INVERSE);

  VALUE lenv = LONG2NUM(_vf->length);
  VALUE inverse = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VdspArrayNativeResource *_inverse = get_vdsp_array_native_resource(inverse);

  vDSP_ztocD(&z, 1, (DSPDoubleComplex *)_inverse->v.d, 2, _vf->halflength);

  return rb_assoc_new(inverse, Qnil);
}

VALUE double_fft_inverse_z(VALUE fft, VALUE real, VALUE imag)
{
  real = rb_funcall(real, rb_intern("clone"), 0);
  imag = rb_funcall(imag, rb_intern("clone"), 0);

  VdspFFTNativeResource *_vf = get_vdsp_fft_native_resource(fft);
  VdspArrayNativeResource *_real = get_vdsp_array_native_resource(real);
  VdspArrayNativeResource *_imag = get_vdsp_array_native_resource(imag);

  DSPDoubleSplitComplex z;
  z.realp = _real->v.d;
  z.imagp = _imag->v.d;

  vDSP_fft_zipD(_vf->setup.d, &z, 1, _vf->log2n, FFT_INVERSE);

  return rb_assoc_new(real, imag);
}

VALUE rb_double_fft_inverse(VALUE self, VALUE real, VALUE imag)
{
  VdspFFTNativeResource *_vf = get_vdsp_fft_native_resource(self);
  VdspArrayNativeResource *_real = get_vdsp_array_native_resource(real);
  VdspArrayNativeResource *_imag = get_vdsp_array_native_resource(imag);

  if (_vf->length==_real->length) {
    return double_fft_inverse_z(self, real, imag);

  } else if (_vf->halflength==_real->length) {
    return double_fft_inverse_r(self, real, imag);

  } else {
    rb_raise(rb_eArgError, "wrong length: Vdsp::FFT=%ld Vdsp::Array(real)=%ld Vdsp::Array(imag)=%ld", _vf->length, _real->length, _imag->length);
  }
}

VALUE rb_double_fft_forward_polar(int argc, const VALUE *argv, VALUE self)
{
  VALUE ri = rb_double_fft_forward(argc, argv, self);
  VALUE real = RARRAY_AREF(ri, 0);
  VALUE imag = RARRAY_AREF(ri, 1);

  VdspArrayNativeResource *_real = get_vdsp_array_native_resource(real);
  VdspArrayNativeResource *_imag = get_vdsp_array_native_resource(imag);

  DSPDoubleSplitComplex z;
  z.realp = _real->v.d;
  z.imagp = _imag->v.d;

  VALUE lenv = LONG2NUM(_real->length*2);
  VALUE tmp = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VdspArrayNativeResource *_tmp = get_vdsp_array_native_resource(tmp);

  vDSP_ztocD(&z, 1, (DSPDoubleComplex *)_tmp->v.d, 2, _real->length);
  vDSP_polarD(_tmp->v.d, 2, _tmp->v.d, 2, _real->length);
  vDSP_ctozD((DSPDoubleComplex *)_tmp->v.d, 2, &z, 1, _real->length);

  return rb_assoc_new(real, imag);
}

VALUE rb_double_fft_inverse_polar(VALUE self, VALUE mag, VALUE pha)
{
  VdspArrayNativeResource *_mag = get_vdsp_array_native_resource(mag);
  VdspArrayNativeResource *_pha = get_vdsp_array_native_resource(pha);

  DSPDoubleSplitComplex z;
  z.realp = _mag->v.d;
  z.imagp = _pha->v.d;

  VALUE lenv = LONG2NUM(_mag->length*2);
  VALUE tmp = rb_class_new_instance(1, &lenv, rb_cDoubleArray);
  VdspArrayNativeResource *_tmp = get_vdsp_array_native_resource(tmp);

  vDSP_ztocD(&z, 1, (DSPDoubleComplex *)_tmp->v.d, 2, _mag->length);
  vDSP_rectD(_tmp->v.d, 2, _tmp->v.d, 2, _mag->length);
  vDSP_ctozD((DSPDoubleComplex *)_tmp->v.d, 2, &z, 1, _mag->length);

  return rb_double_fft_inverse(self, mag, pha);
}


// Vdsp static method

// c[i] = a[i]
VALUE rb_double_copy(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Length _n = NUM2LONG(n);

  if (_a.stride==1 && _c.stride==1) {
    memmove(_c.res0->v.d+_c.offset, _a.res0->v.d+_a.offset, sizeof(double) * _n);
  } else {
    double _b = 0;
    vDSP_vsaddD(
      _a.res0->v.d+_a.offset, _a.stride,
      &_b,
      _c.res0->v.d+_c.offset, _c.stride,
      _n);
  }

  return c;
}

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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vsaddD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    _c.res0->v.d+_c.offset, _c.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vaddD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    _c.res0->v.d+_c.offset, _c.stride,
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
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vsubD(
    _b.res0->v.d+_b.offset, _b.stride,
    _a.res0->v.d+_a.offset, _a.stride,
    _c.res0->v.d+_c.offset, _c.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vsmulD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    _c.res0->v.d+_c.offset, _c.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vmulD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    _c.res0->v.d+_c.offset, _c.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vsdivD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    _c.res0->v.d+_c.offset, _c.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_svdivD(
    &_a,
    _b.res0->v.d+_b.offset, _b.stride,
    _c.res0->v.d+_c.offset, _c.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vaddsubD(
    _i0.res0->v.d+_i0.offset, _i0.stride,
    _i1.res0->v.d+_i1.offset, _i1.stride,
    _o0.res0->v.d+_o0.offset, _o0.stride,
    _o1.res0->v.d+_o1.offset, _o1.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vdivD(
    _b.res0->v.d+_b.offset, _c.stride,
    _a.res0->v.d+_a.offset, _a.stride,
    _c.res0->v.d+_c.offset, _c.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vasmD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    &_c,
    _d.res0->v.d+_d.offset, _d.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vamD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _c.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _d.res0->v.d+_d.offset, _d.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vsbsmD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    &_c,
    _d.res0->v.d+_d.offset, _d.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vsbmD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _c.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _d.res0->v.d+_d.offset, _d.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vmsaD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    &_c,
    _d.res0->v.d+_d.offset, _d.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vsmaD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    _c.res0->v.d+_c.offset, _c.stride,
    _d.res0->v.d+_d.offset, _d.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vmaD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _c.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _d.res0->v.d+_d.offset, _d.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vmsbD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _c.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _d.res0->v.d+_d.offset, _d.stride,
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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vsmsmaD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    _c.res0->v.d+_c.offset, _c.stride,
    &_d,
    _e.res0->v.d+_e.offset, _e.stride,
    _n);

  return e;
}

// e[i] = (a[i] + b[i]) * (c[i] + d[i])
VALUE rb_double_vaam(int argc, const VALUE *argv, VALUE cls)
{
  rb_check_arity(argc, 16, 16);

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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vaamD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _d.res0->v.d+_d.offset, _d.stride,
    _e.res0->v.d+_e.offset, _e.stride,
    _n);

  return e;
}

// e[i] = (a[i] * b[i]) - (c[i] * d[i])
VALUE rb_double_vmmsb(int argc, const VALUE *argv, VALUE cls)
{
  rb_check_arity(argc, 16, 16);

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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vmmsbD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _d.res0->v.d+_d.offset, _d.stride,
    _e.res0->v.d+_e.offset, _e.stride,
    _n);

  return e;
}

// e[i] = (a[i] - b[i]) * (c[i] - d[i])
VALUE rb_double_vsbsbm(int argc, const VALUE *argv, VALUE cls)
{
  rb_check_arity(argc, 16, 16);

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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vsbsbmD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _d.res0->v.d+_d.offset, _d.stride,
    _e.res0->v.d+_e.offset, _e.stride,
    _n);

  return e;
}

// e[i] = (a[i] + b[i]) * (c[i] - d[i])
VALUE rb_double_vasbm(int argc, const VALUE *argv, VALUE cls)
{
  rb_check_arity(argc, 16, 16);

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

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vasbmD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _d.res0->v.d+_d.offset, _d.stride,
    _e.res0->v.d+_e.offset, _e.stride,
    _n);

  return e;
}

// Vector Generation

VALUE rb_double_vramp(
  VALUE cls,
  VALUE a,
  VALUE b,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _c;

  double _a = NUM2DBL(a);
  double _b = NUM2DBL(b);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vrampD(
    &_a,
    &_b,
    _c.res0->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

VALUE rb_double_vgen(
  VALUE cls,
  VALUE a,
  VALUE b,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _c;

  double _a = NUM2DBL(a);
  double _b = NUM2DBL(b);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vgenD(
    &_a,
    &_b,
    _c.res0->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

VALUE rb_double_vrampmul(
  VALUE cls,
  VALUE i, VALUE i_offset, VALUE i_stride,
  VALUE start,
  VALUE step,
  VALUE o, VALUE o_offset, VALUE o_stride,
  VALUE n)
{
  VdspArrayParam _i;
  VdspArrayParam _o;

  array_param(&_i, i, i_offset, i_stride);
  double _start = NUM2DBL(start);
  double _step = NUM2DBL(step);
  array_param(&_o, o, o_offset, o_stride);

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vrampmulD(
    _i.res0->v.d+_i.offset, _i.stride,
    &_start,
    &_step,
    _o.res0->v.d+_o.offset, _o.stride,
    _n
  );

  return o;
}

VALUE rb_double_vrampmul2(
  VALUE cls,
  VALUE i0, VALUE i1, VALUE i_offset, VALUE i_stride,
  VALUE start,
  VALUE step,
  VALUE o0, VALUE o1, VALUE o_offset, VALUE o_stride,
  VALUE n)
{
  VdspArrayParam _i;
  VdspArrayParam _o;

  array_param2(&_i, i0, i1, i_offset, i_stride);
  double _start = NUM2DBL(start);
  double _step = NUM2DBL(step);
  array_param2(&_o, o0, o1, o_offset, o_stride);

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vrampmul2D(
    _i.res0->v.d+_i.offset, _i.res1->v.d+_i.offset, _i.stride,
    &_start,
    &_step,
    _o.res0->v.d+_o.offset, _o.res1->v.d+_o.offset, _o.stride,
    _n
  );

  return rb_assoc_new(o0, o1);
}

VALUE rb_double_vrampmuladd(
  VALUE cls,
  VALUE i, VALUE i_offset, VALUE i_stride,
  VALUE start,
  VALUE step,
  VALUE o, VALUE o_offset, VALUE o_stride,
  VALUE n)
{
  VdspArrayParam _i;
  VdspArrayParam _o;

  array_param(&_i, i, i_offset, i_stride);
  double _start = NUM2DBL(start);
  double _step = NUM2DBL(step);
  array_param(&_o, o, o_offset, o_stride);

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vrampmuladdD(
    _i.res0->v.d+_i.offset, _i.stride,
    &_start,
    &_step,
    _o.res0->v.d+_o.offset, _o.stride,
    _n
  );

  return o;
}

VALUE rb_double_vrampmuladd2(
  VALUE cls,
  VALUE i0, VALUE i1, VALUE i_offset, VALUE i_stride,
  VALUE start,
  VALUE step,
  VALUE o0, VALUE o1, VALUE o_offset, VALUE o_stride,
  VALUE n)
{
  VdspArrayParam _i;
  VdspArrayParam _o;

  array_param2(&_i, i0, i1, i_offset, i_stride);
  double _start = NUM2DBL(start);
  double _step = NUM2DBL(step);
  array_param2(&_o, o0, o1, o_offset, o_stride);

  vDSP_Length _n = NUM2LONG(n);

  vDSP_vrampmuladd2D(
    _i.res0->v.d+_i.offset, _i.res1->v.d+_i.offset, _i.stride,
    &_start,
    &_step,
    _o.res0->v.d+_o.offset, _o.res1->v.d+_o.offset, _o.stride,
    _n
  );

  return rb_assoc_new(o0, o1);
}

VALUE rb_double_vgenp(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n,
  VALUE m)
{
  VdspArrayParam _a;
  VdspArrayParam _b;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  array_param(&_c, c, c_offset, c_stride);

  vDSP_Length _n = NUM2LONG(n);
  vDSP_Length _m = NUM2LONG(m);

  vDSP_vgenpD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _n,
    _m
  );

  return c;
}

VALUE rb_double_vtabi(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE s1,
  VALUE s2,
  VALUE c,
  VALUE m,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  double _s1 = NUM2DBL(s1);
  double _s2 = NUM2DBL(s2);
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);

  vDSP_Length _m = NUM2LONG(m);
  vDSP_Length _n = NUM2LONG(n);

  array_param(&_d, d, d_offset, d_stride);

  vDSP_vtabiD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_s1,
    &_s2,
    _c->v.d,
    _m,
    _d.res0->v.d+_d.offset, _d.stride,
    _n
  );

  return d;
}

VALUE rb_double_blkman_window(VALUE cls, VALUE c, VALUE n, VALUE flag)
{
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);
  vDSP_Length _n = NUM2LONG(n);
  int _flag = (int)NUM2LONG(flag);

  vDSP_blkman_windowD(_c->v.d, _n, _flag);

  return c;
}

VALUE rb_double_hamm_window(VALUE cls, VALUE c, VALUE n, VALUE flag)
{
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);
  vDSP_Length _n = NUM2LONG(n);
  int _flag = (int)NUM2LONG(flag);

  vDSP_hamm_windowD(_c->v.d, _n, _flag);

  return c;
}

VALUE rb_double_hann_window(VALUE cls, VALUE c, VALUE n, VALUE flag)
{
  VdspArrayNativeResource *_c = get_vdsp_array_native_resource(c);
  vDSP_Length _n = NUM2LONG(n);
  int _flag = (int)NUM2LONG(flag);

  vDSP_hann_windowD(_c->v.d, _n, _flag);

  return c;
}

VALUE rb_double_vclr(
  VALUE cls,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _c;

  array_param(&_c, c, c_offset, c_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vclrD(
    _c.res0->v.d+_c.offset, _c.stride,
    _n
  );

  return c;
}

VALUE rb_double_vfill(
  VALUE cls,
  VALUE a,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _c;

  double _a = NUM2DBL(a);
  array_param(&_c, c, c_offset, c_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vfillD(
    &_a,
    _c.res0->v.d+_c.offset, _c.stride,
    _n
  );

  return c;
}

VALUE rb_double_maxv(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_maxvD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_maxmgv(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_maxmgvD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_minv(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_minvD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_minmgv(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_minmgvD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_meanv(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_meanvD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_meamgv(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_meamgvD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_measqv(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_measqvD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_mvessq(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_mvessqD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_rmsqv(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_rmsqvD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_sve(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_sveD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_svemg(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_svemgD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_svesq(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_svesqD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_sve_svesq(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _sum;
  double _sum_of_squares;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_sve_svesqD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_sum,
    &_sum_of_squares,
    _n
  );

  return rb_assoc_new(DBL2NUM(_sum), DBL2NUM(_sum_of_squares));
}

VALUE rb_double_svs(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE n)
{
  VdspArrayParam _a;
  double _c;

  array_param(&_a, a, a_offset, a_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_svsD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_c,
    _n
  );

  return DBL2NUM(_c);
}

VALUE rb_double_vswap(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b, VALUE b_offset, VALUE b_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _b;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_b, b, b_offset, b_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vswapD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    _n
  );

  return rb_assoc_new(a, b);
}

VALUE rb_double_vtmerg(
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
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vtmergD(
    _a.res0->v.d+_a.offset, _a.stride,
    _b.res0->v.d+_b.offset, _b.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _n
  );

  return c;
}

VALUE rb_double_biquad(
  VALUE cls,
  VALUE biquad,
  VALUE x, VALUE x_offset, VALUE x_stride,
  VALUE y, VALUE y_offset, VALUE y_stride,
  VALUE n)
{
  VdspArrayParam _x;
  VdspArrayParam _y;

  VdspBiquadNativeResource *_b = get_vdsp_biquad_native_resource(biquad);
  array_param(&_x, x, x_offset, x_stride);
  array_param(&_y, y, y_offset, y_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_biquadD(
    _b->setup.d,
    _b->delay.d,
    _x.res0->v.d+_x.offset, _x.stride,
    _y.res0->v.d+_y.offset, _y.stride,
    _n);

  return y;
}

VALUE rb_double_vabs(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_c, c, c_offset, c_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vabsD(
    _a.res0->v.d+_a.offset, _a.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

VALUE rb_double_vnabs(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_c, c, c_offset, c_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vnabsD(
    _a.res0->v.d+_a.offset, _a.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

VALUE rb_double_vneg(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_c, c, c_offset, c_stride);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vnegD(
    _a.res0->v.d+_a.offset, _a.stride,
    _c.res0->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

VALUE rb_double_vclip(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_d, d, d_offset, d_stride);
  double _b = NUM2DBL(b);
  double _c = NUM2DBL(c);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vclipD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    &_c,
    _d.res0->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

VALUE rb_double_vclipc(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _d;
  vDSP_Length _n_low = 0;
  vDSP_Length _n_high = 0;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_d, d, d_offset, d_stride);
  double _b = NUM2DBL(b);
  double _c = NUM2DBL(c);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vclipcD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    &_c,
    _d.res0->v.d+_d.offset, _d.stride,
    _n,
    &_n_low,
    &_n_high);

  return rb_ary_new3(3, d, LONG2NUM(_n_low), LONG2NUM(_n_high));
}

VALUE rb_double_viclip(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_d, d, d_offset, d_stride);
  double _b = NUM2DBL(b);
  double _c = NUM2DBL(c);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_viclipD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    &_c,
    _d.res0->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

VALUE rb_double_vthr(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_c, c, c_offset, c_stride);
  double _b = NUM2DBL(b);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vthrD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    _c.res0->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

VALUE rb_double_vlim(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_d, d, d_offset, d_stride);
  double _b = NUM2DBL(b);
  double _c = NUM2DBL(c);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vlimD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    &_c,
    _d.res0->v.d+_d.offset, _d.stride,
    _n);

  return d;
}

VALUE rb_double_vthres(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c, VALUE c_offset, VALUE c_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _c;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_c, c, c_offset, c_stride);
  double _b = NUM2DBL(b);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vthresD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    _c.res0->v.d+_c.offset, _c.stride,
    _n);

  return c;
}

VALUE rb_double_vthrsc(
  VALUE cls,
  VALUE a, VALUE a_offset, VALUE a_stride,
  VALUE b,
  VALUE c,
  VALUE d, VALUE d_offset, VALUE d_stride,
  VALUE n)
{
  VdspArrayParam _a;
  VdspArrayParam _d;

  array_param(&_a, a, a_offset, a_stride);
  array_param(&_d, d, d_offset, d_stride);
  double _b = NUM2DBL(b);
  double _c = NUM2DBL(c);
  vDSP_Length _n = NUM2LONG(n);

  vDSP_vthrscD(
    _a.res0->v.d+_a.offset, _a.stride,
    &_b,
    &_c,
    _d.res0->v.d+_d.offset, _d.stride,
    _n);

  return d;
}


// Init

void Init_vdsp()
{
  // Vdsp
  rb_mVdsp = rb_define_module("Vdsp");
  rb_define_const(rb_mVdsp, "HALF_WINDOW", LONG2NUM(vDSP_HALF_WINDOW));
  rb_define_const(rb_mVdsp, "FULL_WINDOW", LONG2NUM(0));

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
  rb_include_module(rb_cDoubleArray, rb_mEnumerable);
  rb_define_method(rb_cDoubleArray, "initialize", rb_double_array_initialize, 1);
  rb_define_private_method(rb_cDoubleArray, "initialize_copy", rb_double_array_initialize_copy, 1);
  rb_define_singleton_method(rb_cDoubleArray, "create", rb_double_array_create, 1);
  rb_define_method(rb_cDoubleArray, "to_da", rb_double_array_to_da, 0);
  rb_define_method(rb_cDoubleArray, "+", rb_double_array_plus, 1);
  rb_define_method(rb_cDoubleArray, "-", rb_double_array_minus, 1);
  rb_define_method(rb_cDoubleArray, "*", rb_double_array_mul, 1);
  rb_define_method(rb_cDoubleArray, "/", rb_double_array_div, 1);
  rb_define_method(rb_cDoubleArray, "[]", rb_double_array_aref, -1);
  rb_define_method(rb_cDoubleArray, "[]=", rb_double_array_aset, -1);
  rb_define_method(rb_cDoubleArray, "each", rb_double_array_each, 0);
  rb_define_method(rb_cDoubleArray, "to_a", rb_double_array_get_values, 0);
  rb_define_method(rb_cDoubleArray, "coerce", rb_double_array_coerce, 1);
  rb_define_method(rb_cDoubleArray, "slice", rb_double_array_aref, -1);
  rb_define_method(rb_cDoubleArray, "resize", rb_double_array_resize, 1);
  rb_define_method(rb_cDoubleArray, "concat", rb_double_array_concat, -1);
  rb_define_method(rb_cDoubleArray, "abs", rb_double_array_abs, 0);
  rb_define_method(rb_cDoubleArray, "abs!", rb_double_array_abs_bang, 0);
  rb_define_method(rb_cDoubleArray, "nabs", rb_double_array_nabs, 0);
  rb_define_method(rb_cDoubleArray, "nabs!", rb_double_array_nabs_bang, 0);
  rb_define_method(rb_cDoubleArray, "negative", rb_double_array_negative, 0);
  rb_define_alias(rb_cDoubleArray,  "neg", "negative");
  rb_define_method(rb_cDoubleArray, "negative!", rb_double_array_negative_bang, 0);
  rb_define_alias(rb_cDoubleArray,  "neg!", "negative!");

  // Vdsp::DoubleArray Vector Generation
  rb_define_singleton_method(rb_cDoubleArray, "vramp", rb_double_array_vramp, 3);
  rb_define_singleton_method(rb_cDoubleArray, "vgen", rb_double_array_vgen, 3);
  rb_define_singleton_method(rb_cDoubleArray, "blkman_window", rb_double_array_blkman_window, -1);
  rb_define_singleton_method(rb_cDoubleArray, "hamm_window", rb_double_array_hamm_window, -1);
  rb_define_singleton_method(rb_cDoubleArray, "hann_window", rb_double_array_hann_window, -1);

  // Vdsp::DoubleArray Vector Clear and Fill Functions
  rb_define_method(rb_cDoubleArray, "clear!", rb_double_array_clear_bang, 0);
  rb_define_method(rb_cDoubleArray, "fill!", rb_double_array_fill_bang, 1);

  // Vdsp::DoubleArray Vector Extrema Calculation
  rb_define_method(rb_cDoubleArray, "maxv", rb_double_array_maxv, 0);
  rb_define_method(rb_cDoubleArray, "maxmgv", rb_double_array_maxmgv, 0);
  rb_define_method(rb_cDoubleArray, "minv", rb_double_array_minv, 0);
  rb_define_method(rb_cDoubleArray, "minmgv", rb_double_array_minmgv, 0);

  // Vdsp::DoubleArray Vector Average Calculation
  rb_define_method(rb_cDoubleArray, "meanv", rb_double_array_meanv, 0);
  rb_define_method(rb_cDoubleArray, "meamgv", rb_double_array_meamgv, 0);
  rb_define_method(rb_cDoubleArray, "measqv", rb_double_array_measqv, 0);
  rb_define_method(rb_cDoubleArray, "mvessq", rb_double_array_mvessq, 0);
  rb_define_method(rb_cDoubleArray, "rmsqv", rb_double_array_rmsqv, 0);

  // Vdsp::DoubleArray Vector Summation
  rb_define_method(rb_cDoubleArray, "sve", rb_double_array_sve, 0);
  rb_define_method(rb_cDoubleArray, "svemg", rb_double_array_svemg, 0);
  rb_define_method(rb_cDoubleArray, "svesq", rb_double_array_svesq, 0);
  rb_define_method(rb_cDoubleArray, "sve_svesq", rb_double_array_sve_svesq, 0);
  rb_define_method(rb_cDoubleArray, "svs", rb_double_array_svs, 0);

  // Vdsp::Biquad
  rb_mVdspBiquad = rb_define_module_under(rb_mVdsp, "Biquad");
  rb_define_method(rb_mVdspBiquad, "sections", rb_vdsp_biquad_sections, 0);
  rb_define_method(rb_mVdspBiquad, "alloc_sections", rb_vdsp_biquad_alloc_sections, 0);

  // Vdsp::DoubleBiquad
  rb_cDoubleBiquad = rb_define_class_under(rb_mVdsp, "DoubleBiquad", rb_cObject);
  rb_include_module(rb_cDoubleBiquad, rb_mVdspBiquad);
  rb_define_method(rb_cDoubleBiquad, "initialize", rb_double_biquad_initialize, 1);
  rb_define_singleton_method(rb_cDoubleBiquad, "create", rb_double_biquad_create, 1);
  rb_define_method(rb_cDoubleBiquad, "coefficients=", rb_double_biquad_set_coefficients, 1);
  rb_define_method(rb_cDoubleBiquad, "coefficients", rb_double_biquad_get_coefficients, 0);
  rb_define_method(rb_cDoubleBiquad, "apply", rb_double_biquad_apply, 1);

  // Vdsp::FFT
  rb_mVdspFFT = rb_define_module_under(rb_mVdsp, "FFT");
  rb_define_const(rb_mVdspFFT, "Radix2", LONG2NUM(kFFTRadix2));
  rb_define_const(rb_mVdspFFT, "Radix3", LONG2NUM(kFFTRadix3));
  rb_define_const(rb_mVdspFFT, "Radix5", LONG2NUM(kFFTRadix5));

  // Vdsp::DoubleFFT
  rb_cDoubleFFT = rb_define_class_under(rb_mVdsp, "DoubleFFT", rb_cObject);
  rb_include_module(rb_cDoubleFFT, rb_mVdspFFT);
  rb_define_method(rb_cDoubleFFT, "initialize", rb_double_fft_initialize, -1);
  rb_define_method(rb_cDoubleFFT, "forward", rb_double_fft_forward, -1);
  rb_define_method(rb_cDoubleFFT, "inverse", rb_double_fft_inverse, 2);
  rb_define_alias(rb_cDoubleFFT,  "forward_rect", "forward");
  rb_define_alias(rb_cDoubleFFT,  "inverse_rect", "inverse");
  rb_define_method(rb_cDoubleFFT, "forward_polar", rb_double_fft_forward_polar, -1);
  rb_define_method(rb_cDoubleFFT, "inverse_polar", rb_double_fft_inverse_polar, 2);

  // Vdsp::UnsafeDouble
  rb_mUnsafeDouble = rb_define_module_under(rb_mVdsp, "UnsafeDouble");
  rb_define_singleton_method(rb_mUnsafeDouble, "copy", rb_double_copy, 7);

  // Vdsp::UnsafeDouble Vector-based Arithmetic
  rb_define_singleton_method(rb_mUnsafeDouble, "vsadd", rb_double_vsadd, 8);
  rb_define_singleton_method(rb_mUnsafeDouble, "vadd", rb_double_vadd, 10);
  rb_define_singleton_method(rb_mUnsafeDouble, "vsub", rb_double_vsub, 10);
  rb_define_singleton_method(rb_mUnsafeDouble, "vsmul", rb_double_vsmul, 8);
  rb_define_singleton_method(rb_mUnsafeDouble, "vmul", rb_double_vmul, 10);
  rb_define_singleton_method(rb_mUnsafeDouble, "vsdiv", rb_double_vsdiv, 8);
  rb_define_singleton_method(rb_mUnsafeDouble, "svdiv", rb_double_svdiv, 8);
  rb_define_singleton_method(rb_mUnsafeDouble, "vaddsub", rb_double_vaddsub, 13);
  rb_define_singleton_method(rb_mUnsafeDouble, "vdiv", rb_double_vdiv, 10);
  rb_define_singleton_method(rb_mUnsafeDouble, "vasm", rb_double_vasm, 11);
  rb_define_singleton_method(rb_mUnsafeDouble, "vam", rb_double_vam, 13);
  rb_define_singleton_method(rb_mUnsafeDouble, "vsbsm", rb_double_vsbsm, 11);
  rb_define_singleton_method(rb_mUnsafeDouble, "vsbm", rb_double_vsbm, 13);
  rb_define_singleton_method(rb_mUnsafeDouble, "vmsa", rb_double_vmsa, 11);
  rb_define_singleton_method(rb_mUnsafeDouble, "vsma", rb_double_vsma, 11);
  rb_define_singleton_method(rb_mUnsafeDouble, "vma", rb_double_vma, 13);
  rb_define_singleton_method(rb_mUnsafeDouble, "vmsb", rb_double_vmsb, 13);
  rb_define_singleton_method(rb_mUnsafeDouble, "vsmsma", rb_double_vsmsma, 12);
  rb_define_singleton_method(rb_mUnsafeDouble, "vaam", rb_double_vaam, -1);
  rb_define_singleton_method(rb_mUnsafeDouble, "vmmsb", rb_double_vmmsb, -1);
  rb_define_singleton_method(rb_mUnsafeDouble, "vsbsbm", rb_double_vsbsbm, -1);
  rb_define_singleton_method(rb_mUnsafeDouble, "vasbm", rb_double_vasbm, -1);

  // Vdsp::UnsafeDouble Vector Generation
  rb_define_singleton_method(rb_mUnsafeDouble, "vramp", rb_double_vramp, 6);
  rb_define_singleton_method(rb_mUnsafeDouble, "vgen", rb_double_vgen, 6);
  rb_define_singleton_method(rb_mUnsafeDouble, "vrampmul", rb_double_vrampmul, 9);
  rb_define_singleton_method(rb_mUnsafeDouble, "vrampmul2", rb_double_vrampmul2, 11);
  rb_define_singleton_method(rb_mUnsafeDouble, "vrampmuladd", rb_double_vrampmuladd, 9);
  rb_define_singleton_method(rb_mUnsafeDouble, "vrampmuladd2", rb_double_vrampmuladd2, 11);
  rb_define_singleton_method(rb_mUnsafeDouble, "vgenp", rb_double_vgenp, 11);
  rb_define_singleton_method(rb_mUnsafeDouble, "vtabi", rb_double_vtabi, 11);
  rb_define_singleton_method(rb_mUnsafeDouble, "blkman_window", rb_double_blkman_window, 3);
  rb_define_singleton_method(rb_mUnsafeDouble, "hamm_window", rb_double_hamm_window, 3);
  rb_define_singleton_method(rb_mUnsafeDouble, "hann_window", rb_double_hann_window, 3);

  // Vdsp::UnsafeDouble Vector Clear and Fill Functions
  rb_define_singleton_method(rb_mUnsafeDouble, "vclr", rb_double_vclr, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "vfill", rb_double_vfill, 5);

  // Vdsp::UnsafeDouble Vector Extrema Calculation
  rb_define_singleton_method(rb_mUnsafeDouble, "maxv", rb_double_maxv, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "maxmgv", rb_double_maxmgv, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "minv", rb_double_minv, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "minmgv", rb_double_minmgv, 4);

  // Vdsp::UnsafeDouble Vector Average Calculation
  rb_define_singleton_method(rb_mUnsafeDouble, "meanv", rb_double_meanv, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "meamgv", rb_double_meamgv, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "measqv", rb_double_measqv, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "mvessq", rb_double_mvessq, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "rmsqv", rb_double_rmsqv, 4);

  // Vdsp::UnsafeDouble Vector Summation
  rb_define_singleton_method(rb_mUnsafeDouble, "sve", rb_double_sve, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "svemg", rb_double_svemg, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "svesq", rb_double_svesq, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "sve_svesq", rb_double_sve_svesq, 4);
  rb_define_singleton_method(rb_mUnsafeDouble, "svs", rb_double_svs, 4);

  // Vdsp::UnsafeDouble Copying, Element Swapping, and Merging Functions
  rb_define_singleton_method(rb_mUnsafeDouble, "vswap", rb_double_vswap, 7);
  rb_define_singleton_method(rb_mUnsafeDouble, "vtmerg", rb_double_vtmerg, 10);

  // Vdsp::UnsafeDouble Biquadratic IIR Filters
  rb_define_singleton_method(rb_mUnsafeDouble, "biquad", rb_double_biquad, 8);

  // Vdsp::UnsafeDouble Absolute and Negation Functions
  rb_define_singleton_method(rb_mUnsafeDouble, "vabs", rb_double_vabs, 7);
  rb_define_singleton_method(rb_mUnsafeDouble, "vnabs", rb_double_vnabs, 7);
  rb_define_singleton_method(rb_mUnsafeDouble, "vneg", rb_double_vneg, 7);

  // Vdsp::UnsafeDouble Clipping, Limit, and Threshold Operations
  rb_define_singleton_method(rb_mUnsafeDouble, "vclip", rb_double_vclip, 9);
  rb_define_singleton_method(rb_mUnsafeDouble, "vclipc", rb_double_vclipc, 9);
  rb_define_singleton_method(rb_mUnsafeDouble, "viclip", rb_double_viclip, 9);
  rb_define_singleton_method(rb_mUnsafeDouble, "vthr", rb_double_vthr, 8);
  rb_define_singleton_method(rb_mUnsafeDouble, "vlim", rb_double_vlim, 9);
  rb_define_singleton_method(rb_mUnsafeDouble, "vthres", rb_double_vthres, 8);
  rb_define_singleton_method(rb_mUnsafeDouble, "vthrsc", rb_double_vthrsc, 9);
}
