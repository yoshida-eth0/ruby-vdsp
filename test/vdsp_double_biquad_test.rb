require "test_helper"

class VdspDoubleBiquadTest < Minitest::Test
  def setup
    freq = 2000.0
    q = 1.0 / Math.sqrt(2.0)
    samplerate = 44100.0

    omega = 2.0 * Math::PI * freq / samplerate
    alpha = Math.sin(omega) / (2.0 * q)

    @a0 = 1.0 + alpha
    @a1 = -2.0 * Math.cos(omega)
    @a2 = 1.0 - alpha
    @b0 = (1.0 - Math.cos(omega)) / 2.0 
    @b1 = 1.0 - Math.cos(omega)
    @b2 = (1.0 - Math.cos(omega)) / 2.0 
  end

  def test_create
    coef = Vdsp::Biquad::Coefficient.new(@b0/@a0, @b1/@a0, @b2/@a0, @a1/@a0, @a2/@a0)
    biquad1 = Vdsp::DoubleBiquad.new(1)
    biquad2 = Vdsp::DoubleBiquad.create(coef)

    assert_instance_of(Vdsp::Biquad::Coefficient, coef)
    assert_instance_of(Vdsp::DoubleBiquad, biquad1)
    assert_instance_of(Vdsp::DoubleBiquad, biquad2)
  end

  def test_sections
    coef = Vdsp::Biquad::Coefficient.new(@b0/@a0, @b1/@a0, @b2/@a0, @a1/@a0, @a2/@a0)

    a = Vdsp::DoubleBiquad.new(4)
    b = Vdsp::DoubleBiquad.create(coef)
    c = Vdsp::DoubleBiquad.create([coef])
    d = Vdsp::DoubleBiquad.create([coef, coef])

    assert_equal 0, a.sections
    assert_equal 1, b.sections
    assert_equal 1, c.sections
    assert_equal 2, d.sections
  end

  def test_alloc_sections
    coef = Vdsp::Biquad::Coefficient.new(@b0/@a0, @b1/@a0, @b2/@a0, @a1/@a0, @a2/@a0)

    a = Vdsp::DoubleBiquad.new(4)
    b = Vdsp::DoubleBiquad.create(coef)
    c = Vdsp::DoubleBiquad.create([coef])
    d = Vdsp::DoubleBiquad.create([coef, coef])

    assert_equal 4, a.alloc_sections
    assert_equal 1, b.alloc_sections
    assert_equal 1, c.alloc_sections
    assert_equal 2, d.alloc_sections
  end

  def test_get_coefficients
    coef = Vdsp::Biquad::Coefficient.new(@b0/@a0, @b1/@a0, @b2/@a0, @a1/@a0, @a2/@a0)
    biquad = Vdsp::DoubleBiquad.create(coef)
    coefs = biquad.coefficients

    assert_equal 1, coefs.size
    assert_in_delta coef.b0, coefs[0].b0
    assert_in_delta coef.b1, coefs[0].b1
    assert_in_delta coef.b2, coefs[0].b2
    assert_in_delta coef.a1, coefs[0].a1
    assert_in_delta coef.a2, coefs[0].a2
  end

  def test_set_coefficients
    coef = Vdsp::Biquad::Coefficient.new(@b0/@a0, @b1/@a0, @b2/@a0, @a1/@a0, @a2/@a0)
    biquad = Vdsp::DoubleBiquad.new(1)
    biquad.coefficients = coef
    coefs = biquad.coefficients

    assert_equal 1, biquad.sections
    assert_in_delta coef.b0, coefs[0].b0
    assert_in_delta coef.b1, coefs[0].b1
    assert_in_delta coef.b2, coefs[0].b2
    assert_in_delta coef.a1, coefs[0].a1
    assert_in_delta coef.a2, coefs[0].a2
  end

  def test_apply
    # TODO
  end
end
