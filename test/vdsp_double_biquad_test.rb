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
    biquad = Vdsp::DoubleBiquad.new(coef)

    assert_instance_of(Vdsp::Biquad::Coefficient, coef)
    assert_instance_of(Vdsp::DoubleBiquad, biquad)
  end

  def test_sections
    coef = Vdsp::Biquad::Coefficient.new(@b0/@a0, @b1/@a0, @b2/@a0, @a1/@a0, @a2/@a0)

    a = Vdsp::DoubleBiquad.new(coef)
    b = Vdsp::DoubleBiquad.new([coef])
    c = Vdsp::DoubleBiquad.new([coef, coef])

    assert_equal 1, a.sections
    assert_equal 1, b.sections
    assert_equal 2, c.sections
  end

  def test_coefficients
    coef = Vdsp::Biquad::Coefficient.new(@b0/@a0, @b1/@a0, @b2/@a0, @a1/@a0, @a2/@a0)
    biquad = Vdsp::DoubleBiquad.new(coef)
    coefs = biquad.coefficients

    assert_equal 1, coefs.size
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
