require "test_helper"

class VdspDoubleFFTTest < Minitest::Test
  def setup
    @freq = 440.0
    @samplerate = 44100
    @len = 1024

    sinwave = @len.times.map {|i|
      Math.sin(@freq / @samplerate * 2 * Math::PI * i)
    }
    @real = Vdsp::DoubleArray.create(sinwave).hann_window
    @imag = Vdsp::DoubleArray.new(@len)
  end

  def test_forward
    # TODO
  end

  def test_inverse
    fft = Vdsp::DoubleFFT.new(@len)

    # real
    r_real, r_imag = fft.forward(@real)
    r_real, r_imag = fft.inverse(r_real, r_imag)

    assert_in_delta @real[10], r_real[10]
    assert_in_delta @real[20], r_real[20]
    assert_in_delta @real[30], r_real[30]

    assert_nil r_imag

    # complex
    z_real, z_imag = fft.forward(@real, @imag)
    z_real, z_imag = fft.inverse(z_real, z_imag)

    assert_in_delta @real[10], z_real[10]
    assert_in_delta @real[20], z_real[20]
    assert_in_delta @real[30], z_real[30]

    assert_in_delta 0.0, z_imag[10]
    assert_in_delta 0.0, z_imag[20]
    assert_in_delta 0.0, z_imag[30]
  end

  def test_forward_polar
    # TODO
  end

  def test_inverse_polar
    fft = Vdsp::DoubleFFT.new(@len)

    # real
    r_magnitude, r_phase = fft.forward_polar(@real)
    r_real, r_imag = fft.inverse_polar(r_magnitude, r_phase)

    assert_in_delta @real[10], r_real[10]
    assert_in_delta @real[20], r_real[20]
    assert_in_delta @real[30], r_real[30]

    assert_nil r_imag

    # complex
    z_magnitude, z_phase = fft.forward_polar(@real, @imag)
    z_real, z_imag = fft.inverse_polar(z_magnitude, z_phase)

    assert_in_delta @real[10], z_real[10]
    assert_in_delta @real[20], z_real[20]
    assert_in_delta @real[30], z_real[30]

    assert_in_delta 0.0, z_imag[10]
    assert_in_delta 0.0, z_imag[20]
    assert_in_delta 0.0, z_imag[30]
  end
end
