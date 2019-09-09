require "test_helper"

class VdspTest < Minitest::Test
  def test_that_it_has_a_version_number
    refute_nil ::Vdsp::VERSION
  end

  def test_vsadd
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = 5.0
    c = Vdsp::DoubleArray.new(3)
    Vdsp::vsaddD(a, 1, b, c, 1, a.length)

    assert_equal 6.0, c[0]
    assert_equal 7.0, c[1]
    assert_equal 8.0, c[2]
  end
end
