require "test_helper"

class VdspDoubleArrayTest < Minitest::Test
  def test_create
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = Vdsp::DoubleArray.new(3)
    c = a.to_da

    assert_instance_of(Vdsp::DoubleArray, a)
    assert_instance_of(Vdsp::DoubleArray, b)
    assert_instance_of(Vdsp::DoubleArray, c)
  end

  def test_length
    a = Vdsp::DoubleArray.new(3)
    assert_equal 3, a.length
    assert_equal 3, a.size
  end

  def test_to_a
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    assert_equal [1.0, 2.0, 3.0], a.to_a
  end

  def test_operator_add
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])

    assert_equal [2.0, 4.0, 6.0], (a + a).to_a
    assert_equal [3.0, 4.0, 5.0], (a + 2.0).to_a
    assert_equal [3.0, 4.0, 5.0], (2.0 + a).to_a
  end

  def test_operator_sub
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])

    assert_equal [0.0, 0.0, 0.0], (a - a).to_a
    assert_equal [-1.0, 0.0, 1.0], (a - 2.0).to_a
    assert_equal [1.0, 0.0, -1.0], (2.0 - a).to_a
  end

  def test_operator_mul
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])

    assert_equal [1.0, 4.0, 9.0], (a * a).to_a
    assert_equal [2.0, 4.0, 6.0], (a * 2.0).to_a
    assert_equal [2.0, 4.0, 6.0], (2.0 * a).to_a
  end

  def test_operator_div
    a = Vdsp::DoubleArray.create([1.0, 2.0, 4.0])

    assert_equal [1.0, 1.0, 1.0], (a / a).to_a
    assert_equal [0.5, 1.0, 2.0], (a / 2.0).to_a
    assert_equal [2.0, 1.0, 0.5], (2.0 / a).to_a
  end
end
