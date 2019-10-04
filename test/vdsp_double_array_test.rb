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

  def test_clone
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = a.clone
    c = a.dup

    assert a!=b
    assert a!=c
    assert b!=c

    assert_equal a.to_a, b.to_a
    assert_equal a.to_a, c.to_a
    assert_equal a.length, b.length
    assert_equal a.length, c.length

    a[0] = 4.0
    assert a[0]!=b[0]
    assert a[0]!=c[0]
  end

  def test_length
    a = Vdsp::DoubleArray.new(3)
    assert_equal 3, a.length
    assert_equal 3, a.size
  end

  def test_each
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])

    yield_values = []
    a.each {|v|
      yield_values << v
    }
    assert_equal a.to_a, yield_values
    assert_equal a.to_a, a.each.to_a
  end

  def test_enumerable
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    b = a.map{|v| v * 2}

    assert_equal [2.0, 4.0, 6.0], b
  end

  def test_to_a
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    assert_equal [1.0, 2.0, 3.0], a.to_a
  end

  def test_aref
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    assert_equal 1.0, a[0]
    assert_equal 2.0, a[1]
    assert_equal 3.0, a[2]
  end

  def test_aset
    a = Vdsp::DoubleArray.new(3)
    a[0] = 1.0
    a[1] = 2.0
    a[2] = 3.0

    assert_equal 1.0, a[0]
    assert_equal 2.0, a[1]
    assert_equal 3.0, a[2]
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

  def test_vramp
    # TODO
  end

  def test_vgen
    # TODO
  end

  def test_blkman_window
    flag = Vdsp::FULL_WINDOW
    a = Vdsp::DoubleArray.blkman_window(4, flag)

    assert a[0]<a[1]
    assert a[1]<a[2]
    assert a[2]>a[3]
  end

  def test_hamm_window
    flag = Vdsp::FULL_WINDOW
    a = Vdsp::DoubleArray.hamm_window(4, flag)

    assert a[0]<a[1]
    assert a[1]<a[2]
    assert a[2]>a[3]
  end

  def test_hann_window
    flag = Vdsp::FULL_WINDOW
    a = Vdsp::DoubleArray.hann_window(4, flag)

    assert a[0]<a[1]
    assert a[1]<a[2]
    assert a[2]>a[3]
  end

  def test_clear_bang
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a.clear!

    assert_equal 0.0, a[0]
    assert_equal 0.0, a[1]
    assert_equal 0.0, a[2]
  end

  def test_fill_bang
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    a.fill!(4.0)

    assert_equal 4.0, a[0]
    assert_equal 4.0, a[1]
    assert_equal 4.0, a[2]
  end

  def test_maxv
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.maxv

    assert_equal 3.0, c
  end

  def test_maxmgv
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.maxmgv

    assert_equal 3.0, c
  end

  def test_minv
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.minv

    assert_equal 1.0, c
  end

  def test_minmgv
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.minmgv

    assert_equal 1.0, c
  end

  def test_meanv
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.meanv

    assert_equal 2.0, c
  end

  def test_meamgv
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.meamgv

    assert_equal 2.0, c
  end

  def test_measqv
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.measqv

    assert_in_delta 4.666666666666667, c
  end

  def test_mvessq
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.mvessq

    assert_in_delta 4.666666666666667, c
  end

  def test_rmsqv
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.rmsqv

    assert_in_delta 2.160246899469287, c
  end

  def test_sve
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.sve

    assert_equal 6.0, c
  end

  def test_svemg
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.svemg

    assert_equal 6.0, c
  end

  def test_svesq
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.svesq

    assert_equal 14.0, c
  end

  def test_sve_svesq
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    sum, sum_of_squares = a.sve_svesq

    assert_equal 6.0, sum
    assert_equal 14.0, sum_of_squares
  end

  def test_svs
    a = Vdsp::DoubleArray.create([1.0, 2.0, 3.0])
    c = a.svs

    assert_equal 14.0, c
  end

  def test_biquad
    # TODO
  end

  def test_abs
    a = Vdsp::DoubleArray.create([-1.0, -2.0, 3.0])
    a = a.abs

    assert_equal 1.0, a[0]
    assert_equal 2.0, a[1]
    assert_equal 3.0, a[2]
  end

  def test_abs_bang
    a = Vdsp::DoubleArray.create([-1.0, -2.0, 3.0])
    a.abs!

    assert_equal 1.0, a[0]
    assert_equal 2.0, a[1]
    assert_equal 3.0, a[2]
  end

  def test_nabs
    a = Vdsp::DoubleArray.create([-1.0, -2.0, 3.0])
    a = a.nabs

    assert_equal -1.0, a[0]
    assert_equal -2.0, a[1]
    assert_equal -3.0, a[2]
  end

  def test_nabs_bang
    a = Vdsp::DoubleArray.create([-1.0, -2.0, 3.0])
    a.nabs!

    assert_equal -1.0, a[0]
    assert_equal -2.0, a[1]
    assert_equal -3.0, a[2]
  end

  def test_negative
    a = Vdsp::DoubleArray.create([-1.0, -2.0, 3.0])
    a = a.negative

    assert_equal 1.0, a[0]
    assert_equal 2.0, a[1]
    assert_equal -3.0, a[2]
  end

  def test_negative_bang
    a = Vdsp::DoubleArray.create([-1.0, -2.0, 3.0])
    a.negative!

    assert_equal 1.0, a[0]
    assert_equal 2.0, a[1]
    assert_equal -3.0, a[2]
  end
end
